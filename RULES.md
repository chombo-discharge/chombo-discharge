# Nearest-Neighbor Particle Merge — Algorithm Rules

Design specification for a distributed, MPI-safe, nearest-neighbor particle-merging algorithm.
Worked out through iterative design and validated with a Python/mpi4py prototype (6 domain
topologies in 2D and 3D, 10 particle configurations each, 6 adversarially-randomized local
processing orders each, up to 8 consecutive rounds each — 360 runs, all passing weight
conservation, particle-ID uniqueness, a per-particle lineage ledger audit, and determinism under
randomized processing order). This document is the rule set that implementation should follow;
it is not itself an implementation.

## Data model

- **Mandatory fields** the core matching/decision engine touches: position, weight, a globally
  unique ID, and owning rank. That's it — the engine is payload-agnostic.
- **Type-specific reconciliation data** (e.g. `energy` for `ItoParticle`) rides along as an
  **opaque payload** — never inspected by the matching engine, only unpacked by a caller-supplied
  combine/reconcile callback when a pair is finalized. Not a static per-type column declaration
  (deliberately not reusing HDF5's `h5PartColumns` mechanism for this — different merge
  algorithms legitimately want different fields, which only a caller-supplied lambda can express,
  mirroring the existing `BinaryParticleReconcile<P>` pattern). That callback should be a genuine
  template parameter, not `std::function`, since it sits on a hot path.
- **Owning rank is always read from the particle's own data** — single source of truth, never
  inferred from MPI transport/source-rank metadata.

## Sequencing (hard invariant)

- **Ghost-fill happens exactly once per round, strictly before any merge decision is made** —
  never reused/cached from earlier in a timestep, never done "after" some local pass. This is the
  fix for an early "compression toward patch interior" bug: every particle's nearest-neighbor
  search must always run against the *full* (local + ghost) view from the start. Restricting the
  search for some particles first (e.g. same-patch only) reintroduces exactly that bias, because
  boundary particles would get locked into worse matches before ever being allowed to look across
  a boundary.
- One ghost-fill serves the whole round — trivial-tier resolution (below) doesn't need ghosts at
  all, so no second fill is required.

## Crowding trigger

- Per-cell live count uses an **unclamped**, position-derived cell key (`floor(pos/dx)`), not an
  index clamped to any one patch's box — so a ghost particle buckets into its true origin cell,
  not a boundary cell of the receiving patch.
- A candidate pair is only eligible to merge if **both participants' cells** currently hold more
  than N particles — checked dynamically at the moment a pair is actually processed, not from a
  stale snapshot taken at candidate-generation time (a cell can drop below N mid-pass as earlier
  pairs in the same round consume particles).
- This rule is deliberately conservative: a cell can legitimately have >N particles whose true
  nearest neighbors are all in non-crowded cells, and the rule correctly refuses to merge them
  with an uncrowded partner. In a static test (no particle motion between rounds) this can show up
  as a permanent, *correct* standoff — not a bug, just something a real timestep's particle drift
  would very likely break.

## Boundary-exposed classification

- A particle is "boundary-exposed" iff it appears in **any** of its patch's own outgoing
  ghost-shipment targets — i.e. it would actually be shipped to at least one neighbor. Read
  directly off the existing ghost-mask machinery; no new geometry needed.
- This is the criterion that decides fast-path eligibility (below) — **not** "same rank" and
  **not** "same `ParticleSoA`." Both were tried and found unsound: a particle can be exposed to
  (visible from) a different patch owned by the very same rank, so "same rank" doesn't imply safe;
  a particle can be exposed while its partner happens to sit in the same `ParticleSoA`, so "same
  patch" doesn't either.

## Trivial tier (immediate, zero-coordination resolution)

- Eligible **iff neither participant is boundary-exposed** — regardless of which patch either one
  lives in.
- All eligible candidate edges (across every patch a rank owns) are processed in one **globally
  distance-sorted order**, with a deterministic secondary tiebreak on participant IDs — required
  because exact distance ties (e.g. true mutual nearest neighbors) are common enough that leaving
  them to arbitrary iteration/insertion order produces real, observed nondeterminism that cascades
  into unrelated downstream decisions.
- At the moment a pair is popped for processing (not at generation time): re-check both
  participants are still alive (lazy invalidation — skip if either was already consumed by a
  closer pair), re-check the local partner hasn't **in the meantime** committed to its own
  outgoing (non-trivial) proposal — a purely geometric "exposed" check alone is insufficient,
  since a non-exposed particle's true nearest neighbor can still be a foreign ghost even though
  the particle itself is unreachable externally; missing this dynamic check was a real
  double-consumption bug — and re-check both cells still exceed the crowding threshold.
- On commit: remove both participants, decrement their cells' live counts, record the merge.
  Left-unmatched candidates are not retried within the same round.

## Propose (for anything the trivial tier can't safely resolve)

- Every remaining over-threshold, unconsumed particle whose true best candidate involves an
  exposed participant emits **exactly one** outgoing proposal to whoever it found as its nearest
  neighbor.
- **Routing rule: always send to the target's owner, full stop — never by any ID comparison.** (An
  earlier "route to owner of the lower-ID particle" framing was found to be subtly wrong: it can
  send the two directions of a genuinely mutual pair to two different, uncoordinated judges.)
- If the target's owner is the same rank (even a different patch), this is a cheap in-process
  handoff, not an MPI message — but it still goes through the *exact same* pooled judge logic
  below, never a shortcut.
- The proposer is marked as having its own outgoing commitment (needed by the judge's
  reject-overlap rule).

## Judge

Per rank, pooling everything it owns is targeted by — same-rank in-process handoffs and genuinely
cross-rank MPI-received proposals, uniformly, in one accounting pass:

- Group incoming proposals by target.
- **A target with no outgoing commitment of its own**: eligible winner = closest incoming
  proposal (distance, with ID tiebreak).
- **A target that itself has an outgoing commitment** (i.e. it's simultaneously proposing to
  something else): normally reject *every* incoming proposal unconditionally — it's "spoken for."
  **Exception: if its own outgoing target is itself one of its current incoming proposers** — a
  genuinely mutual best-match pair, where both sides' independent, full-information searches agree
  — accept it. Blindly rejecting this case (the original, simpler rule) causes a **permanent
  deadlock**: a stable mutual pair proposes to each other identically every round and gets
  rejected every round forever, since nothing changes to break the symmetry. Both sides detect the
  same mutual match independently, so a deterministic tiebreak (**lower-ID side performs the
  accept**, higher-ID side falls back to the normal reject) is required — otherwise both sides
  would accept and double-merge the same two particles.
- **All eligible target decisions (plain-argmin and mutual-match alike) are resolved in one
  globally distance-sorted pass (with an ID tiebreak), not by iterating targets in arbitrary
  order.** Accepting one target decrements its cell's live count, which can push a *different*
  target sharing that cell below threshold — processing targets in whatever order a shuffled
  proposal list happened to produce made the outcome order-dependent. Sorting closest-first with a
  dynamic per-target threshold recheck (mirroring the trivial tier's discipline exactly) fixes
  this, and also closes a related gap where the mutual-match path had skipped that dynamic recheck
  entirely.
- **No N-way corner limitation**: this handles any number of simultaneous suitors for one
  particle, not just pairwise boundaries — every proposal naming a given target funnels to the
  same judge regardless of how many neighbors it comes from, so a corner where 3+ patches meet is
  resolved by the exact same logic as an ordinary two-patch boundary, with no special-casing.
  Verified concretely by the 8-rank 3D corner and mixed-ownership scenarios in the prototype sweep.

## Verdict and placement

- Judge replies accept/reject to every proposer (in-process or MPI, matching how the proposal
  arrived). Accept ⇒ judge creates the merged particle (weighted-centroid position,
  weighted-average opaque payload, summed weight, fresh globally-unique ID) and removes its
  target; the winning proposer's owner removes its own source particle on receiving "accepted."
  Reject ⇒ no action, reconsidered next round.
- **Placement applies uniformly to every newly created particle, trivial-tier or judge-decided
  alike**: a local point-location check against the (globally replicated) domain decomposition;
  append directly if it lands in one of the creating rank's own patches, otherwise scatter it to
  whichever rank actually owns that region. Needed even for trivial-tier (same-rank) merges, since
  the union of two patches isn't generally convex (T-junctions, corners, refinement boundaries).

## Single pass, caller loops for convergence

The whole pipeline runs once per round; it never internally retries a rejected or unmatched
particle within the same round — that's deferred to the next round's fresh ghost-fill and fresh
candidate generation. Deliberate, accepted conservatism: not globally optimal in one pass, but
nothing incorrect happens.

### Optional: iterate the local (trivial) tier to convergence within a round

A caller-selectable option (`a_iterateLocalTierToConvergence`, default `false`) partially relaxes
this for the LOCAL portion only — the cross-patch propose/judge/verdict exchange is unaffected,
still exactly once per call, always.

- **The gap this closes**: `findNearestNeighborCandidates()` fixes each query particle's
  nearest-neighbor edge once, at the start of a pass. If that candidate gets consumed by a
  different, closer edge before this one is reached in the sorted processing order, the query
  particle is left unmatched for the *rest* of that pass — even if a different, still-alive
  partner exists nearby. There is no second query within the same pass.
- **Measured severity (Python prototype)**: a single isolated, non-exposed cell seeded at 64
  particles against a threshold of 16 — i.e. trivial tier alone has everything it needs to fully
  drain the cell, no cross-patch coordination required at all — completed only ~12 merges in one
  un-iterated pass. That's far below the 24 merges needed to reach 16, and even below the
  **32-merge hard ceiling a single non-iterated pass can ever reach regardless of match quality**
  (each participant gets exactly one query per pass, so one pass can never do better than n/2
  merges — going from 64 all the way to 16 needs at least two full halvings, i.e. at least two
  *something* — either two external timesteps, or, with this option, two internal local
  iterations within the same call).
- **The fix**: when enabled, the candidate-search-plus-trivial-tier-resolve step repeats — reusing
  the SAME ghost-fill and the SAME per-particle boundary-exposed flags (neither depends on which
  local particles are still alive, so neither needs recomputing) — until one full local iteration
  commits zero further merges. This is cheap (no new ghost-fill, no new MPI) precisely because
  nothing about local re-matching requires new data, only a fresh look at who's still alive.
- **Measured effect turned out much smaller than the reasoning above predicts**: enabling this
  option in the same 64-vs-16 test recovered only ~5% more merges, converging in 3-4 internal
  iterations (the vast majority of the gain happens in the very first iteration; each subsequent
  one contributes only a handful more before hitting zero). Re-running the search only helps a
  particle whose actual BEST choice changes between iterations, and that turns out to be rare —
  most particles left unmatched after the first pass keep re-discovering the exact same
  unavailable candidate every time, because nothing about that candidate's situation changes
  either. The dominant loss mechanism is a different, *stable* condition that iteration cannot
  fix at all — see the next section, which is what was actually needed to close most of the gap.
- **Tradeoff**: more local re-search compute per call, in exchange for a real but modest reduction
  in how many external timesteps are needed to clear a severe local backlog. Left off by default
  so existing documented behavior doesn't change silently. Still worth having — it composes with
  the fallback-candidate option below rather than competing with it — just not the primary lever.

### Optional: retry with a fallback candidate when the first choice is blocked

A second, independent caller-selectable option (`a_maxFallbackCandidates`, default `0`, forwarded
to `resolveTrivialTier()`) turned out to be the one that actually closes most of the single-round
completion gap.

- **What it does**: when a query particle's chosen candidate is blocked — already consumed, OR
  already has its own outgoing commitment, OR exposed/a ghost — and the query particle itself is
  not exposed (if it is, no local candidate could ever help, so don't bother retrying):
  immediately re-query for the next-nearest still-available candidate, excluding every candidate
  already tried for this edge, up to this many extra attempts, before giving up. Resolved
  immediately at the point of failure rather than deferred to preserve strict global
  distance-sort order across every particle's retries — a deliberate simplification, justified
  because a retry only ever fires for an edge whose first choice already failed, so it was never
  competing on equal footing with anyone's first attempt anyway. A particle that exhausts every
  fallback still proposes (if it proposes at all) using its TRUE, original nearest neighbor —
  fallbacks only affect local trivial-tier eligibility, never what gets proposed cross-patch.
- **Why this is the dominant lever, not iteration**: instrumenting exactly why edges failed to
  commit in the same 64-vs-16 test found "candidate already has its own outgoing commitment" to
  be the single largest loss category — larger than genuine boundary exposure. A particle whose
  own best candidate happens to be exposed or a ghost gets marked "busy"; any OTHER particle whose
  nearest neighbor is that same busy particle then loses its edge too, even though the busy
  particle was never actually consumed — it's just reserved. This is a *stable* condition (the
  busy particle's own true nearest neighbor doesn't change), which is exactly why iterating the
  whole search doesn't fix it — only giving the BLOCKED particle a fallback candidate does.
- **A real bug found and fixed along the way**: implementing and measuring this surfaced a
  genuine discrepancy between the design as documented (`resolveTrivialTier()`: a query whose
  candidate turns out to already have its own outgoing commitment should be marked as having an
  outgoing commitment itself, so it gets a chance at `judgeProposals()`'s existing mutual-match
  rule) and what the prototype actually did (silently dropped it — no trivial merge, no proposal,
  no outcome at all for the round). This bug predates this option; it was present in the original
  `run_scenario.py` too, just never surfaced because nothing had previously exercised this path
  hard enough to make its effect visible. Fixing it in isolation (still `max_fallback_candidates`
  = 0) initially made completion measurably *worse* (`n_trivial` dropped from the buggy 7258 to
  7176) — correctly marking more particles "busy" makes the blocking cascade above fire *more*
  often, not less. The fallback option is what recovers from that regression and then some.
  Genuinely-stale candidates (already fully consumed, not just busy) are deliberately excluded
  from this fix — there's no meaningful id to propose to there, and marking the query particle
  "busy" in that case would incorrectly block others from claiming it, since it hasn't actually
  committed to anything.
- **Measured result** (same 64-vs-16 test, corrected baseline):

  | `max_fallback_candidates` | trivial merges | final count |
  |---|---|---|
  | 0 | 7176 | 26248 |
  | 1 | 7672 | 25838 |
  | 2 | 7771 | 25789 |
  | 3 | 7824 | 25759 |
  | 5 | 7871 | 25735 |
  | 10 | 7914 | 25707 |

  A single fallback candidate recovers most of the achievable gain (+496 merges, +6.9%); pushing
  from 1 to 10 adds only another +242 (+3.2%). `max_fallback_candidates = 1` or `2` is the
  practical choice — an open-ended search buys very little beyond that.
- **Tradeoff**: a handful of extra nearest-neighbor queries per blocked particle (bounded by the
  cap), in exchange for a real, substantial improvement in single-round completion — cheaper and
  more effective than `a_iterateLocalTierToConvergence` for the specific problem that motivated
  both options. The two compose freely and don't interfere with each other.

### Optional: a physical distance cap on merge eligibility

A third, independent caller-selectable option (`a_maxCellDistance`, default `std::nullopt` =
disabled) addresses a different concern entirely from the two above: nothing before this point
puts any limit on how physically far apart two merged particles can be — the search always finds
the TRUE nearest neighbor across the whole local+ghost neighborhood, however far that turns out to
be, and merges it if eligible. For a PIC-style code this is a real physical concern: merging two
particles separated by several grid cells represents a large, potentially non-physical
displacement of charge/mass, not the local coarse-graining a merge is meant to be.

- **What it does**: caps merge eligibility at a whole-number-of-cells Chebyshev (max-per-axis)
  distance between the query's and candidate's cell — same cell = 0, any of the 8 (2D) / 26 (3D)
  Moore-neighborhood-adjacent cells (including diagonals) = 1, and so on. Chebyshev distance on
  cell *indices* was chosen deliberately over a raw Euclidean position cutoff: "closest grid cell
  is OK, two cells is not" is a statement about which cells are adjacent, and a Euclidean cutoff
  would inconsistently include or exclude diagonal neighbors depending on exactly where within
  each cell the two particles happen to sit. The comparison uses the SAME unclamped,
  position-derived cell key already used for the crowding trigger (`floor(pos/dx)`) — not either
  particle's own patch-local cell index — which is what makes the comparison meaningful at all
  across a patch boundary.
- **Where it's enforced**: a candidate beyond the cap is treated as ineligible at the exact same
  point in `resolveTrivialTier()` as a stale or busy candidate — before trivial-tier eligibility
  is even checked. Critically, unlike the busy/exposed/ghost cases, a too-far candidate does
  **not** trigger a proposal either: the pair is exactly as physically invalid via
  `judgeProposals()` as it would be trivially, so there's nothing valid to hand a judge. If
  `a_maxFallbackCandidates` is enabled, the query gets a fallback retry for a closer candidate
  instead; otherwise it's simply left unmatched this round, with no side effects on anyone else's
  processing (same treatment as a stale candidate).
- **Verified directly** (Python prototype) with a hand-placed pair 2 cells apart and nothing else
  nearby, so they're unambiguously each other's only candidate: uncapped, they merge trivially at
  the exact weighted centroid. Capped at 1 cell, both survive completely untouched — no merge, no
  proposal, positions unchanged. Capped at exactly 2 cells (the boundary itself), they merge
  again — confirming the cap is inclusive (`<=`, not `<`).
- **Composes freely** with both options above — independent concern, no interaction either was
  found to have with it.

## What's proven vs. what's empirically validated

- **Mass conservation, no double-consumption, order-independence**: verified two ways — exact
  global-sum weight conservation, and (stronger) a per-particle lineage ledger confirming every
  original particle traces to exactly one final particle with matching weight, across 360 runs (6
  topologies × 10 configs × 6 shuffle-orderings × up to 8 rounds), all under adversarially
  randomized local processing order.
- **The mutual-pair fix specifically resolves 2-party reciprocal deadlocks** — confirmed both by
  direct reproduction and by the full sweep passing afterward.
- **Longer proposal cycles (3+ particles, none pairwise mutual)**: there is a rigorous argument for
  why they cannot occur — under full mutual visibility, a nearest-neighbor functional graph cannot
  contain a cycle of length ≥3 (the max-weight edge in any such cycle immediately contradicts the
  nearest-neighbor property of one of its endpoints), so only mutual pairs or non-cyclic chains are
  possible, which the rules above already cover. Checked both mathematically and by a deliberate
  adversarial attempt to break it via asymmetric ghost visibility (mismatched patch geometry
  designed to make A see B, B see C, C see A non-reciprocally) — no deadlock found. This should be
  read as strong theoretical plus empirical support, not an exhaustive proof covering every
  possible asymmetric-visibility configuration in the real system — worth staying alert to during
  real integration, not something to treat as fully closed.

## Spatial aliasing — empirical findings (separate investigation, Python prototype)

Beyond the correctness properties above, a follow-up investigation checked whether repeated
merging distorts the *spatial* distribution of surviving particles — the classic failure mode
here is a cell-based merge scheme that snaps merged particles toward cell centers, stamping a
periodic comb onto the density field (visible as spurious spectral peaks at the cell frequency).
Method: repeatedly (1) inject fresh uniformly-sampled particles, (2) merge (old and new together),
(3) repeat for many generations, then compare the surviving particles' spatial power spectrum
against a never-merged uniform reference of the same count.

- **Our algorithm has no cell-frequency aliasing.** Confirmed at multiple harmonics, including a
  deliberately "safe" high-order one (17 cells/patch, 11×11 patches, so the cell frequency is the
  17th patch harmonic and shares no small-integer relationship with the patch tiling or domain
  size) — spike factor 1.16×, z=1.59, indistinguishable from noise. This makes sense mechanically:
  merge position is always a true weighted pairwise centroid, never snapped to any grid point.
- **Our algorithm *does* have a patch-frequency artifact, but only tied to the true, finite
  simulation domain boundary — not the interior.** A finite (non-periodic) domain showed a strong,
  unambiguous spike at the patch tiling frequency (z up to ~54) because corner/edge patches have
  fewer real neighbors than interior ones and complete fewer merges per round as a result. Switched
  to periodic BCs (removing the true edge, making every patch topologically identical) at matched
  grid sizes (5×5 through 20×20): the spike collapsed by roughly an order of magnitude and further
  scaling up (10×10, 20×20) showed no growth — consistent with a boundary-of-the-whole-domain
  effect, not something that compounds throughout a large interior.
- **Reference/control points measured for calibration**: a classic cell-based merge that snaps to
  the cell center is the true pathological case (z≈290–330 at the cell frequency — two orders of
  magnitude worse than anything our algorithm produces). A milder scheme — nearest-neighbor
  pairing *restricted to within one cell*, but positioned at the true centroid (not snapped) — does
  show a real but much smaller cell-frequency artifact (z≈12–15 at standard cycle counts), which
  *shrinks* rather than grows with more generations. Our algorithm never restricts its search to a
  single cell, so this doesn't apply to it directly, but it's a useful calibration point for how
  much search-scope restriction alone (independent of a bad position rule) can bias things.

**Methodological note, corrected after a follow-up check**: an earlier version of this note
warned against using a small-integer cells-per-patch ratio (e.g. 4 or 5), reasoning that a
coincidence between the cell frequency and a low patch harmonic would make a real cell-frequency
effect and a patch-harmonic artifact numerically indistinguishable — this was true, but only
because the original coincidence (5 cells/patch, 5×5 finite domain) was checked in a *finite*
domain, where low patch harmonics still carry real power. Re-tested directly at 4 cells/patch on
a much larger *periodic* 40×40 domain: the cell frequency read completely clean (z=−0.11), no
ambiguity at all. The reasoning failure was generalizing from the finite-domain case — under
periodic BCs, patch-frequency power is already near baseline at every harmonic (see above), so
there's no real patch signal left to be confused with, regardless of which harmonic the cell
frequency happens to coincide with. **The coincidence only matters when checking a finite,
non-periodic domain; it is a non-issue under periodic BCs, at any cells-per-patch ratio.** Any
future spatial-aliasing check on a finite domain should still pick a cells-per-patch ratio with no
low-order common structure with the patch grid (or just run under periodic BCs, which sidesteps
the question entirely).

## Explicit scope limitations

- Single pass per round (see above) — a deliberate, accepted conservatism, not a defect.
- The propose/judge/verdict protocol (a two-message round: propose, then verdict) fully replaced
  an earlier idea of needing a separate third confirmation round — that earlier idea is obsolete,
  not a deferred TODO.
- Validated so far against synthetic domain topologies (2D and 3D corners, mixed same-rank/
  cross-rank ownership, a genuine T-junction) in a Python/mpi4py prototype, not yet against the
  real `ParticleSoA`/`AmrMesh`/ghost-fill C++ machinery. Translating this rule set into that
  machinery (including the real `ParticleGhostMask`-driven ghost-fill, real point-to-point/`alltoall`
  routing, and the `NNMergeParticle`/proposal wire format) is the next step.
