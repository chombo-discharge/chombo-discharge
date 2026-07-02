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
