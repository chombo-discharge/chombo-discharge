# Review of PLAN.md revision 5 — findings

**Reviewed** `Prototypes/MirrorDeposition/PLAN.md` at revision 5 · **Branch** `mirror_deposition`
· **Code citations re-resolved against** `dab94cbd` · **All harnesses re-run from this directory.**

This is review 5, following the four adversarial reviews tabulated in PLAN.md §0. It follows the
review order that section prescribes: re-run every table, re-resolve every code citation, follow the
deposit past the deposit, attack §2 and §3, and do not re-derive widths along the normal.

---

## 0. Verdict

Revision 5 is materially sound. **Every number in the file reproduces exactly** — not one table
row is off — and every derivation I could check independently held, including the ones review 4
corrected. Three real findings survive, all of the class §0 predicts: a claim whose *scope* is
narrower than its *statement*.

Act on **Finding 2** first. It is the only one of the three that can produce a wrong answer with a
clean Jacobian and a passing mass check, and the harnesses in this directory are structurally
incapable of detecting it.

---

## 1. What reproduced

### 1.1 Tables

All twelve reproduce exactly.

| Section | Table | Status |
|---|---|---|
| §1.2 | Band mask retention, both kernels | exact |
| §1.2 | Kernel-dependent Chebyshev prediction, four axis cases | exact |
| §1.2 | Worked example at `kappa = 0.05` | exact, all five rows |
| §1.2 | Image-cell reach (`img_reach.py`, new in r5) | exact |
| §2.4 | Reflection-source ladder, seven spheres | exact |
| §3.1 | Curvature error without the Jacobian | exact |
| §3.3 | Curvature from the discrete normal, 5³ | exact |
| §3.3 | Normal-noise sensitivity, 3³/5³/7³, both signs | exact |
| §6.1 | `sum(kappa phi dV)/sum(m)` versus the sliver bin | exact |
| §6.1 | Images-only (zeroing) build | exact |
| §6.1 | Scoring-margin sensitivity | exact |
| §7 | Tight-cavity Jacobian and reflect fractions | exact — but see Finding 3 |

Spot values, for the next reviewer to diff against:

- **§1.2** Chebyshev 2 → 100.000% for both kernels on all ten cases. Chebyshev 1 → CIC 91.596%,
  TSC 87.163%, both on `axis, s = 0.45`. Chebyshev 0 → 5.104% / 5.093%. Max Chebyshev is 2 for
  both. Deepest reflecting particle `d = 2.570` (CIC) / `3.390` (TSC). `axis, s = 0.12` is
  Chebyshev 1 under CIC and 2 under TSC, with Chebyshev-1 retention 99.617%.
- **§1.2** `img_reach.py`: image→cut max 1 for both kernels; image→source max 4 (CIC) / 5 (TSC);
  the "ghost if no remap" column maxes at 5 and 6, which is the `4 + 1` / `5 + 1` in the text.
- **§6.1** The images-only build reads −97.54% to −100.00% in the regular bin and −1.82% to
  −7.69% in the `kappa <= 0.05` bin.
- **§6.1** The "~1.8% deficit at margin 1" is `45xy, s max` (−0.01807) and `diag, s max`
  (−0.01825). There is genuinely no plateau at margin 3.
- **§7** Cavity `Jmax` 179.1 at `R = 3`; worst binned deviation 2.16% (at `R = 3.5`). Reflect
  fractions 62.7 / 71.3 / 91.2 / 98.0%.
- **§8** The self-reported README discrepancy is itself correct: the harness gives
  +0.00686 / +0.01108 / +0.02079 / +0.03559 against the README's −0.007 / −0.011 / −0.019 / −0.025,
  with sample sizes 235→24 matching exactly.

One note on §3.1: the "+117% / +40% / +15.3% / +4.2%" figures quote the **worst** `kappa` bin, which
is the sliver bin at `R = 4, 8, 16` but the `0.15–0.30` bin at `R = 40` (+0.0424 against the
sliver's +0.0219). The numbers are right; the selection rule is not stated.

### 1.2 Counts

Every count in the file is exact.

| Claim | Verified |
|---|---|
| 93 C++ occurrences, 8 files | 49 + 18 + 10 + 10 + 2 + 2 + 1 + 1 = 93 |
| 85 `a_forceIrregNGP` + 4 `forceIrregNGP` + 4 `m_forceIrregDepositionNGP` | exact |
| 7 `m_forceIrregInterpolationNGP` out of scope | exact |
| `irr_ngp_deposition`: 29 occurrences in 29 files, plus `Ito.rst` | exact |
| `eb_ghost` 39×2 / 100×4 | exact |
| `num_ghost` 71×2 / 68×3 | exact |
| `CD_AmrMesh.options` ships 2 for both | exact |
| `ItoSolver.redistribute` 7 false / 23 true (30 total) | exact |
| `ItoSolver.blend_conservation` 0 true | exact |
| `McPhoto.blend_conservation` 14 false / 0 true; no `redistribute` key | exact |
| `CdrGodunov` 7 files / `CdrCTU` 34 files | exact |
| 8 `m_amr->deposit*` callers in 5 files | exact |
| 9 ad-hoc `EBParticleMesh` constructions | exact, all nine lines |
| `AmrMesh::depositParticles<Members...>` has zero callers | exact |
| `AmrMesh.buffer_size` ships 2 | exact |

Beware when re-running these: a git worktree at `.claude/worktrees/issue-696` shadows the whole
tree and doubles every raw `grep -r` count. Filter it out.

### 1.3 Derivations and code facts

Checked independently, not merely re-read:

- **CIC is a top-hat cloud**, `weight = max(0, min(b,L) - max(a,-L))` with `L = 0.5*a_particleWidth`
  (`CD_EBParticleMesh.H:732`); **TSC is triangular** with `L = a_particleWidth` (`:794`). With
  `cicWidth = 1*widthScale` and `tscWidth = 2*widthScale` (`:474-475`), the cloud half-widths are
  `h = 1/2` and `h = 1` — so §1.2's `kappa < 0.25` (CIC) and `kappa < 0.5` (TSC) thresholds follow
  directly.
- **The band derivation.** `d_max = t + h` and `x_p < c + 2t + h` both follow from requiring
  `x_img + h > c`; the particle leaves cell `c+1` when `2t + h > 2`; and `x_p < c + 3` in both
  cases, so Chebyshev 2 is the bound for both kernels. Confirmed.
- **The `kappa = 0.05` worked example.** All five rows recomputed by hand — positions, images,
  cells, Chebyshev distances and overlap weights. Correct.
- **The 1-D image-cell argument.** `x_img >= c - h` and `x_img < a < c+1` give image cell
  `c-1` or `c` for both kernels. Correct.
- **§4.1's estimator.** `phi_0 = 0.5n` under plain CIC and exactly `n` under the mirror, verified by
  integrating the top-hat overlap. The `V = 1.6` (wall at 0.3, particle at 0.4) and `V = 1.9`
  (wall at 0.5, particle at 0.6) counterexamples both recomputed and correct — so the issue's
  per-particle normalization really is wrong.
- **§2.2 and §3.1 signs.** Verified independently against a convex sphere (fluid outside) and a
  concave cavity (fluid inside), from `S = +dn̂` with `n̂` into the fluid: `S = +I/R` convex,
  `-I/R` concave; `d = eta + (1/2) xi^T S xi` reproduces the exact signed distance to second order
  in both; `nhat = n_c + S xi` gives the foot-point normal in both; `J = (R-d)^2/(R+d)^2` convex and
  `(R+d)^2/(R-d)^2` concave, which are the correct shell-area ratios. The numerator sign change past
  `d = 1/c_max` and the denominator zero at `d = -1/c_i` are real and are the same event.
- **`Location::position`** does return a physical position — `ret *= a_dx` — so the §2.1 warning
  about `bndryCentroid` being cell-relative is correct. (Citation range is short; see nit 2.)
- **The valid holder is restored after every deposit.** Not cited in the plan but load-bearing for
  §5.3: `depositHaloCore` copies then clears (`:849`, `:912`); `depositHaloNGPCore` and
  `depositTransitionCore` transfer out and both call
  `particles.transferParticles(particles.getMaskParticles())` (`:988`, `:1294`).
- **`getTransitionMaskWidth`** really is `refRat/2` for CIC and `refRat` for TSC (`:1001-1006`).
- **O7 is exactly right**: `remap()` drops particles with no valid destination and counts them into
  `m_numOutcastLocal`, readable via `getNumberOfOutcastParticles{Local,Global}()`.

---

## 2. Findings

### Finding 1 — surface data can be undefined on a valid band cell at a level's grid boundary

**Severity: high. Not covered by any open item.**

§5.4's exchange-then-extend ordering fixes the *inter-patch* case. It does not fix the
*outside-the-level's-grids* case, and these are different problems with different fixes.

The band mask is built from `EBISBox`, which is valid `eb_ghost` cells out **regardless of whether
any box on that level owns those cells** — that is exactly what
`new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_numEbGhostsCells, &(*m_ebis))`
(`CD_PhaseRealm.cpp:411`) buys. `LevelData::exchange()`, by contrast, fills ghost cells only from
*valid* regions of the same level.

So on a fine level whose grid edge passes within 2 cells of the EB, a valid fine band cell can be
marked — its nearest cut cell is a ghost cut cell, plainly visible through `EBISBox` — while no fine
box ever fits an `S` for that cut cell and no exchange can deliver one. `(x_c, n_c, S_c)` is then
whatever the allocation left behind, and the particle reflects off it.

This is adjacent to **O9** but is not the same item. O9 is about images being *lost* after they are
correctly built. This is about images being *built from uninitialized state*. With
`AmrMesh.buffer_size` shipping 2 and the band being 2, it is not a corner case.

**Fix.** The band-mask predicate must require that the nearest cut cell is one for which surface
data was actually delivered — not merely one that `EBISBox` can see. Add an explicit no-reflect
fallback with a counter, in the same spirit as the crease fallback in §5.4 step 2. State it in §5.4,
because "extend to every band cell from the nearest cut cell" currently reads as total.

### Finding 2 — the tangent frame is part of the stored format, and no harness exercises it

**Severity: high. Silent failure mode. Act on this one first.**

§2.1 stores `S_c` as 3 Reals, "2×2 symmetric, tangent frame". §3.3 says the fit uses "**any**
orthonormal basis perpendicular to `n_i`". §2.2 says the reflection uses "`t1, t2` spanning the
plane perpendicular to `n_c`". Those two bases must be the *same* basis, not merely both
orthonormal. Otherwise `xi = (t1.w, t2.w)` is expressed in one frame while `S` lives in another, and
both `d` and `nhat` come out wrong.

The failure mode is the worst available: `2H = tr S` and `K = det S` are **frame-invariant**, so a
frame mismatch leaves the Jacobian *correct* while the reflection point drifts. §6.1's secondary
mass check would not move. Nothing looks wrong.

**The harnesses cannot catch this, structurally.** `mirror_source.py:quad_reflect` takes
`kap_princ` — a scalar, its own docstring says `k1 = k2 (sphere)` — so it computes `S = c·I`, which
is frame-*invariant*. The fitted anisotropic 2×2 `S` produced by `mirror_discrete_curvature.py` is
never fed into the reflection by any harness in this directory. The two halves of the scheme are
each measured; their composition is not.

Two consequences:

1. **§2.4's ladder validates the quadratic patch only on umbilic surfaces.** "The quadratic patch is
   <= 2.9% everywhere" is measured with an isotropic `S` on seven spheres. The anisotropic path —
   the whole reason `S` is a 2×2 rather than a scalar — is unexercised end to end. §7's
   "a torus can" note frames the gap as being about *borrowing* `(x_c, n_c, S_c)` across cells; the
   anisotropy gap is separate and larger.
2. **The frame rule is currently an accident, not a specification.** Both harnesses do happen to use
   the same rule (`a = e_x` unless `|n_x| > 0.9`, then `e_y`; `t1 = n × a` normalized, `t2 = n × t1`)
   — `mirror_discrete_curvature.py:125-130` and `mirror_source.py:61-65`. So a consistent
   implementation is achievable. But nothing in PLAN.md says it must be, and the rule branches on a
   threshold, which is precisely the kind of thing that fails to reproduce bit-for-bit once `n_c` has
   been through a store, an exchange and a renormalization.

**Fix, pick one.** Either (a) write the `t1`/`t2` construction into §2.1 as a normative part of the
stored format, with the requirement that it be a pure deterministic function of `n_c` invoked
identically at fit and at use; or (b) store `t1` explicitly — 12 comps instead of 9 — and stop
depending on the branch reproducing. Independently of which, **§6.1's curved case must be the torus,
not a sphere**, and it should be checked that it actually produces anisotropic `S` in the sliver
bin.

### Finding 3 — §7's cavity numbers are CIC-only, computed with a band §9 retired

**Severity: medium. Same class as harness bug 5, one level up.**

`mirror_cavity.py` hardcodes `BANDMAX = 1.5*sqrt(3) = 2.598` and selects reflecting particles with
`d <= 1.5*DX*|n|.sum()` — the `3*s_max` distance band that §9 lists as **superseded**. Under the
plan's actual criterion, `reach_cells.py` measures `d` out to **3.390** for TSC.

Re-running the analytic amplification `J = ((R+d)/(R-d))^2` at the real TSC reach:

| cavity `R/dx` | `J` at `d = 2.598` (as tabled) | `J` at `d = 3.390` (TSC) |
|---|---|---|
| 4.0 | 22.2 | **146.8** |
| 3.5 | 45.7 | **3923** |
| 3.0 | 194.0 | **singular** |

So §7's *"Down to `R = 3*dx`, where `J` reaches 179, the worst binned deviation stays under 2.2%"*
is a CIC statement presented as a general one. The deposition runs beneath it exclude the deepest
reflecting particles by the same superseded criterion, so the "under 2.2%" figure is measured on a
sample the real scheme would not produce.

The stated mechanism is also not the harness's mechanism: §7 says the high-`J` images *"land where
every stencil cell is covered and are dropped"*, but the harness drops them by a distance test
before that ever comes into play.

**Fix.** Either re-run `mirror_cavity.py` against the plan's real criterion for both kernels, or
label the row as CIC-only and state the TSC singular radius explicitly. The conclusion
("a NaN risk, not an accuracy risk") probably survives — but the guard in §3.2 has to be sized for
`d = 3.39`, not `2.598`.

**Related, and worth fixing before the next review:** `mirror_band_kernels.py` still prints a
`plan band SHORT by …` verdict column against that same retired `3*s_max` / `4*s_max` formula. It
refutes nothing in revision 5, since the plan no longer states a band as a distance. But the next
reviewer will read eight "SHORT by" verdicts as a live refutation and spend the review on it. Retire
or relabel the column.

---

## 3. Nits

1. **§5.1, "seventeen of them, `:605-695`"** — there are **eighteen** `mutable` members in that
   range: `:605, 610, 615, 622, 630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690,
   695`. The three named in the text (`m_redistributionOp :665`, `m_particleMesh :675`,
   `m_nonConservativeDivergence :695`) are correct, as is `m_irregularCells :600` not being
   `mutable`.

2. **§2.1, "`Location::position` is `CD_LocationImplem.H:37-41`"** — `:37-41` is only the
   `Location::Cell::Boundary` switch case. The `ret *= a_dx` that makes the return a *physical*
   position is at **`:49`**, outside the cited range. Since the citation exists precisely to nail
   the units that revision 3 got wrong, it should read `:37-51`.

3. **§6, "`depositHybrid` (`:2193`)"** — `:2193` is `void`; the function name is at **`:2194`**. The
   `deltaM = (1-kappa)*(dc - kappa*dnc)` line and its *"Remember, `dc` already scaled by kappa"*
   comment are at **`:2234`**, one past the cited `:2225-2233`.

4. **§4.3's "before any deposit" contradicts §5.3's pseudocode**, which deposits the real particles
   first. The pseudocode is the correct one — but for a reason the plan never states: all three
   coarse-fine cores restore the valid holder before returning (`CD_EBAMRParticleMesh.H:912, 988,
   1294`). The operative constraint is *"from the valid holder, never from a mask holder"*, which is
   order-independent. Say that, and cite the restore. As written, the next implementer will reorder
   §5.3 to obey §4.3 and believe they fixed a bug.

5. **§3.1's error figures quote the worst `kappa` bin, not the sliver bin** (§1.1 above). Worth one
   clause, since every other table in the file bins by `kappa` and reports the sliver.

---

## 4. Reproducing this review

From `Prototypes/MirrorDeposition/`:

```
python3 reach_cells.py                          # §1.2 band and retention
python3 img_reach.py                            # §1.2 image reach
python3 mirror_source.py                        # §2.4 ladder        (~1 min)
python3 mirror_discrete_curvature.py fit        # §3.3 fit
python3 mirror_discrete_curvature.py noise      # §3.3 noise sweep
python3 mirror_zeroing.py                       # §6.1 zeroing
python3 mirror_cavity.py                        # §7 cavity          (see Finding 3)
python3 mirror_sphere_ext.py radii              # §3.1, §6.1         (~2 min)
python3 mirror_planar_ext.py margin             # §6.1 margin
python3 mirror_band_kernels.py                  # slow; stale verdict column
```

Code citations were re-resolved with a per-file `sed -n '<line>p'` sweep rather than by text search,
and counts with `grep -r … | grep -v '\.claude/worktrees'`.

**This file is throwaway and must be deleted before the branch merges, along with the rest of this
directory.**
