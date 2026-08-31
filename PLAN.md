# Issue #706 — grad(D) drift correction for the Ito/particle transport

Branch: `claude/706-grad-d-drift` (worktree `.claude/worktrees/issue-706`, base `main` @ 813be50d)

## 1. Problem restatement

`ItoSolver`/`ItoKMC` integrate the Ito-interpretation Euler-Maruyama update

```
X^(k+1) = X^k + dt*v + sqrt(2*D*dt)*N(0,1)
```

whose Fokker-Planck equation is `dn/dt = -div(v*n) + grad^2(D*n)`, i.e. flux `v*n - n*grad(D) - D*grad(n)`.
`CdrSolver` integrates `dn/dt = -div(v*n) + div(D*grad(n))`, i.e. flux `v*n - D*grad(n)`.

The two agree iff `grad(D) = 0`. Recovering the fluid form from an Ito SDE requires the **drift** to be

```
a = v + grad(D)          (scalar isotropic D; for a tensor it would be v + div(D))
```

We take resolution 1 from the issue (fix it properly), enabled by default.

## 2. Where the difficulty is

`ItoKMCGodunovStepper` is semi-implicit in the field. Its update is split as

```
X^(k+1) = X^k + [ explicit displacement, known at time k ]  +  dt * mu^k * E^(k+1)(X^k)
```

The semi-implicit Poisson operator

```
div( (eps + dt*sigma^k) grad(phi^(k+1)) ) = -rho^dagger
```

is exact only because (a) `rho^dagger` is deposited from particles displaced by *everything except* the
`dt*mu*E^(k+1)` term, and (b) that remaining term is **linear in `E^(k+1)`**. Adding `grad(D)` to the particle
velocity columns would put a non-`E`-proportional piece into the implicit half and break the operator.

**`grad(D^k)` is explicit — it is a lagged, field-independent displacement.** So the correction belongs with the
stochastic hop, in the `rho^dagger` displacement. This is exactly the pattern the CDR half of the same stepper
already uses: `computeDiffusionTermCDR` (`CD_ItoKMCGodunovStepperImplem.H:1763`) folds the explicit
`dt*div(D*grad(phi^k))` straight into `m_semiImplicitRhoCDR`. Nothing about the Poisson operator changes.

Concretely:

| term | current | after |
|---|---|---|
| `rho^dagger` particle position | `X^k + hop` | `X^k + hop + dt*grad(D)` |
| `scratch_x/y/z` (explicit displacement) | `hop` | `hop + dt*grad(D)` |
| `stepEulerMaruyamaParticles` | `X = X_old + f*v + g*scratch` | **unchanged** |
| conductivity `sigma` | `sum |Z|*e*w*mu` | **unchanged** |
| Poisson operator | `div((eps + dt*sigma)grad phi)` | **unchanged** |

`f = mobile ? dt : 0`, `g = diffusive ? 1 : 0` already, so a diffusive-but-immobile species (neutral excited
states) picks the correction up for free. That property is the main reason to lump into the hop rather than into
the velocity — see the rejected alternative in §7.

## 3. Where grad(D) comes from

`ItoSolver::m_diffusionFunction` (`CD_ItoSolver.H:1438`) is a 1-component `EBAMRCellData` on the particle realm,
already coarsened and ghost-filled by `ItoKMCStepper::computeDiffusionCoefficients`
(`CD_ItoKMCStepperImplem.H:3277-3288`) immediately before `solver->interpolateDiffusion()`.

`AmrMesh::computeGradient(EBAMRCellData&, const EBAMRCellData&, realm, phase)` (`CD_AmrMesh.H:455`) gives a
`SpaceDim`-component cell-centred gradient through `EBGradient`, which handles cut cells and coarse-fine
boundaries and is the same operator `FieldSolverGMG` uses for `E = -grad(phi)`. `s_eb_gradient` is registered by
default in every `PhaseRealm` (`CD_PhaseRealm.cpp:35`), so no new registration is strictly required — we add an
explicit `registerOperator(s_eb_gradient, m_realm, m_phase)` in `ItoSolver::registerOperators` anyway so the
dependency is not silent. `EBGradient` asserts `a_phi.ghostVect() == m_ghostVector`, which the default
`m_amr->allocate` already satisfies.

Interpolation to particles uses the same path as `ItoSolver::interpolateVelocities(lvl, dit)`
(`CD_ItoSolver.cpp:2506`): `EBParticleMesh::interpolate<D_DECL(&ItoParticle::scratch_x, ...)>` with `m_deposition`
and `m_forceIrregInterpolationNGP`.

## 4. Storage — no new types, existing containers only

* **Per-particle**: reuse the existing `ItoParticle::scratch_x/y/z` columns (`CD_ItoParticle.H:49-53`). Their only
  current consumer is the Godunov hop (`CD_ItoKMCGodunovStepperImplem.H:1748` write, `:1842` read); the
  `m_extendConductivityEB` flag uses the *scalar* `scratch`, a different column. No new payload columns, so no
  extra `SpaceDim*sizeof(ParticleReal)` per particle in a 10^8-particle run, and nothing new to checkpoint.
  The gradient is interpolated into `scratch_x/y/z` and then overwritten in the same particle loop by
  `hop + dt*gradD`, so the lifetime is a handful of instructions and never spans a remap, merge or regrid.
* **On the mesh**: one `SpaceDim` cell-data holder. See open question **A** — per-solver persistent vs one shared
  scratch.

Neither is a new type or a new container kind; both are `EBAMRCellData`/existing payload columns.

## 5. Staged implementation

### Stage 0 — make the mesh diffusion field honest in BrownianWalker (prerequisite)

`BrownianWalkerStepper::setDiffusion` (`CD_BrownianWalkerStepper.cpp:134`) writes
`std::numeric_limits<Real>::max()` into `m_solver->getDiffusionFunction()` and sets the particle coefficients
directly from the scalar `m_diffCo`. Any generic `grad(D)` machinery reading that field would produce garbage.
Replace with `setDiffusionFunction(m_diffCo)` plus `conservativeAverage` + `interpGhostPwl`. With a constant
`m_diffCo` the gradient is then identically zero and Stage 3 is a no-op for every existing BrownianWalker test —
which is also how we show the new code path is neutral.

### Stage 1 — `ItoSolver`: mesh gradient + interpolation

* New option in `CD_ItoSolver.options`, default on:
  `ItoSolver.diffusion_gradient_drift = true   ## Add grad(D) to the drift so the SDE transports D*grad(n)`
* `computeDiffusionGradient()` — `m_amr->computeGradient(...)` on `m_diffusionFunction`, then
  `conservativeAverage` + `interpGhostPwl` on the result. No-op when `!m_isDiffusive` or the flag is off.
* `interpolateDiffusionGradient()` and `interpolateDiffusionGradient(int a_lvl, const DataIndex& a_dit)` —
  interpolate onto `scratch_x/y/z`, mirroring `interpolateVelocities` exactly (including the irregular-cell NGP
  override). Zero the columns when disabled, so callers can add unconditionally.
* `isDiffusionGradientDrift()` query.
* `registerOperator(s_eb_gradient, m_realm, m_phase)` in `registerOperators`.
* Optional: a `grad_dco` entry in `ItoSolver.plt_vars` for debugging. Cheap and worth it while validating.

### Stage 2 — `ItoKMCGodunovStepper`

In `diffuseParticlesEulerMaruyama` (`CD_ItoKMCGodunovStepperImplem.H:1701`):

1. Per solver, before the level loop and only when `diffusive && solver->isDiffusionGradientDrift()`:
   `solver->computeDiffusionGradient()`.
2. Inside the patch loop, before the particle loop: `solver->interpolateDiffusionGradient(lvl, din)`.
3. In the particle loop: read `gradD` from `scratch_x/y/z` into a local `RealVect` **first**, compute
   `hop = diffusionFunction(p, dt)`, then write `disp = hop + dt*gradD` back into `scratch_x/y/z` and append the
   `rho^dagger` particle at `pos + disp`.

`stepEulerMaruyamaParticles` (`:1803`) needs no change. `advanceEulerMaruyama`'s
`remapPointParticles(..., SpeciesSubset::ChargedDiffusive)` already covers every species that now moves further.
Between `diffuseParticlesEulerMaruyama` and `stepEulerMaruyamaParticles` only the *point* particles are remapped,
so the `scratch` columns on the real particles survive untouched.

Note the ordering hazard: `p = leaf.gather(i)` is taken before we overwrite the scratch columns, and
`forwardIsotropicDiffusion` reads `p.vx/vy/vz`, not scratch — so gathering first is safe either way, but the
local-`RealVect`-first ordering must be explicit in the code, not incidental.

### Stage 3 — `BrownianWalker`

`BrownianWalkerStepper::advance` (`:457`) currently applies the drift only under `isMobile()`. Add, under
`isDiffusive()`, a `computeDiffusionGradient()` + `interpolateDiffusionGradient()` before the patch loop and
`pos += dt*gradD` in the Euler step. With Stage 0 in place and a constant `m_diffCo` this is exactly zero, so
`Exec/Tests/BrownianWalker/DriftDiffusion` and the merge tests must be bit-comparable before/after.

### Stage 4 — time step

`ItoSolver::computeAdvectiveDt` bounds `dt` from the particle velocity columns, which will not contain `grad(D)`.
For realistic streamer parameters the omission is not the binding constraint: with `D ~ 0.1 m^2/s` varying over
`~10 um`, `|grad D| ~ 1e4 m/s`, and at `dt = 1 ps` the correction displacement is
`dt*|gradD| / sqrt(2*D*dt) ~ 2e-2` — two percent of the hop the diffusive `dt` rule already bounds. Plan:
document this in the doxygen for the new functions and in `Ito.rst`; optionally add a cheap
`DataOps::getMaxMin` bound on `|grad D|` folded into the Ito advective `dt` in `ItoKMCStepper`. See open
question **D**.

### Stage 5 — verification

This is the part the issue specifically calls out as missing ("invisible in any test with a constant diffusion
coefficient"). The decisive test:

> **Uniform density, spatially varying D, zero drift.** The fluid equation `dn/dt = div(D grad n)` has `n = const`
> as an exact steady state for *any* `D(x)`. The uncorrected Ito scheme does not: particles pile up where `D` is
> small. With the correction, uniformity is preserved to statistical noise.

Home: a new `Exec/Tests/BrownianWalker/VariableDiffusion` (or `Exec/Convergence/BrownianWalker`). Needs
BrownianWalker to support a spatially varying `D` — a linear ramp or a Gaussian blob set through
`DataOps::setValue(..., lambda, probLo, dx, ...)` in `setDiffusion`, with the particle coefficients taken from
`solver->interpolateDiffusion()` instead of `setParticleDiffusion(m_diffCo)`. Measure the RMS deviation of the
deposited density from its initial uniform value; it should stay at the sqrt(N) noise floor with the correction
and grow linearly in `t` without it. See open question **E** on whether to build this now.

Secondary checks:
* `Exec/Examples/ItoKMC/ComparisonCdrPlasma` — rerun both models with the tabulated `D(E/N)` and confirm the head
  positions converge. Update the "will not agree exactly" paragraph in its `README.md:43-47`, which currently
  tells the reader to freeze `D` and points at #706.
* All existing `Exec/Tests/ItoKMC/*` regression inputs run to completion with no NaN/assert.

### Stage 6 — documentation

* `Docs/Sphinx/source/Solvers/Ito.rst` — equation `ito_diffusion` (line 9) currently states the plain Ito SDE.
  State the drift as `V + grad(D)`, explain the Fokker-Planck equivalence to the CDR flux, and document the new
  option.
* `Docs/Sphinx/source/Applications/ItoKMC.rst` — the semi-implicit split, and that `grad(D)` sits in the explicit
  half.
* Doxygen on `ItoKMCPhysics::isotropicDiffusion` (`CD_ItoKMCPhysicsImplem.H:908`) — it returns the stochastic hop
  only; the drift correction is applied by the stepper.
* **`literalinclude` drift.** `Ito.rst` pulls absolute `:lines:` ranges out of `CD_ItoSolver.H` and
  `CD_ItoParticle.H`. Adding declarations shifts every range below the insertion point. Ranges at risk, in file
  order: `691-697`, `714-719`, `791-797`, `809-814`, `824-829`, `862-868`, `917-928`, `954-964`, `1052-1057`,
  `1076-1081`. Adding the new methods next to `interpolateDiffusion` (`CD_ItoSolver.H:829`) shifts everything from
  `862` down. Follow CLAUDE.md §4 and rebuild the docs:
  `cd Docs/Sphinx && python3 -m sphinx -W --keep-going -b html source build/html`.

## 6. Known limitations to state explicitly (not fix)

* **Energy-parametric D (LEA).** `ItoSolver::updateDiffusion()` sets `D` from the particle energy rather than from
  a mesh field. `m_diffusionFunction` is then not the source of truth and its gradient is meaningless. In practice
  nothing in the tree calls `updateDiffusion()` — `ItoKMCStepper` always goes through `interpolateDiffusion()` —
  so this is a documented limitation, not a live bug. (A route exists via `depositDiffusivity`, but it is noisy
  and out of scope.)
* **Anisotropic D.** The correction is `grad(D)` for scalar isotropic `D`; a tensor would need `div(D)`.
* **Near the EB.** `E` has a boundary layer at the electrode, so `D(E/N)` and hence `grad(D)` can be large and
  poorly resolved in cut cells. The correction is a physical drift, so this is not an instability, but it is worth
  watching in the first `ComparisonCdrPlasma` run. If it misbehaves, the mitigation is a limiter on
  `|dt*grad(D)|` (as a fraction of `dx` or of the hop) — deliberately *not* in the initial plan, to avoid adding a
  knob before there is evidence it is needed.
* **Discrete, not exact, agreement.** The corrected particle scheme reproduces `div(D grad n)` in the continuum
  limit; the CDR solver's face-centred `D` discretisation differs at O(dx). Statistical noise and superparticle
  merging remain the dominant difference in `ComparisonCdrPlasma`.
* **Enabling by default changes results.** Every existing ItoKMC run with a tabulated `D(E/N)` moves. Must be
  called out in the PR description.

## 7. Alternative considered and rejected

**Put `grad(D)` into the particle velocity columns** (`interpolateVelocities` returns `mu*E + grad(D)`), so it is
literally the drift velocity, `computeAdvectiveDt` sees it for free, and BrownianWalker needs no change at all.
Rejected because:

1. `interpolateVelocities` early-returns on `!m_isMobile`, and `stepEulerMaruyamaParticles` uses
   `f = mobile ? dt : 0`. A **diffusive-but-immobile** species — neutral excited states in ItoKMC — would silently
   lose the correction. Fixing that means allocating `m_velocityFunction` for immobile solvers, running
   `interpolateVelocities` for them, and changing `f`; that is a wider blast radius than the whole rest of this
   change.
2. The Godunov stepper still needs `grad(D)` *separately* for the `rho^dagger` displacement, so the term would
   have to be interpolated twice and carefully *not* double-counted between the explicit and implicit halves.
3. `setItoVelocityFunctions` only fills the velocity function for `Z != 0`, so neutral diffusive species are
   excluded twice over.

The `dt` benefit (point §4 Stage 4) is the only thing lost, and it is a two-percent effect on a constraint the
diffusive rule already dominates.

## 8. Open questions

**A. Mesh storage for the gradient.** Two options, both existing containers:
   * *Per-solver persistent* `ItoSolver::m_diffusionGradient` (`SpaceDim`, allocated only when diffusive and the
     flag is on). Same shape as the existing `m_velocityFunction`. Self-contained, plottable, works for
     BrownianWalker without extra plumbing. Costs `SpaceDim` components per diffusive species — for a 10-species
     3-D ItoKMC run, ~30 extra doubles/cell.
   * *One shared scratch*, reused species-by-species. `ItoKMCStepper::m_particleScratchD` (`CD_ItoKMCStepper.H:787`)
     already exists and is exactly the right shape; BrownianWalker would need its own. O(1) memory.
   **Recommendation: per-solver persistent**, for the self-containment and plottability, unless the memory matters
   to you — the ItoKMC species count is the deciding factor and you know it better than I do.

**B. Confirm the scratch-column reuse.** §4 reuses `ItoParticle::scratch_x/y/z` rather than adding
   `gradD_x/y/z` payload columns. This keeps the per-particle footprint unchanged but means two consumers share
   one scratch group within a single function body. Dedicated columns would be self-documenting at the cost of
   `SpaceDim` reals per particle everywhere. **Recommendation: reuse scratch.**

**C. Where the option lives.** `ItoSolver.options` (per-solver, next to `mobility_interp`, and automatically
   available to BrownianWalker) vs. an `ItoKMCGodunovStepper` option. **Recommendation: `ItoSolver.options`.**

**D. Time-step bound.** Document only, or also fold a mesh `max|grad D|` into the Ito advective `dt`?
   **Recommendation: document now, add the bound only if a run shows it matters.**

**E. Verification test scope.** Build the varying-`D` BrownianWalker test in this PR (it needs a spatially varying
   `D` mode in BrownianWalker, which is real but contained work), or land the fix first and add the test
   separately? **Recommendation: build it in this PR** — the issue's own argument is that the defect survived
   because no test varies `D`, and without it we have no evidence the sign and magnitude are right.

**F. Default on.** You said enabled by default. Confirming that you accept that every existing ItoKMC result with
   a tabulated `D(E/N)` shifts, and that this goes in the PR description.
