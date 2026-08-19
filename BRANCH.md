# Branch: `claude/cdrplasma-itokmc-comparison`

A side-by-side comparison of the `CdrPlasma` and `ItoKMC` models on the same streamer problem, and the
bug hunt that came out of it. The two models initially disagreed on streamer velocity, radius, and head
field; they now agree on velocity and radius.

The example lives in `Exec/Examples/ItoKMC/ComparisonCdrPlasma`.

---

## 1. The example

A rod-needle electrode at 35 kV in a grounded 8 cm x 8 cm box of atmospheric air, 2D only, run with

* `Physics/CdrPlasma` -- `CdrCTU` transport, `McPhoto` radiative transfer, deterministic reaction integrator.
* `Physics/ItoKMC` -- `ItoSolver` transport, `McPhoto` radiative transfer, kinetic Monte Carlo chemistry.

selected by a single switch, `Comparison.model = cdr_plasma | ito_kmc`. Geometry, mesh, voltage, transport
data, chemistry, and initial particles are shared. See the example's `README.md` for the file layout and how
to run it.

**The electron diffusion coefficient is deliberately held constant** (`0.1373223` m^2/s, the tabulated value
at E/N = 150 Td) in *both* chemistry files. This is not a modelling preference -- see section 4.

---

## 2. Bugs found

### Fixed on this branch

**Photon production was not accumulated across chemistry substeps** -- upstream issue
[#702](https://github.com/chombo-discharge/chombo-discharge/issues/702).

`CdrPlasmaJSON::integrateReactions` sub-cycles the reactive advance when the transport step exceeds
`CdrPlasmaJSON.chemistry_dt`. The plasma densities accumulate in place across substeps; the photon
production was **assigned** rather than accumulated, so only the final substep survived. The photon
production handed to the radiative transfer solver was therefore too small by a factor
`numSteps = ceil(dt/chemistry_dt)`.

Dormant at the default `chemistry_dt = 1.E99` (`numSteps == 1`), and it grows as `chemistry_dt` is lowered
to resolve the chemistry -- so lowering `chemistry_dt` improved the density integration and degraded the
photon production at the same time. At the settings used here (`chemistry_dt = 1.E-12`, dt ~ 5-9 ps) it was
a factor 5-9.

Fixed by changing `=` to `+=` at `CD_CdrPlasmaJSON.cpp:6081` (Euler), `:6126` (RK2) and `:6193` (RK4).
`photonProduction` is already zero-initialised by the caller at `:5235`, and all four call sites of the
three integrators are inside `integrateReactions`, so no other change is needed.

**This was the cause of the radius and head-field discrepancy.** With CdrPlasma emitting 5-9x too few
photons, it had far too little plasma around the streamer head, giving a narrower, more focused channel
and a higher head field (549 Td vs 410 Td for ItoKMC). Radii and velocities agree after the fix.

### Filed, not fixed

**`KMCSolver::advanceTau` never applies the tau-leaping condition** -- upstream issue
[#701](https://github.com/chombo-discharge/chombo-discharge/issues/701).

`advanceTau` bounds its substep only by the remaining splitting interval and refines only when the
resulting state comes out invalid. The leap condition `tau = eps*X_i/|sum_j nu_ij a_j|` is never applied,
even though the codebase already computes it in `getNonCriticalTimeStep` and `advanceHybrid` uses it
correctly. Consequences: the whole splitting interval is integrated in a single leap unless a population
goes negative, and **`ItoKMCJSON.prop_eps` has no effect** on the `explicit_euler`, `midpoint`, `prc` and
`implicit_euler` selectors. The `hybrid_*` selectors are unaffected. See the issue for the fix.

---

## 3. New feature: `parent` particle placement

`ItoKMCJSON`'s `"particle placement"` gained a `"parent"` method alongside `centroid`, `random` and
`downstream`. It draws one of the particles already in the cell, with probability proportional to its
weight, and puts the new particle on top of it.

The motivation is that `random` and `centroid` each inject a sub-grid transport error, with opposite signs:

* `random` scatters children uniformly over the cell, uncorrelated with their parents -- a variance
  injection of `dx^2/12` per dimension on the newly created fraction, i.e. an effective diffusivity
  `D_place = f*dx^2/(24*dt)`. At 7 AMR levels this was 1-3x the physical D at the streamer head. It
  broadens the front.
* `centroid` puts every child at the cell centre, dragging newly created charge backwards each step. It
  slows the front.

`parent` has neither. Implementation notes: parents are snapshotted before any append (otherwise the first
child becomes an eligible parent for its siblings); the cumulative weights are walked linearly to avoid a
per-cell allocation inside an OpenMP region; and an inherited position is by construction already on the
fluid side of the EB, so no cut-cell rejection sampling is needed.

**Known limitation.** The parent is drawn from the *same species* as the child. That is exactly right for
electrons -- the parent electron is the ionizer -- but for ions it is a proxy: a new O2+ inherits from an
existing O2+ in the cell rather than from the electron that created it. Photoionization products land in
cells with no particles of their species and fall back to `random`. Drawing the parent from the *reactant*
species would need the reaction's reactant list plumbed into `reconcileParticles`.

---

## 4. Why the diffusion coefficient is constant in both chemistry files

`ItoKMCPhysics::isotropicDiffusion` returns `sqrt(2*D*dt)*N(0,1)` and the stepper applies
`X += dt*v + hop`. That is Ito-interpretation Euler-Maruyama, whose Fokker-Planck equation is

    dn/dt = -div(v n) + grad^2(D n)      i.e. flux  v*n - grad(D*n)

whereas `CdrSolver` integrates

    dn/dt = -div(v n) + div(D grad n)    i.e. flux  v*n - D*grad(n)

These differ by `div(n grad D)` and **coincide exactly for constant D**. To get the fluid form out of an Ito
SDE the drift must be `v + grad(D)` (equivalently the Hanggi-Klimontovich / anti-Ito interpretation); there
is no grad-D term anywhere in `Source/ItoDiffusion` or `Physics/ItoKMC`.

This is *not* established as a bug. Which form is physically right when D varies over the front is beyond
the local-field approximation both models rest on, and the density-gradient expansion that defines D assumes
a locally uniform field, where the two agree. What is true is that the two solvers integrate different PDEs,
undocumented, and that the fluid convention is `D grad n`. Freezing D removes the ambiguity so the
comparison measures the models rather than their diffusion conventions. In air at 1 atm the bias was
|grad D|/|v| ~ 1-3% at the head, signed so the Ito swarm lags.

---

## 5. Numerical settings, and why they are what they are

Both models are numerically diffusive at the front, in different amounts and with different convergence
rates, and the settings in `example.inputs` were chosen to push both terms below the physical diffusion:

| term | source | scaling | at dx = 4.883e-6 |
|---|---|---|---|
| CdrPlasma | MUSCL/minmod front numerical diffusion, `0.5*v*dx*(1-CFL)` | `~dx` | 0.126 m^2/s at cfl 0.9 |
| ItoKMC | uniform secondary placement (+ merging) | `~dx^2` | 0.10 m^2/s |
| physical | | -- | 0.137 m^2/s |

Two consequences that are easy to get backwards:

* **`cfl` is not an accuracy dial in the usual direction.** The `(1-CFL)` factor is the O(dt) temporal error
  of the explicit step cancelling part of the O(dx) upwind diffusion. *Lowering* the CFL number *increases*
  the net numerical diffusion, saturating at `0.5*v*dx` as dt -> 0. `cfl = 0.9` is deliberate.
* **A dt-refinement study is not a convergence study for streamer velocity.** It converges, but to a
  semi-discrete answer that still carries `0.5*v*dx`. Only refining dx removes it. This is why the example
  runs at 8 AMR levels.

`CdrPlasmaJSON.chemistry_dt = 1.E-12` sub-cycles the CDR chemistry (~5 substeps) so it is not integrating
the reactions in one coarse leap. Note this is exactly the setting that exposed issue #702.

---

## 6. Audited and found equivalent -- do not re-derive

The two chemistry files and the two models were compared line by line. All of the following were verified
equivalent and are *not* sources of discrepancy:

* **Chemistry inputs.** Gas density, all reaction rates, the 1e-6 scale on the three-body attachment
  (correct, given the BOLSIG note that the 3-body cross section is normalised to gas density in cm^-3),
  photon generation efficiency and quenching, photoionization efficiency, initial seeds.
* **Quenching pressure.** Both parsers form `pq/(pq+p)` with `m_gasPressure` in **Pa** -- `CdrPlasmaJSON`
  converts its atm input at `CD_CdrPlasmaJSON.cpp:420`. This is the trap the differing input units set up,
  and it is handled correctly.
* **Townsend coefficients.** `"auto"` (ItoKMC, derived from the reaction list) vs tabulated (CdrPlasma)
  agree to 0.1% above ~130 Td, which is where the `refine_alpha` threshold bites. The meshes are not
  systematically different.
* **Field coupling.** rho-dagger sign and completeness (immobile ions are included, gated only on `Z != 0`);
  conductivity `Qe*sum|Z|*n*mu` using `abs(Z)` on both sides so species add rather than cancel; the
  semi-implicit operator `eps_r += dt*sigma/eps0`; identical cell-to-face averaging; same operator splitting
  order (transport+field -> reactions -> RTE).
* **Photon model.** Emission position (uniform in cell), direction, `partitionParticleWeights(num, 128)`,
  the `validCells` guard, `randomExponential(kappa)` transport, CIC deposition, and the photon count itself
  (both Poisson with the same mean, both via `Random::getPoisson`). Volume conventions agree --
  `pow(dx, SpaceDim)` in `McPhoto::computeNumPhysicalPhotons`, in `CdrPlasmaJSON`'s discrete-photon sampling
  and in `ItoKMCStepper`'s per-cell counts.
* **Absorption coefficients.** Sampled directly from both `RtSpecies` objects, 2e6 draws each: means differ
  by 0.52 sigma and agree with the analytic log-uniform mean; the endpoints `chi_min*p` and `chi_max*p`
  agree to 1 part in 1e6. The distributions are the same distribution.

---

## 7. Latent issues found but not fixed

None of these affect the current runs, but all are real:

* **Photoionization `"efficiency"` has different meanings in the two parsers.** In `ItoKMCJSON` it is a
  *relative branching ratio* among competing photo-reactions (`m_photoPathways` normalises by the sum), so
  with a single reaction it has no effect at all. In `CdrPlasmaJSON` it is an *absolute* multiplier. Porting
  a chemistry file between the models will silently change the answer. Both default to 1.0 here.
* **`ItoKMCJSON` silently ignores `"scale"` on a photoionization reaction**; `CdrPlasmaJSON` honours it.
* **Equal-vs-proportional weight partition.** With `ItoKMCJSON.increment_weights = true`,
  `reconcileParticles` distributes new physical particles *equally* across the computational particles in a
  cell instead of in proportion to their weights, transferring growth from heavy particles to light ones.
  Additionally, when `diff <= numParticles` only the first `diff` particles in SoA order are incremented.
  Dormant under equal-weight mergers such as `reinitialize`, where the equal split *is* the proportional one.
* **Hardcoded 3D volume in a 2D-capable path.** `CD_ItoKMCJSON.cpp:3072/3084` uses `std::pow(dx, 3)` in the
  `"ppc threshold"` grid factor. No reaction in this example uses that field.
* **`CdrPlasmaGodunovStepper` drops the diffusion term from the semi-implicit rho-dagger.** The
  `dt*sum(Z*q*div(D grad n))` contribution is `#if 0`-ed at `CD_CdrPlasmaGodunovStepper.cpp:1381-1389`, disabled
  in `dc9bfa73` (#208, 2022-03-21) -- ten months before the ItoKMC model existed. The adjacent comment's
  minus sign is also wrong; plus is correct, and matches what `ItoKMCGodunovStepper::computeDiffusionTermCDR`
  does.
* **`llround` quantises the KMC state.** `CD_ItoKMCStepperImplem.H:4205` rounds the per-cell count to an
  integer and discards the fractional part on write-back, even though `FPR = Real`. For *CDR* species inside
  `ItoKMC` this quantises the density to multiples of `1/dx^SpaceDim` and zeroes anything below half that.
  Irrelevant for Ito species, whose per-cell weights are large.
* **Dead code and dead options.** `ItoSolver::randomGaussian()` (and with it `ItoSolver.normal_max`) is never
  called -- the live hop uses the untruncated `Random::getNormal01()` in `ItoKMCPhysics::isotropicDiffusion`,
  and the `.options` comment misdescribes it as an exponential distribution. `McPhoto.random_kappa` is never
  read. `ItoKMCStepper::reconcilePhotoionization()` (the standalone overload at `:4664`) is never called; if
  it were wired up alongside the per-cell calls it would double-count photoionization.
* **`"forward isotropic"` diffusion is not a Brownian increment.** `forwardIsotropicDiffusion` reflects the
  component antiparallel to v, giving a non-zero mean displacement `sqrt(D/(pi*dt))` that *diverges* as
  dt -> 0. It is not a consistent discretization of anything. Not used here.

---

## 8. Investigated and ruled out

For the record, so these are not re-run: superparticle noise (estimated ~5% of physical D, and it scales as
`1/N_p`, which was tested); particle merging algorithm; secondary particle placement; deposition scheme
(both `cic`); time step (the Ito run is limited by `'Physics'` at ~3.7 ps, *below* the CDR CFL step);
double-counted photoionization; initial particle seeds; and the KMC midpoint propagator itself, which is a
correct Gillespie midpoint leap whose deterministic limit reproduces `1 + u + u^2/2`.
