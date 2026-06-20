# How particles are used in chombo-discharge

A survey of `Source/Particle` and every `ParticleContainer` user, grouped into
functional categories. Each category lists the operations and the *unique
requirements* it places on the new SoA core.

Particle types in the wild today (all derive from `GenericParticle<M,N>`, an AoS
struct with `position`, `particleID`, `rankID`, `M` reals, `N` RealVects):

| Type | Layout | Fields |
|------|--------|--------|
| `PointParticle` | `GenericParticle<1,0>` | weight |
| `Photon` | `GenericParticle<2,1>` | weight, kappa, velocity |
| `ItoParticle` | `GenericParticle<5,3>` | weight, mobility, diffusion, energy, tmp + oldPos, velocity, tmp |
| `TracerParticle<M,N>` | `GenericParticle<M,N>` | weight + user scratch |

Main users: `ItoSolver`/`ItoLayout`, `McPhoto`, `TracerParticleSolver`,
`BrownianWalker`, `DischargeInception`, `ItoKMC*`, plus `AmrMesh` (deposition glue).

---

## 1. Storage & AMR bookkeeping (`ParticleContainer`)

Operations: `define`/`regrid`/`preRegrid`, per-level `ParticleData<P>` holding
`List<P>` per box; buffer particles (grown grids), mask particles (halo near
refinement boundaries), cache particles (regrid), outcast particles. `remap()`,
`transferParticles`, `addParticles` (from `List`, `BinFab`, or another container).

Unique requirements:
- The SoA must be the per-box leaf storage that replaces `List<P>` inside the
  Chombo `ParticleData`/`BinFab` hierarchy — or we replace that hierarchy too.
- Cheap **move/transfer** between containers (today: list splice). SoA: move the
  column vectors.
- Multiple parallel holders of the *same* particle type (valid / buffer / mask /
  cache / outcast).

## 2. Particle ↔ mesh transfer (`EBParticleMesh`, `EBAMRParticleMesh`, deposition)

Operations: NGP/CIC/TSC **deposit** of a particle scalar/vector field onto an
`EBCellFAB`; **interpolate** a mesh field onto a particle field; coarse-fine
deposition, non-conservative + hybrid deposition in cut cells.

Current idiom (the crux):
```cpp
meshInterp.deposit<P, const Real&, &P::weight>(rho, particles, NGP, 1.0);
meshInterp.interpolate<P, RealVect&, &P::velocity>(particles, E, CIC);
```
i.e. a **pointer-to-member-function** selects *which field* to deposit/interp.

Unique requirements:
- Hot inner loop over all particles in a patch reading `position()` + one field.
  This is the loop that most wants SoA contiguity / vectorization.
- Need a compile-time way to name "which column" replacing `&P::weight`
  (prototype: column index/tag).
- Deposit reads one field; interpolate writes one field in place.

## 3. Per-cell sorting & kinetic Monte Carlo (`organizeParticlesByCell`, ItoKMC)

Operations: sort patch particles into a `BinFab<P>` (per-cell lists), run KMC
reactions per cell, merge/split, then re-flatten to patch storage.

Unique requirements:
- Fast **bucket by cell index** within a patch. SoA: either an index permutation
  per cell, or contiguous ranges + offsets (CSR-style), rather than physically
  moving particles per cell.
- Add/remove particles per cell during reactions (sources & sinks).

## 4. Merge / split / population control (`ParticleManagement`, `KDNode`)

Operations: KD-tree equal-weight partition & split
(`partitionAndSplitEqualWeightKD<P, &P::weight, &P::position>`),
`removePhysicalParticles`, `deleteParticles` (weight threshold),
`partitionParticleWeights`, with a `BinaryParticleReconcile<P>` hook to fix up
child-particle fields after a split.

Unique requirements:
- Frequent **insertion and deletion** of particles mid-array (swap-and-pop is fine
  since order is not semantically meaningful, except for determinism — see §7).
- Splitting copies a particle then tweaks a couple of fields → needs gather/scatter
  of a single particle and a per-field "reconcile" callback.
- KD-tree needs `weight` (mutable) and `position` (const) — again the "named field"
  problem.

## 5. Sampling / initialization (`ParticleManagement::draw*`)

Operations: draw N particles from sphere/box/Gaussian/custom distributions,
partitioned across MPI ranks.

Unique requirements: bulk append of freshly-constructed particles; only
`position` + `weight` are mandatory at creation.

## 6. Reductions & queries (`ParticleOps`, container counters)

Operations: `sum<P, &P::weight>`, particles-per-cell (physical & computational),
local/global counts (valid/outcast/mask), `removeParticles(functor)`,
`transferParticles(functor)`, `setValue<P, &P::field>`, `setData(functor)`.

Unique requirements:
- Whole-column reductions (SoA shines: `sum(weights())`).
- Predicate-based remove/transfer (functor over a particle view).
- Set-a-whole-column-to-a-value.
- Hybrid OpenMP reductions over particles (current code declares custom OMP
  reductions over `List<P>`).

## 7. Communication & I/O (linearization)

Operations:
- **MPI**: `linearOut`/`linearIn` pack the *full* particle (incl. `particleID`,
  `rankID`) to/from a byte buffer; `ParticleOps::scatterParticles` ships particles
  between ranks during remap.
- **HDF5**: `H5size`/`H5linearOut`/`H5linearIn` pack only the `Real` components
  (no IDs) so they plug into Chombo's particle checkpoint API.
- **Determinism**: `operator<` gives a lexicographic ordering; `sortParticles()`
  and `resetParticleIDs()` are used for reproducible output.

Unique requirements:
- **Two** linearization flavors: full (MPI) and Real-subset (HDF5).
- Must be cheap and ideally **generic** — see the dedicated note in
  `LINEARIZATION.md`.
- Deterministic ordering / stable IDs for checkpoint-restart and regression.

---

## 8. ItoDiffusion + ItoKMC kinetics (the most demanding consumer)

Surveyed `Source/ItoDiffusion/CD_ItoSolver.{cpp,Impl}` and
`Physics/ItoKMC/CD_ItoKMCStepper*`. These drive most of the per-cell particle churn
and impose the tightest requirements.

- **Reaction reconciliation** (`reconcileParticles`): after the KMC chemistry step
  computes new physical-particles-per-cell per species, each cell must *create* or
  *destroy* particles to match. Done **per cell** on cell-sorted `BinFab<P>` storage,
  for **many containers at once** (one `BinFab<ItoParticle>*` per Ito species, one
  `BinFab<PointParticle>*` per CDR photoionization product, photons per species).
  New particles are sampled at random positions **inside the cell, cut-cell-aware**
  (uses `EBISBox`/`CellInfo`/valid-cell mask). Runs under `#pragma omp parallel for`
  over patches.
  → Requirements: cheap per-cell append/remove; bulk creation; lockstep iteration
    over several SoA containers keyed by the same cell; cut-cell geometry per cell.

- **Superparticle merge/split** (`ItoSolver::makeSuperparticles`): per cell, reduce
  to a target particle count `m_particlesPerCell[level]`. Two strategies:
  equal-weight **KD-tree** partition/split and **BVH reinitialize**. Uses `CellInfo`
  and a per-split **reconcile** callback to fix child fields. Invoked every
  `m_mergeInterval` steps and on regrid.
  → Requirements: gather a cell's particles into a working set, run a tree algorithm
    over `(position, weight)`, write back a different number of particles; per-field
    reconcile hook on split.

- **Covered-particle removal** (`removeCoveredParticles`): remove particles inside
  the EB / covered region with an `EBRepresentation` + tolerance, across many species.
  → Requirements: geometry-predicate removal (swap-and-pop), applied per container.

- **Secondary emission**: `m_secondaryParticles` (ItoParticle) and
  `m_secondaryPhotons` (Photon) are *transient side-containers* filled at the EB,
  then merged/transferred into the bulk containers and cleared.
  → Requirements: many short-lived containers of the same type; cheap
    transfer/merge/clear.

- **Per-cell accounting**: physical PPC = Σ weights per cell (reactive vs. total);
  computational PPC = count per cell; mean energy per cell. Mixed regular/cut-cell
  loops.
  → Requirements: per-cell reductions over a chosen scalar column with cut-cell
    awareness.

- **Velocity/mobility interpolation then in-place rescale** (`ItoSolver`): interpolate
  a mesh field into `velocity`, then `velocity *= mobility` per particle — an
  in-place update touching two columns. (Captured in the prototype's
  `testInterpolateLoop`.)

These confirm that SoA storage must support: (i) a **cell-sorted view** that many
containers share, (ii) fast per-cell add/remove and variable output counts,
(iii) a per-particle "working set" gather for tree-based merge/split, and (iv) a
per-split field-reconcile callback.

---

## Cross-cutting "named field" problem

Categories 2, 4, 6 all rely today on `&P::someField` (pointer-to-member or
pointer-to-member-function) as a *compile-time* selector of which quantity to
operate on. SoA removes member functions, so the new core needs an equally cheap
compile-time field selector. The prototype answers this with **column indices**
exposed by `ParticleTraits<P>` (e.g. `Traits::velocityIndex`); a tag-type scheme is
an alternative. This is the single most important API decision.
