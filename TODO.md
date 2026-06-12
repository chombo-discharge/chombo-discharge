# DataOps Reorganization — Checkpoint List

## Checkpoint 1 — Literalincludes with line subsets

Files referenced via `.. literalinclude::` **with** `:lines:` or `:emphasize-lines:` in the
Sphinx documentation. These must be audited and updated whenever the source lines shift due to
reorganization, before merging a PR.

Each entry lists the source file and a brief note on what code segment is pulled in.

---

### Physics/AdvectionDiffusion/

**`CD_AdvectionDiffusionStepper.H`** — *Applications/AdvectionDiffusionModel.rst*
- `:lines: 215-220` — general initial data specification
- `:lines: 222-227` — user-specified velocity field
- `:lines: 229-234` — diffusion coefficient setting

**`CD_AdvectionDiffusionStepper.cpp`** — *Source/TimeStepper.rst*
- `:lines: 136-162` — `setupSolvers` implementation example
- `:lines: 175-185` — `registerOperators` implementation example
- `:lines: 187-196` — `allocate` implementation example
- `:lines: 198-218` — `initialData` implementation example
- `:lines: 221-230` — `writeCheckpointData` implementation example
- `:lines: 234-243` — `readCheckpointData` implementation example
- `:lines: 265-275` — `getNumberOfPlotVariables` implementation example
- `:lines: 277-286` — `getPlotVariableNames` implementation example
- `:lines: 288-303` — `writePlotData` implementation example
- `:lines: 305-350` — `computeDt` implementation example
- `:lines: 352-454` — `advance` implementation example
- `:lines: 456-469` — `synchronizeSolverTimes` implementation example
- `:lines: 471-480` — `preRegrid` implementation example
- `:lines: 482-501` — `regrid` implementation example

**`CD_AdvectionDiffusionStepper.options`** — *Applications/AdvectionDiffusionModel.rst*
- `:lines: 11-13` — basic input options
- `:lines: 14` — velocity field adjustment flag
- `:lines: 15` — omega parameter
- `:lines: 23-27` — general runtime adjustments

---

### Physics/DischargeInception/

**`CD_DischargeInceptionStepper.H`** — *Applications/DischargeInceptionModel.rst*
- `:lines: 73-93` — model implementation overview
- `:lines: 310-407` — user API functions

---

### Physics/Electrostatics/

**`CD_FieldStepper.H`** — *Applications/ElectrostaticsModel.rst*
- `:lines: 236-242` — space charge specification function
- `:lines: 244-249` — surface charge specification function

---

### Physics/ItoKMC/

**`CD_ItoKMCGodunovStepper.options`** — *Applications/ItoKMC.rst*
- `:lines: 22-35` — input variable limits section

---

### Physics/MeshODE/

**`CD_MeshODEStepper.H`** — *Applications/MeshODEModel.rst*
- `:lines: 37-48` — class template specification

---

### Source/AmrMesh/

**`CD_AmrMesh.H`** — *Source/MeshData.rst*, *Source/Particles.rst*
- `:lines: 98-116` — general allocation (MeshData.rst)
- `:lines: 447-460` — gradient computation (MeshData.rst)
- `:lines: 701-708,733-741,768-780` — multiple allocation signatures (MeshData.rst)
- `:lines: 1175-1182` — ghost cell update signatures (MeshData.rst)
- `:lines: 1184-1197` — single grid level ghost cell update (MeshData.rst)
- `:lines: 1264-1283` — fine-grid interpolation (MeshData.rst)
- `:lines: 221-228` — `ParticleContainer` allocation (Particles.rst)
- `:lines: 945-981` — particle deposition (Particles.rst)
- `:lines: 983-1001` — mesh data interpolation to particles (Particles.rst)
- `:lines: 1115-1173` — particle intersection algorithm (Particles.rst)

**`CD_AmrMesh.options`** — *Source/AmrMesh.rst*
- `:emphasize-lines: 9-15,21-26` — runtime-adjustable options highlighted in template

---

### Source/CellTagger/

**`CD_CellTagger.H`** — *Source/CellTagger.rst*
- `:lines: 50-63,71-81` + `:emphasize-lines: 6-7,13-14,24-25` — header excerpt with key members highlighted
- `:lines: 89-111` + `:emphasize-lines: 5-6,12-13,22-23` — required function implementations highlighted

---

### Source/ConvectionDiffusionReaction/

**`CD_CdrMultigrid.H`** — *Solvers/CDR.rst*
- `:lines: 289-296` — pure function implementations

**`CD_CdrSolver.H`** — *Solvers/CDR.rst*
- `:lines: 122-135,146-158` — implicit diffusion advance methods
- `:lines: 160-210` — divergence operator approximations
- `:lines: 212-220` — function call requirement
- `:lines: 600-668` — mesh data retrieval functions

---

### Source/Driver/

**`CD_Driver.H`** — *Source/Driver.rst*
- `:lines: 40-50` — constructor signature

**`CD_Driver.options`** — *Source/Driver.rst*
- `:emphasize-lines: 4,8-13,17-18,20-21,24-27,29-33` — runtime-adjustable options highlighted

**`CD_TimeStepper.H`** — *Source/TimeStepper.rst*
- `:lines: 149-157` — function signature
- `:lines: 161-168` — level-by-level data reading function
- `:lines: 186-194` — function signature
- `:lines: 234-241` — function signature
- `:lines: 294-301` — function signature
- `:lines: 303-309` — function signature
- `:lines: 311-330` — load balance realm function

**`CD_TimeStepper.cpp`** — *Source/TimeStepper.rst*
- `:lines: 107-137` — default `loadBalanceBoxes` implementation

---

### Source/Electrostatics/

**`CD_FieldSolver.H`** — *Solvers/Electrostatics.rst*
- `:lines: 112-122` — pure member function
- `:lines: 386-392` — voltage setting function
- `:lines: 406-412` — member function

**`CD_FieldSolverGMG.options`** — *Solvers/Electrostatics.rst*
- `:emphasize-lines: 4,12-14,16-37` — runtime-adjustable options highlighted

---

### Source/Geometry/

**`CD_ComputationalGeometry.H`** — *Source/ComputationalGeometry.rst*
- `:lines: 49-61` — electrode and dielectric retrieval functions
- `:lines: 163-201` + `:emphasize-lines: 4,24,29` — data members with mandatory members highlighted

**`CD_Dielectric.H`** — *Source/ComputationalGeometry.rst*
- `:lines: 37-51` — `Dielectric` constructors

**`CD_Electrode.H`** — *Source/ComputationalGeometry.rst*
- `:lines: 34-41` — `Electrode` constructor

---

### Source/ItoDiffusion/

**`CD_ItoParticle.H`** — *Solvers/Ito.rst*
- `:lines: 40` — `GenericParticle` derivation declaration

**`CD_ItoSolver.H`** — *Solvers/Ito.rst*
- `:lines: 92-97` — particle merging function
- `:lines: 204-217,228-240` — deposited quantity functionality
- `:lines: 293-299` — bulk particle density deposition
- `:lines: 333-347` — general deposition method
- `:lines: 414-429` — relevant function
- `:lines: 631-637` — particle retrieval
- `:lines: 654-659` — data fetching
- `:lines: 731-737,749-755` — particle velocity interpolation
- `:lines: 763-768` — function signature
- `:lines: 791-797` — particle splitting and merging
- `:lines: 864-875` — particle remapping functions
- `:lines: 899-909` — combined advection and diffusion
- `:lines: 986-991` — routine implementation
- `:lines: 1010-1015` — diffusion time step signatures

**`CD_ItoSpecies.H`** — *Solvers/Ito.rst*
- `:lines: 36-43` — `ItoSpecies` constructor
- `:lines: 105-117` — construction-time or explicit population
- `:lines: 152-160` — mandatory data members

---

### Source/KineticMonteCarlo/

**`CD_KMCSingleState.H`** — *Solvers/KineticMonteCarlo.rst*
- `:lines: 95-100` — required member function

**`CD_KMCSingleStateReaction.H`** — *Solvers/KineticMonteCarlo.rst*
- `:lines: 75-113` — required member functions

**`CD_KMCSolver.H`** — *Solvers/KineticMonteCarlo.rst*
- `:lines: 42-63,79-84` — KMC solver implementation
- `:lines: 112-122` — hybrid solver parameter setting
- `:lines: 413-423` — hybrid advance method

---

### Source/MeshODESolver/

**`CD_MeshODESolver.H`** — *Solvers/MeshODESolver.rst*
- `:lines: 23-28` — class template
- `:lines: 31-34,121-127` — source term computation
- `:lines: 41-46` — full constructor
- `:lines: 99-112` — data setting via function
- `:lines: 114-120` — component-by-component source setting
- `:lines: 137-143` — mesh component fetching
- `:lines: 151-156` — source term fetching
- `:lines: 165-172` — regrid data preparation
- `:lines: 174-182` — regridded data retrieval

---

### Source/Particle/

**`CD_GenericParticle.H`** — *Source/Particles.rst*
- `:lines: 68-76` — particle type template
- `:lines: 165-173,184-193` — Real and RealVect variable fetching

**`CD_ParticleContainer.H`** — *Source/Particles.rst*
- `:lines: 166-172` — function signatures
- `:lines: 216-220` — particle deletion
- `:lines: 284-289` — particle masking and copying
- `:lines: 459-463` — remapping function

**`CD_ParticleManagement.H`** — *Solvers/Ito.rst*
- `:lines: 34-42` — `ParticleMerger` alias

---

### Source/RadiativeTransfer/

**`CD_EddingtonSP1.H`** — *Solvers/RTE.rst*
- `:lines: 81-91` — solution advance member function
- `:lines: 117-127` — member function specification

**`CD_EddingtonSP1.options`** — *Solvers/RTE.rst*
- `:emphasize-lines: 4,6-9,20-33` — runtime-adjustable options highlighted

**`CD_EddingtonSP1DomainBc.H`** — *Solvers/RTE.rst*
- `:lines: 44-47` — `a_function` argument definition

**`CD_McPhoto.H`** — *Solvers/RTE.rst*
- `:lines: 408-441` — particle retrieval functions

**`CD_RtSolver.H`** — *Solvers/RTE.rst*
- `:lines: 260-279` — variable setting functions

**`CD_RtSpecies.H`** — *Solvers/RTE.rst*
- `:lines: 50-56` — lightweight class description

---

### Source/SurfaceODESolver/

**`CD_SurfaceODESolver.H`** — *Solvers/SurfaceODESolver.rst*
- `:lines: 23-28` — class template
- `:lines: 33-44` — solver constructors
- `:lines: 191-203` — direct data setting
- `:lines: 213-218` — mesh data fetching
- `:lines: 227-254` — right-hand side setting
- `:lines: 281-289` — regrid preparation
- `:lines: 291-299` — regridded data retrieval
- `:lines: 301-337` — function signatures

---

### Source/TracerParticles/

**`CD_TracerParticle.H`** — *Solvers/TracerParticles.rst*
- `:lines: 26-33` — `TracerParticle` inheritance
- `:lines: 52-78` — data member accessibility

**`CD_TracerParticleSolver.H`** — *Solvers/TracerParticles.rst*
- `:lines: 23-37` — solver template specification
- `:lines: 57-62` — full constructor initialization
- `:lines: 147-152` — velocity field setting
- `:lines: 186-191` — particle deposition
- `:lines: 193-198` — scalar field interpolation
- `:lines: 200-204` — particle velocity computation
- `:lines: 281-293` — solver particle retrieval

---

### Source/Utilities/

**`CD_DataOps.H`** — *Solvers/CDR.rst*
- `:lines: 1176-1188` — `DataOps` lambda functions used in CDR solver example

**`CD_DischargeIO.H`** — *Source/Particles.rst*
- `:lines: 141-163` — `writeH5Part` function signature

**`CD_LookupTable1D.H`** — *Utilities/LookupTable.rst*
- `:lines: 30-31` — class template declaration
- `:lines: 62-68` — `addData` member function
- `:lines: 77-85` — column swapping function
- `:lines: 87-93` — column scaling function
- `:lines: 95-104` — data range truncation
- `:lines: 106-118` — range strategy setting
- `:lines: 120-127` — table regularization
- `:lines: 129-144` — data interpolation
- `:lines: 174-200` — debugging output functions

**`CD_Random.H`** — *Utilities/RandomNumbers.rst*
- `:lines: 72-118` — pre-defined distributions
- `:lines: 120-127` — random number drawing routine

---

*Total: 44 source files, ~139 line-subset literalinclude directives across 25 RST files.*

---

## Checkpoint 2 — AmrMesh extensions for pre-built VoFIterators

### Background

`DataOps` currently constructs `VoFIterator` objects on-the-fly inside every function call by
querying `EBISBox::getIrregIVS(box)` and building a new iterator from the result. This incurs
allocation and traversal cost at every call site. The fix is to move VoFIterator ownership into
`AmrMesh` so that iterators are built once (at regrid time) and handed to callers via
`getVofIterator(realm, phase)`.

`AmrMesh` already stores a VoFIterator per realm+phase that covers **all cut-cells** (irregular
cells, whether singly or multiply cut), exposed via:

```cpp
Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
getVofIterator(const std::string& a_realm, const phase::which_phase a_phase) const;
```

Two things are missing:

1. A VoFIterator that covers only **multi-valued** (multiply-cut) cells — those cells for which
   `EBISBox::getMultiCells(box)` is non-empty and which require a separate pass when the
   `getSingleValuedFAB` path already covers the singly-cut VoF.

2. Both iterators must survive regrid, i.e., be rebuilt in `regridOperators`.

### Required changes to `CD_AmrMesh.H` / `CD_AmrMesh.cpp` / `CD_Realm.H` / `CD_Realm.cpp`

#### 2a. New public accessor on `AmrMesh`

Add a new getter parallel to `getVofIterator`:

```cpp
/** @brief Get VoFIterators that iterate only over multiply-cut cells (cells with
 *         more than one VoF).
 *  @param[in] a_realm Realm name
 *  @param[in] a_phase Phase (gas or solid)
 *  @return Per-level LayoutData of VoFIterators restricted to multi-valued cells.
 */
Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
getMultiCutVofIterator(const std::string& a_realm, const phase::which_phase a_phase) const;
```

Return type and calling convention must be identical to `getVofIterator` so callers can be
written symmetrically.

#### 2b. Storage in `Realm`

`Realm` already stores the all-cut-cell VoFIterators in a private member (pattern:
`m_vofIter[phase][level]` or similar). Add a parallel member for multi-cut cells:

```
m_multiCutVofIter  —  same type as existing m_vofIter
```

#### 2c. Population during `regridOperators`

In the regrid path where the all-cut iterator is built using
`ebisbox.getIrregIVS(validBox)`, add a second construction using
`ebisbox.getMultiCells(validBox)`:

```cpp
// Existing (all cut-cells):
VoFIterator& allCutIter = (*m_vofIter[phase][lvl])[din];
allCutIter.define(ebisbox.getIrregIVS(validBox), ebgraph);

// New (multiply-cut cells only):
VoFIterator& multiCutIter = (*m_multiCutVofIter[phase][lvl])[din];
multiCutIter.define(ebisbox.getMultiCells(validBox), ebgraph);
```

`EBISBox::getMultiCells(box)` returns an `IntVectSet` containing only cells that have more
than one VoF (multi-valued cells). This is the correct iteration space for the second pass.

#### 2d. Pre-regrid cleanup

In the `preRegrid` path, clear `m_multiCutVofIter` alongside the existing VoFIterator cleanup.

#### Summary of files to modify

The actual VoFIterator storage lives in `PhaseRealm` (`m_vofIter` at line 504 of
`CD_PhaseRealm.H`). `Realm` is a thin dispatcher that keys into
`m_realms[phase]->getVofIterator()`, and `AmrMesh` dispatches by realm name into a `Realm`.
The delegation chain is:

```
AmrMesh::getVofIterator(realm, phase)
  → Realm::getVofIterator(phase)
    → m_realms[phase]->getVofIterator()   // PhaseRealm
      → PhaseRealm::m_vofIter             // actual storage
```

The new `getMultiCutVofIterator` must follow the same chain:

| File | Change |
|------|--------|
| `Source/AmrMesh/CD_PhaseRealm.H` | Add `m_multiCutVofIter` data member; add `getMultiCutVofIterator()` declaration |
| `Source/AmrMesh/CD_PhaseRealm.cpp` | Populate `m_multiCutVofIter` in `regridOperators` using `getMultiCells`; clear in `preRegrid` |
| `Source/AmrMesh/CD_Realm.H` | Add `getMultiCutVofIterator(phase)` declaration |
| `Source/AmrMesh/CD_Realm.cpp` | Implement as `return m_realms[a_phase]->getMultiCutVofIterator();` |
| `Source/AmrMesh/CD_AmrMesh.H` | Add `getMultiCutVofIterator(realm, phase)` declaration |
| `Source/AmrMesh/CD_AmrMesh.cpp` | Implement; delegate to `Realm` (same pattern as `getVofIterator`) |

---

## Checkpoint 3 — DataOps signature changes

### Design principle

`EBAMRCellData` carries a realm name but **no** phase. The phase is known only to the caller
(solver or time stepper). Therefore, any `DataOps` function that needs pre-built VoFIterators
must receive them from the caller rather than fetching them itself.

Preferred signature convention: **pass the iterator directly**, not `AmrMesh + realm + phase`.
This keeps `DataOps` independent of `AmrMesh` and makes the dependency explicit at the call
site. The caller does:

```cpp
auto& vofIter = m_amr->getVofIterator(m_realm, m_phase);
DataOps::someFunction(data, ..., vofIter);
```

#### Functions that need NO change (already use Chombo built-in EBCellFAB arithmetic)

These delegate entirely to `EBCellFAB`'s own operators, which handle all VoFs internally.
No `VoFIterator` argument needed:

| Function | Mechanism |
|----------|-----------|
| `incr(LD<EBCellFAB>, LD<EBCellFAB>, Real)` | `EBCellFAB::plus()` |
| `scale(LD<EBCellFAB>, Real)` | `EBCellFAB::mult()` |
| `multiply(LD<EBCellFAB>, LD<EBCellFAB>)` | `EBCellFAB::operator*=` |
| `multiplyScalar(LD<EBCellFAB>, LD<EBCellFAB>)` | `EBCellFAB::mult()` |
| `axby(LD<EBCellFAB>, ...)` | `EBCellFAB::axby()` |
| `divide(LD<EBCellFAB>, LD<EBCellFAB>, int, int)` | `EBCellFAB::divide()` |
| `divideByScalar(LD<EBCellFAB>, LD<EBCellFAB>)` | calls `divide()` |
| `plus(LD<EBCellFAB>, LD<EBCellFAB>, ...)` | `EBCellFAB::plus()` |
| `copy(EBAMRCellData/LD<EBCellFAB>, ...)` | `localCopyTo` |
| `setValue(LD<EBCellFAB>, Real)` / `setValue(LD<EBCellFAB>, Real, int)` | `EBCellFAB::setVal()` |

#### Functions that need NO change because VoFIterators are self-contained in the data holder

These operate on `BaseIVFAB<Real>` whose `getIVS()` / `getEBGraph()` are embedded in the
data holder itself — no realm or phase information is needed:

| Function |
|----------|
| `incr(LD<BaseIVFAB>, LD<BaseIVFAB>, Real)` |
| `incr(LD<EBCellFAB>, LD<BaseIVFAB>, Real)` |
| `incr(LD<BaseIVFAB>, LD<EBCellFAB>, Real)` |
| `floor(LD<BaseIVFAB<Real>>, Real)` |
| `roof(LD<BaseIVFAB<Real>>, Real)` |
| `multiply(LD<BaseIVFAB>, LD<BaseIVFAB>)` |
| `multiplyScalar(LD<BaseIVFAB>, LD<BaseIVFAB>)` |
| `scale(LD<BaseIVFAB<Real>>, Real)` |
| `getMaxMinNorm(Real&, Real&, LD<BaseIVFAB<Real>>)` |
| `setValue(LD<BaseIVFAB<Real>>, function, ...)` |

#### Functions requiring a new `phase::which_phase` parameter (or pre-built VoFIterator)

These currently call `ebisbox.getIrregIVS(box)` to build their own `VoFIterator`. They must
instead accept a `const LayoutData<VoFIterator>&` (level-level) or
`const Vector<RefCountedPtr<LayoutData<VoFIterator>>>&` (AMR-level):

**Full cut-cell iterator required** (kernel is EB-geometry-sensitive, or regular kernel
explicitly skips cut-cells via `isRegular` / `volFrac`):

| Function | Reason |
|----------|--------|
| `averageFaceToCell(LD<EBCellFAB>, LD<EBFluxFAB>, ProblemDomain)` | Irregular kernel calls `ebisbox.getFaces(vof, dir, sit())` — geometry-specific |
| `kappaScale(LD<EBCellFAB>)` | Only has irregular kernel; uses `ebisbox.volFrac(vof)` |
| `kappaSum(Real&, LD<EBCellFAB>, int)` | Regular kernel guards `isRegular(iv)`, skipping all cut-cells; volFrac weighting in irregular kernel |
| `norm(LD<EBCellFAB>, int, int)` | Regular kernel guards `isRegular(iv)`, so all cut-cells must go through VoFIter |

**Multi-cut iterator sufficient** (kernel is EB-geometry-independent; `getSingleValuedFAB` +
BoxLoops::loop handles regular cells and singly-cut VoF 0; only the extra VoFs in
multiply-cut cells need a separate pass):

| Function | Notes |
|----------|-------|
| `floor(LD<EBCellFAB>, Real)` | Operates on full region including ghosts — AmrMesh iterator covers only valid box; ghost multi-cut cells must be handled separately or the full VoFIterator used |
| `roof(LD<EBCellFAB>, Real)` | Same note as `floor` |
| `max(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>)` | Geometry-independent; multi-cut VoFIter for extra VoFs |
| `compute(LD<EBCellFAB>, function)` | Same function applied to all cells; multi-cut VoFIter |
| `dotProduct(EBCellFAB, EBCellFAB, EBCellFAB, Box)` | Identical kernel for regular and cut; multi-cut VoFIter |
| `divideFallback(LD<EBCellFAB>, LD<EBCellFAB>, Real/LD)` | Clone prevents interference; same logic in both kernels; multi-cut VoFIter |
| `getMaxMin(LD<EBCellFAB>, int)` | No `isRegular` guard; regular kernel already covers VoF 0 of cut-cells; multi-cut VoFIter for extra VoFs only |
| `getMaxMinNorm(LD<EBCellFAB>)` | Covered-cell mask skips covered cells; regular kernel handles cut-cells; multi-cut VoFIter |
| `setInvalidValue(EBAMRCellData, Vector<int>, Real)` | Same value for all cells; multi-cut VoFIter |
| `setValue(LD<EBCellFAB>, function<Real(RealVect)>, ...)` | Same function applied everywhere; multi-cut VoFIter |
| `setValue(LD<EBCellFAB>, function<RealVect(RealVect)>, ...)` | Same |
| `vectorLength(EBCellFAB, EBCellFAB, Box)` | Geometry-independent length computation; multi-cut VoFIter |
| `vectorLength2(EBCellFAB, EBCellFAB, Box)` | Same |

**Special cases:**

- `filterSmooth(LD<EBCellFAB>, ...)` — builds a grown `IrregIVS` (irregCells grown by
  `a_stride - 1`), which cannot be pre-built in AmrMesh without knowing the stride. This
  function should **keep its local VoFIterator construction**; it is inherently stride-dependent.

- `MFCellFAB` variants (`setValue(LD<MFCellFAB>, ...)`, `kappaScale(LD<MFCellFAB>)`,
  `squareRoot(LD<MFCellFAB>)`) — these iterate internally over phases and build per-phase
  VoFIterators. Passing phase-keyed iterators from outside would require a vector-of-phases
  argument. Defer to a follow-on task; keep local VoFIterator construction for now.
  Note: `squareRoot(LD<MFCellFAB>)` already uses `getMultiCells(box)` internally — a sign
  that this design was partially anticipated.

#### Signature convention for changed functions

At the AMR level, add a trailing parameter:

```cpp
static void setValue(
    EBAMRCellData&                              a_lhs,
    const std::function<Real(const RealVect)>&  a_function,
    const RealVect&                             a_probLo,
    const Vector<Real>&                         a_dx,
    const int                                   a_comp,
    const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);
```

At the level, add a trailing parameter:

```cpp
static void setValue(
    LevelData<EBCellFAB>&                       a_lhs,
    const std::function<Real(const RealVect)>&  a_function,
    const RealVect                              a_probLo,
    const Real                                  a_dx,
    const int                                   a_comp,
    const LayoutData<VoFIterator>&              a_vofIter);
```

The iterator is used directly in `BoxLoops::loop((*a_vofIter)[din], kernel)` replacing the
on-the-fly `VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph)` construction.

---

## Checkpoint 4 — BoxLoops internals reorganization in DataOps

### Principle

In functions that use the pattern

```cpp
BoxLoops::loop<D_DECL(1,1,1)>(box, regularKernel);   // operates on getSingleValuedFAB
BoxLoops::loop(vofit, irregularKernel);               // operates on EBCellFAB via VolIndex
```

the `BoxLoops::loop<>` over the box already covers:
- all regular cells
- the singly-cut (single-valued) cut-cell VoF (VoF index 0) via `getSingleValuedFAB`

The `VoFIterator` is then redundant for singly-cut cells and only needed for the **extra VoFs**
in multiply-cut cells. When the kernel is EB-geometry-independent (same computation whether
accessed by IntVect or VolIndex), the `VoFIterator` can be replaced by a multi-cut iterator.

### Category 1 — BoxLoops::loop alone; no VoFIterator pass at all

These functions already use Chombo built-in `EBCellFAB` arithmetic which handles all VoFs
internally. The BoxLoops loop is not even used in these:

| Function | Why no irregular pass needed |
|----------|------------------------------|
| `incr(LD<EBCellFAB>, LD<EBCellFAB>)` | `EBCellFAB::plus()` walks all VoFs |
| `scale(LD<EBCellFAB>)` | `EBCellFAB::mult()` walks all VoFs |
| `multiply(LD<EBCellFAB>, LD<EBCellFAB>)` | `operator*=` walks all VoFs |
| `multiplyScalar(LD<EBCellFAB>, LD<EBCellFAB>)` | `EBCellFAB::mult()` walks all VoFs |
| `axby(LD<EBCellFAB>, ...)` | `EBCellFAB::axby()` walks all VoFs |
| `divide(LD<EBCellFAB>, LD<EBCellFAB>, int, int)` | `EBCellFAB::divide()` walks all VoFs |
| `plus(LD<EBCellFAB>, LD<EBCellFAB>, ...)` | `EBCellFAB::plus()` walks all VoFs |
| `setValue(LD<EBCellFAB>, Real)` | `EBCellFAB::setVal()` walks all VoFs |
| All `BaseIVFAB` variants | IVS-based VoFIterator is self-contained; single pass is correct by construction |

### Category 2 — BoxLoops::loop + multi-cut VoFIterator pass

The regular kernel via `getSingleValuedFAB` handles regular cells and VoF 0 of all cut-cells.
A second pass with the **multi-cut** VoFIterator (from `AmrMesh::getMultiCutVofIterator`)
handles the extra VoFs in multiply-cut cells. The two kernels are identical in logic.

| Function | Note |
|----------|------|
| `max(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>)` | Geometry-independent pointwise max |
| `compute(LD<EBCellFAB>, function)` | Same transform applied everywhere |
| `dotProduct(EBCellFAB, EBCellFAB, EBCellFAB, Box)` | Sum of products — same in all cells |
| `divideFallback(LD<EBCellFAB>, LD<EBCellFAB>, Real/LD)` | Clone pattern ensures correctness; kernels compute same value |
| `getMaxMin(LD<EBCellFAB>, int)` | No `isRegular` guard; BoxLoops sees VoF 0; multi-cut catches extra VoFs |
| `getMaxMinNorm(LD<EBCellFAB>)` | Covered-cell mask; cut-cells already processed by regular loop; multi-cut for extra VoFs |
| `setInvalidValue(EBAMRCellData, ...)` | Same constant written everywhere |
| `setValue(LD<EBCellFAB>, function<Real(RealVect)>, ...)` | Same spatial function applied everywhere |
| `setValue(LD<EBCellFAB>, function<RealVect(RealVect)>, ...)` | Same |
| `vectorLength(EBCellFAB, EBCellFAB, Box)` | EB-geometry-independent length formula |
| `vectorLength2(EBCellFAB, EBCellFAB, Box)` | Same |

**Caveat — ghost cells in `floor` and `roof`:** These two functions use `lhs.getRegion()`
(the full allocated region including ghost cells) rather than the valid box. AmrMesh
VoFIterators only cover the valid box. Options:
  - (a) Accept a full-cut VoFIterator for these two functions (simplest)
  - (b) Fall back to local VoFIterator construction over the full region (preserves current semantics, but defeats the purpose for ghost multi-cut cells)
  - (c) Assert/document that multi-cut cells in ghost regions are pathological and use the multi-cut iterator over the valid box only

Recommend option (a) for now: pass the **all-cut** VoFIterator from AmrMesh for `floor` and
`roof`, accepting that VoF 0 of singly-cut cells gets a redundant but harmless write.

### Category 3 — Full cut-cell VoFIterator required

The irregular kernel is fundamentally different from the regular kernel (uses EB geometry), OR
the regular kernel explicitly guards `isRegular(iv)` to skip all cut-cells. The pre-built
all-cut VoFIterator from `AmrMesh::getVofIterator` must be passed in.

| Function | Reason |
|----------|--------|
| `averageFaceToCell(LD<EBCellFAB>, LD<EBFluxFAB>, ...)` | Irregular kernel calls `ebisbox.getFaces(vof, dir, sit())` — face connectivity is geometry-specific and cannot be reproduced from cell-center data alone |
| `kappaScale(LD<EBCellFAB>)` | No regular kernel; only irregular kernel exists, multiplying by `ebisbox.volFrac(vof)` — intrinsically EB-specific |
| `kappaSum(Real&, LD<EBCellFAB>, int)` | Regular kernel guards `ebisbox.isRegular(iv)` → skips all cut-cells; irregular kernel provides volFrac-weighted sum over cut-cells only |
| `norm(LD<EBCellFAB>, int, int)` | Regular kernel guards `ebisbox.isRegular(iv)` → skips all cut-cells; irregular kernel handles them |

### Category 4 — Keep local VoFIterator construction (cannot pre-build)

| Function | Reason |
|----------|--------|
| `filterSmooth(LD<EBCellFAB>, Real, int, bool)` | Builds `IrregIVS` grown by `a_stride - 1`; the iteration space is stride-dependent and cannot be stored in AmrMesh without knowing the stride at regrid time |
| All `MFCellFAB` variants | Build per-phase VoFIterators inside a phase loop; refactoring requires passing a phase-indexed vector of iterators — deferred |
| `setInvalidValue` (the `VoFIterator(IntVectSet(overlapBox), ...)` pattern) | Built from a dynamically-computed overlap box; cannot be pre-built |
