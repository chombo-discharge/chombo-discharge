# DataOps Reorganization — Checkpoint List

## Status — update at the start/end of each session

**Branch:** `dataops` | **Last updated:** 2026-06-12 | **No code written yet.**

### Planning — complete (do not re-litigate these decisions)

- [x] Checkpoint 1 — Literalinclude inventory compiled (reference list; act on it in Checkpoint 6.1)
- [x] Checkpoint 2 — AmrMesh/PhaseRealm extension design finalized
- [x] Checkpoint 3 — DataOps signature change design finalized
- [x] Checkpoint 4 — BoxLoops internal reorganization design finalized
- [x] Checkpoint 5 — Call site inventory compiled

### Implementation — not started

- [ ] **2-code** — Write `PhaseRealm` / `Realm` / `AmrMesh` changes (6 files; see §2 task list)
- [ ] **3+4-code** — Write `DataOps` signature + implementation changes (see §3+4 task list)
- [ ] **5-sites** — Update all call sites in Checkpoint 5 (check off each site as updated)
- [ ] **6-cleanup** — Pre-merge cleanup (do last)

**Resume here:** `2-code` — start with `Source/AmrMesh/CD_PhaseRealm.H`.

> **Keep this current.** When finishing a session, update "Resume here" and the last-updated
> date so the next session (on any machine) can orient instantly.

---

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

### 2-code task list

Check off each file once the code change is written, reviewed, and compiles cleanly.

- [ ] `Source/AmrMesh/CD_PhaseRealm.H` — add `m_multiCutVofIter` member + `getMultiCutVofIterator()` declaration
- [ ] `Source/AmrMesh/CD_PhaseRealm.cpp` — populate in `regridOperators` (`getMultiCells`); clear in `preRegrid`
- [ ] `Source/AmrMesh/CD_Realm.H` — add `getMultiCutVofIterator(phase)` declaration
- [ ] `Source/AmrMesh/CD_Realm.cpp` — implement one-liner delegation to `PhaseRealm`
- [ ] `Source/AmrMesh/CD_AmrMesh.H` — add `getMultiCutVofIterator(realm, phase)` declaration
- [ ] `Source/AmrMesh/CD_AmrMesh.cpp` — implement delegation to `Realm`

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

### 3+4-code task list

One checkbox per function that changes signature and/or internal loop. Check off once both
`CD_DataOps.H` (declaration) and `CD_DataOps.cpp` (implementation) are updated and compile.
Functions marked **no change** and **keep local** are omitted — they need no edits.

**Full cut-cell VoFIterator** — add `const Vector<RefCountedPtr<LayoutData<VoFIterator>>>&` parameter:

- [ ] `averageCellToFace` — H + cpp
- [ ] `averageCellVelocityToFaceVelocity` — H + cpp
- [ ] `averageFaceToCell` — H + cpp
- [ ] `floor(EBAMRCellData / LevelData<EBCellFAB>, Real)` — H + cpp (full-cut; ghost-region caveat)
- [ ] `roof(EBAMRCellData / LevelData<EBCellFAB>, Real)` — H + cpp (full-cut; ghost-region caveat)
- [ ] `kappaScale(LevelData<EBCellFAB>)` — H + cpp
- [ ] `kappaSum` — H + cpp
- [ ] `norm` — H + cpp
- [ ] `volumeScale` — H + cpp

**Multi-cut VoFIterator** — add `const Vector<RefCountedPtr<LayoutData<VoFIterator>>>&` parameter
(pass `AmrMesh::getMultiCutVofIterator` at call sites):

- [ ] `compute` — H + cpp
- [ ] `divideFallback` — H + cpp
- [ ] `dotProduct` — H + cpp
- [ ] `getMaxMin(LevelData<EBCellFAB>, int)` — H + cpp
- [ ] `getMaxMinNorm(LevelData<EBCellFAB>)` — H + cpp
- [ ] `max` — H + cpp
- [ ] `setInvalidValue(EBAMRCellData, Vector<int>, Real)` — H + cpp
- [ ] `setValue(LevelData<EBCellFAB>, function<Real(RealVect)>, ...)` — H + cpp
- [ ] `setValue(LevelData<EBCellFAB>, function<RealVect(RealVect)>, ...)` — H + cpp
- [ ] `vectorLength(EBCellFAB, EBCellFAB, Box)` — H + cpp
- [ ] `vectorLength2(EBCellFAB, EBCellFAB, Box)` — H + cpp

---

## Checkpoint 5 — Call sites of DataOps functions

All external call sites of `DataOps::` across `Source/`, `Physics/`, `Exec/`, and
`Geometries/` (excluding `CD_DataOps.H`, `CD_DataOps.cpp`, `CD_DataOpsImplem.H`).

Call sites are grouped by function. Each function header notes the change category from
Checkpoints 3 and 4:

- **No change** — signature unchanged; call site requires no update.
- **Multi-cut VoFIter** — caller must pass `AmrMesh::getMultiCutVofIterator(realm, phase)`.
- **Full-cut VoFIter** — caller must pass `AmrMesh::getVofIterator(realm, phase)`.
- **Keep local** — `DataOps` retains on-the-fly VoFIterator construction; no call-site change.

Check off each site after verifying or updating the call.

---

### `DataOps::averageCellToFace` — Full-cut VoFIter (EB-specific irregular kernel)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:251
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:1666
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2961
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2972
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5265
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3453
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1202
- [ ] Source/AmrMesh/CD_AmrMesh.cpp:1404
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:540
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:558

---

### `DataOps::averageCellVelocityToFaceVelocity` — Full-cut VoFIter (EB-specific irregular kernel)

- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:384

---

### `DataOps::averageFaceToCell` — Full-cut VoFIter (calls `ebisbox.getFaces`)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2854
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1849
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1855

---

### `DataOps::axby` — No change (delegates to `EBCellFAB::axby`)

- [ ] Source/Elliptic/CD_EBHelmholtzOp.cpp:649

---

### `DataOps::compute` — Multi-cut VoFIter (geometry-independent transform)

- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3143
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3144
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3158
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3159
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3210
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3215

---

### `DataOps::computeMinValidBox` — No change (pure box/particle geometry, no VoFIterator)

- [ ] Physics/ItoKMC/CD_ItoKMCPhysicsImplem.H:701
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:4483
- [ ] Physics/ItoKMC/PlasmaModels/ItoKMCJSON/CD_ItoKMCJSON.cpp:3553
- [ ] Source/AmrMesh/CD_CellInfo.cpp:42
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:665
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1331
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:1104
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:1217

---

### `DataOps::copy` — No change (delegates to `localCopyTo`)

- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:383
- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:424
- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionTagger.cpp:76
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:484
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2644
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2836
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2845
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4596
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:405
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:406
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:599
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:850
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1272
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:535
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:549
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:670
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:685
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:839
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:842
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:857
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:861
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:888
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:934
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:941
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:950
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:954
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:984
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:991
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1070
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1084
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1091
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:2028
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:2038
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:2058
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:2068
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:265
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:275
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:276
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:289
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:299
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:300
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:322
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1724
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1751
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3140
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3141
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3209
- [ ] Physics/DischargeInception/CD_DischargeInceptionTagger.cpp:182
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1996
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2066
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3892
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3902
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5262
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5272
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:242
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:243
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:252
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:336
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:579
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:580
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1155
- [ ] Source/ConvectionDiffusionReaction/CD_CdrCTU.cpp:200
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:203
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:232
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:238
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:289
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:302
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:331
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:337
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:235
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:260
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2475
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2727
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:197
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:289
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:304
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:361
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:609
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:697
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:1148
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:402
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:464
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:597
- [ ] Source/TracerParticles/CD_TracerParticleSolverImplem.H:302

---

### `DataOps::divideByScalar` — No change (delegates to `divide`, which uses `EBCellFAB::divide`)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:205
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4256
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:1153

---

### `DataOps::divideFallback` — Multi-cut VoFIter (geometry-independent; clone pattern)

- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1325
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1372
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1839
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2090
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:254
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2469
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1557
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1758
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1789
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1821

---

### `DataOps::dotProduct` — Multi-cut VoFIter (geometry-independent sum of products)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4211
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1156
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1165
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1339
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:301

---

### `DataOps::filterSmooth` — Keep local (stride-dependent grown IVS; cannot pre-build)

- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:610
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:619
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1173
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1312

---

### `DataOps::floor` (EBCellFAB variant) — Full-cut VoFIter (operates over full ghost region)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:266
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2588
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:900
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1332
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1419
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1428
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1675
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1685
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:936
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1002
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1544
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1581
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1699
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1763
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5210
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1760

---

### `DataOps::getMaxMin` — Multi-cut VoFIter (no `isRegular` guard; BoxLoops covers VoF 0)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:3226
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:337
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1657
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1658
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2373
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2374
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1326
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1373
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1844
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5229
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:261
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:271
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:272
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2960

---

### `DataOps::getMaxMinNorm` — Multi-cut VoFIter (covered-cell mask; regular loop handles cut-cells)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaFieldTagger.cpp:114
- [ ] Physics/CdrPlasma/CD_CdrPlasmaFieldTagger.cpp:115
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:3013
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4262
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2488
- [ ] Physics/ItoKMC/CD_ItoKMCFieldTaggerImplem.H:114
- [ ] Physics/ItoKMC/CD_ItoKMCFieldTaggerImplem.H:115
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2097

---

### `DataOps::incr` — No change (delegates to `EBCellFAB::plus` / self-contained `BaseIVFAB`)

- [ ] Exec/Convergence/AdvectionDiffusion/C1/main.cpp:100
- [ ] Exec/Convergence/AdvectionDiffusion/C2/main.cpp:72
- [ ] Exec/Convergence/CdrPlasma/C1/main.cpp:140
- [ ] Exec/Convergence/Electrostatics/C1/main.cpp:103
- [ ] Exec/Convergence/Electrostatics/C2/main.cpp:108
- [ ] Exec/Convergence/Electrostatics/C3/main.cpp:107
- [ ] Exec/Convergence/RadiativeTransfer/C1/main.cpp:100
- [ ] Exec/Convergence/RadiativeTransfer/C2/main.cpp:74
- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:384
- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:388
- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:389
- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:432
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:149
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:209
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:263
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:335
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:336
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2592
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2838
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2860
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:3955
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4043
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4130
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:852
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1139
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1226
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1235
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1241
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1242
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1250
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1286
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1291
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1295
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1296
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1317
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1357
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1368
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1412
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:496
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:517
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:840
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:850
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:862
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:873
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:877
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:889
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:894
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:899
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:935
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:955
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:987
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1007
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1008
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1062
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1072
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1113
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1120
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1140
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1169
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1820
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:314
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:315
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:330
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1140
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1141
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1146
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1147
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1154
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1155
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1163
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1164
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1588
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1589
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1595
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1596
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1708
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1725
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1729
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1730
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5287
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5288
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1919
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1930
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1980
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1998
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3968
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5157
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5166
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5189
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5204
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5665
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:253
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1138
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1158
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1211
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1254
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1259
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1295
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1351
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1352
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1666
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1668
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1752
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1757
- [ ] Physics/MeshODE/CD_MeshODEStepperImplem.H:425
- [ ] Physics/MeshODE/CD_MeshODEStepperImplem.H:446
- [ ] Physics/MeshODE/CD_MeshODEStepperImplem.H:449
- [ ] Physics/MeshODE/CD_MeshODEStepperImplem.H:450
- [ ] Source/ConvectionDiffusionReaction/CD_CdrGodunov.cpp:249
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:205
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:293
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:307
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1228
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1884
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:1674
- [ ] Source/Elliptic/CD_EBHelmholtzOp.cpp:637
- [ ] Source/Elliptic/CD_MFHelmholtzOp.cpp:355
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:715
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:111
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:731
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:761
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:841
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:856

---

### `DataOps::kappaScale` — Full-cut VoFIter (only irregular kernel; uses `volFrac`)

- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:420
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1747
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:204
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:290
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:295
- [ ] Source/Elliptic/CD_EBHelmholtzOp.cpp:1503
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:616
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:717

---

### `DataOps::kappaSum` — Full-cut VoFIter (regular kernel guards `isRegular`; volFrac in irregular)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4219

---

### `DataOps::multiply` — No change (delegates to `EBCellFAB::operator*=`)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:206
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3146
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3147
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3211
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1156
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2952
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2970

---

### `DataOps::multiplyScalar` — No change (delegates to `EBCellFAB::mult`)

- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2593
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2645
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2837
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2857
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5189
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1997
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2067
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2879
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2728

---

### `DataOps::norm` — Full-cut VoFIter (regular kernel guards `isRegular`)

- [ ] Exec/Convergence/AdvectionDiffusion/C1/main.cpp:103
- [ ] Exec/Convergence/AdvectionDiffusion/C1/main.cpp:104
- [ ] Exec/Convergence/AdvectionDiffusion/C1/main.cpp:105
- [ ] Exec/Convergence/AdvectionDiffusion/C2/main.cpp:76
- [ ] Exec/Convergence/AdvectionDiffusion/C2/main.cpp:77
- [ ] Exec/Convergence/AdvectionDiffusion/C2/main.cpp:78
- [ ] Exec/Convergence/CdrPlasma/C1/main.cpp:143
- [ ] Exec/Convergence/CdrPlasma/C1/main.cpp:144
- [ ] Exec/Convergence/CdrPlasma/C1/main.cpp:145
- [ ] Exec/Convergence/Electrostatics/C1/main.cpp:106
- [ ] Exec/Convergence/Electrostatics/C1/main.cpp:107
- [ ] Exec/Convergence/Electrostatics/C1/main.cpp:108
- [ ] Exec/Convergence/Electrostatics/C2/main.cpp:111
- [ ] Exec/Convergence/Electrostatics/C2/main.cpp:112
- [ ] Exec/Convergence/Electrostatics/C2/main.cpp:113
- [ ] Exec/Convergence/Electrostatics/C3/main.cpp:110
- [ ] Exec/Convergence/Electrostatics/C3/main.cpp:111
- [ ] Exec/Convergence/Electrostatics/C3/main.cpp:112
- [ ] Exec/Convergence/RadiativeTransfer/C1/main.cpp:103
- [ ] Exec/Convergence/RadiativeTransfer/C1/main.cpp:104
- [ ] Exec/Convergence/RadiativeTransfer/C1/main.cpp:105
- [ ] Exec/Convergence/RadiativeTransfer/C2/main.cpp:77
- [ ] Exec/Convergence/RadiativeTransfer/C2/main.cpp:78
- [ ] Exec/Convergence/RadiativeTransfer/C2/main.cpp:79
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1145
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1146

---

### `DataOps::plus` — No change (delegates to `EBCellFAB::plus`)

- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5985

---

### `DataOps::scale` — No change (delegates to `EBCellFAB::mult`)

- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp:421
- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionTagger.cpp:92
- [ ] Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp:392
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:154
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:215
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2865
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1331
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:863
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1009
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1057
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1748
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2640
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2824
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5145
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1934
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2002
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2818
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2849
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5984
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1163
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1673
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:562
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1843
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1872
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2953
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:198
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:290
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:305
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:362
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:885
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:913
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:944
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:979
- [ ] Source/Elliptic/CD_EBHelmholtzOp.cpp:657
- [ ] Source/Elliptic/CD_MFHelmholtzOp.cpp:363
- [ ] Source/Elliptic/CD_MFHelmholtzOpFactory.cpp:133
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:610
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:1062
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:1154

---

### `DataOps::setCoveredValue` — No change (loops via `isCovered` test, no VoFIterator)

- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionTagger.cpp:82
- [ ] Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp:177
- [ ] Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp:641
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:161
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4505
- [ ] Physics/CdrPlasma/CD_CdrPlasmaTagger.cpp:179
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5421
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1093
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:4788
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:4789
- [ ] Physics/ItoKMC/CD_ItoKMCTaggerImplem.H:179
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1229
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1938
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:3102
- [ ] Source/Driver/CD_Driver.cpp:2678
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1611
- [ ] Source/MeshODESolver/CD_MeshODESolverImplem.H:692
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:653
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:664
- [ ] Source/RadiativeTransfer/CD_RtSolver.cpp:370
- [ ] Source/TracerParticles/CD_TracerParticleSolverImplem.H:582

---

### `DataOps::setInvalidValue` — Keep local (overlap box computed dynamically at call time)

- [ ] Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp:640

---

### `DataOps::setValue` — Mixed: `(LD, Real)` overloads need no change; `(LD, function, ...)` overloads need multi-cut VoFIter

- [ ] Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp:170
- [ ] Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp:171
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:135
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:188
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:242
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:243
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:262
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:395
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:402
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:1655
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2591
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2602
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2654
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2766
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:2808
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:3943
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4031
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4118
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:4255
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:598
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:652
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:653
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:654
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:655
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1127
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1273
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.cpp:1347
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:494
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:515
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:712
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:982
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1006
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1016
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1060
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1112
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1119
- [ ] Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp:1813
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:194
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:195
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:206
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:207
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:208
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:209
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:210
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:211
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:222
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:223
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:224
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:225
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:226
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:313
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1139
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1145
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1153
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1162
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1172
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1177
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1349
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1354
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1587
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1594
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2204
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2513
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:2514
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3069
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3070
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3072
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3073
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3202
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3896
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:3897
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:4119
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:4120
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5262
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5263
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:5286
- [ ] Physics/DischargeInception/CD_DischargeInceptionTagger.cpp:80
- [ ] Physics/Electrostatics/CD_FieldStepperImplem.H:229
- [ ] Physics/Electrostatics/CD_FieldStepperImplem.H:230
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:649
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1125
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1774
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1905
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1965
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2089
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2952
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:2953
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3442
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3497
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3617
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3722
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3905
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3906
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:4826
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:4893
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5199
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5275
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5276
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5622
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5623
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:5829
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:200
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:201
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1120
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1132
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1143
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1193
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1194
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1242
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1280
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1281
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1648
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H:1746
- [ ] Physics/TracerParticle/CD_TracerParticleStepperImplem.H:370
- [ ] Physics/TracerParticle/CD_TracerParticleStepperImplem.H:421
- [ ] Source/ConvectionDiffusionReaction/CD_CdrCTU.cpp:195
- [ ] Source/ConvectionDiffusionReaction/CD_CdrGodunov.cpp:236
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:181
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:261
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:355
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:356
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:405
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:515
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:516
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:550
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:588
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:634
- [ ] Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp:675
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:230
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:255
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:280
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:281
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:290
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:291
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:292
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:308
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:309
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:310
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:325
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:326
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:327
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:328
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:329
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:378
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:427
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:688
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1196
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1227
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1260
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1468
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1577
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1578
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1579
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1594
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1595
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1596
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1620
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1663
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1677
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1722
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1738
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1782
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:1883
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2332
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2468
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2716
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2974
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:3127
- [ ] Source/Driver/CD_Driver.cpp:2124
- [ ] Source/Driver/CD_Driver.cpp:2304
- [ ] Source/Driver/CD_Driver.cpp:2404
- [ ] Source/Driver/CD_Driver.cpp:2453
- [ ] Source/Driver/CD_Driver.cpp:2656
- [ ] Source/Driver/CD_Driver.cpp:2883
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:152
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:153
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:154
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:155
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:156
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:385
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:386
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:387
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:388
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:470
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:484
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:498
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:511
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:1063
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:1064
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:1065
- [ ] Source/Electrostatics/CD_FieldSolver.cpp:1673
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:286
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:815
- [ ] Source/Electrostatics/CD_FieldSolverGMG.cpp:816
- [ ] Source/Elliptic/CD_EBHelmholtzOp.cpp:740
- [ ] Source/Elliptic/CD_MFHelmholtzOp.cpp:371
- [ ] Source/Elliptic/CD_MFHelmholtzOpFactory.cpp:173
- [ ] Source/Elliptic/CD_MFHelmholtzOpFactory.cpp:176
- [ ] Source/Elliptic/CD_MFHelmholtzOpFactory.cpp:346
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1123
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1664
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:1904
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:2104
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:2116
- [ ] Source/MeshODESolver/CD_MeshODESolverImplem.H:184
- [ ] Source/MeshODESolver/CD_MeshODESolverImplem.H:269
- [ ] Source/MeshODESolver/CD_MeshODESolverImplem.H:540
- [ ] Source/Particle/CD_EBAMRParticleMeshImplem.H:87
- [ ] Source/Particle/CD_EBAMRParticleMeshImplem.H:151
- [ ] Source/Particle/CD_EBAMRParticleMeshImplem.H:252
- [ ] Source/Particle/CD_EBAMRParticleMeshImplem.H:338
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:503
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:504
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:505
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:696
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:1034
- [ ] Source/RadiativeTransfer/CD_EddingtonSP1.cpp:1035
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:65
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:717
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:728
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:908
- [ ] Source/RadiativeTransfer/CD_McPhoto.cpp:1310
- [ ] Source/RadiativeTransfer/CD_RtSolver.cpp:223
- [ ] Source/RadiativeTransfer/CD_RtSolver.cpp:234
- [ ] Source/RadiativeTransfer/CD_RtSolver.cpp:248
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:376
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:389
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:438
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:451
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:730
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:760
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:840
- [ ] Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H:855
- [ ] Source/TracerParticles/CD_TracerParticleSolverImplem.H:370
- [ ] Source/TracerParticles/CD_TracerParticleSolverImplem.H:424

---

### `DataOps::sgn` — No change (simple sign function, no VoFIterator)

- [ ] (no external call sites found)

---

### `DataOps::squareRoot` — Keep local (MFCellFAB variant; per-phase VoFIterator loop deferred)

- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1157
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1166
- [ ] Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H:1340
- [ ] Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp:2954

---

### `DataOps::vectorLength` — Multi-cut VoFIter (geometry-independent length formula)

- [ ] Physics/AdvectionDiffusion/CD_AdvectionDiffusionTagger.cpp:81
- [ ] Physics/CdrPlasma/CD_CdrPlasmaFieldTagger.cpp:73
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:185
- [ ] Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp:204
- [ ] Physics/ItoKMC/CD_ItoKMCFieldTaggerImplem.H:83
- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:1836
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:249
- [ ] Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H:250
- [ ] Source/ItoDiffusion/CD_ItoSolver.cpp:2220

---

### `DataOps::volumeScale` — Full-cut VoFIter (uses `volFrac`; geometry-specific irregular kernel)

- [ ] Physics/ItoKMC/CD_ItoKMCStepperImplem.H:3964
- [ ] Source/TracerParticles/CD_TracerParticleSolverImplem.H:375

---

## Checkpoint 6 — Pre-merge cleanup

Work through these steps in order before opening the PR. Each item is a gate: do not proceed
to the next until the current one is green.

### 6.1 — Fix broken literalincludes

For every file listed in Checkpoint 1, open the corresponding RST file and verify that the
`:lines:` / `:emphasize-lines:` range still selects the intended code block. Update any range
that has shifted due to the reorganization.

- [ ] Run `pre-commit run check-literalincludes --all-files` and confirm zero failures.
- [ ] Run `doxygen Docs/doxygen.conf` and confirm zero warnings (or run
  `pre-commit run doxygen-check --all-files`).

### 6.2 — Generate benchmark references on `main` and run the test suite

Benchmark reference files must be produced from an unmodified `main` build so that any
regression is attributable solely to this PR.

- [ ] Check out `main`, build, and run each regression test listed in `Exec/Tests/` to
  produce reference output (HDF5 plot files or convergence tables).
- [ ] Check out this branch, build with identical flags, and run the same tests.
- [ ] Diff the outputs. Confirm that no benchmark values change beyond floating-point
  round-off (identical grid, same flags ⇒ bit-identical output is the target).
- [ ] Confirm that the CI compilation jobs (`Linux-GNU`, `Linux-oneAPI`) pass.

### 6.3 — Pre-commit hooks

- [ ] Run `pre-commit run --all-files` and fix every failure before continuing.
  Key hooks to watch: `clang-format`, `reuse`, `codespell`, `format-input-files`.
- [ ] If `clang-tidy` is enabled locally, run it and resolve any new warnings introduced by
  this PR (existing suppressions are acceptable; new ones must be justified).

### 6.4 — New-file audit

At this point, enumerate **every file that exists on this branch but not on `main`** (use
`git diff --name-only --diff-filter=A main`) and ask the user, file by file, whether each
new file should be:

- **Kept** — it is a genuine deliverable and should be merged.
- **Deleted** — it is scaffolding, scratch work, or a planning artefact (e.g. `TODO.md`)
  that must not appear in the merged history.

Delete the unwanted files, stage the deletion, and amend or add a new commit before the
final push. Do **not** merge until every new file has been explicitly classified.

- [ ] Remove the `## 6. Current work in progress` section from `CLAUDE.md` (added on this
  branch to help Claude pick up context across sessions; must not appear on `main`).
