# TODO — Vectorization of `BoxLoops::loop(Box, Functor)`

## Instructions for Claude (read first)

This branch (`loop_factor`) is a large PR. Maintain this file as the authoritative work
list. Check items off as they are completed and leave short notes inline where a decision
was made (e.g. "cannot vectorize because …").

This session consists of:

1. **Rewrite functors for vectorization.** Rewrite the functors passed to
   `BoxLoops::loop(Box, Functor)` at *all* call sites so that the innermost loop
   auto-vectorizes, **or** investigate and document *why* a given functor cannot vectorize
   (e.g. data dependencies, indirect addressing, EB/irregular access, function calls that
   are not inlinable). Record the reason inline in this file next to the item.

2. **Normalize raw `BoxIterator` loops.** For every hand-written `for (BoxIterator …)` loop,
   *either* inline it into a direct triple/`D_TERM` index loop, *or* rewrite it to the
   standard `BoxLoops::loop<D_DECL(1,1,1)>(box, kernel)` constructor — whichever keeps the
   call site readable and consistent with the surrounding code.

3. **Documentation audit before merge.** Before the PR is merged, run a documentation audit
   over all `.rst` and Doxygen files. Specifically, run `CheckDocs.py` to list files changed
   on this branch that are referenced by `.. literalinclude::` directives, and verify that
   every such literalinclude still points to the *same* code block as it did on `main`
   (line ranges / `:start-after:` / `:end-before:` anchors must still bracket the intended
   code). Fix any drift.

### Reference

- Loop definitions: `Source/AmrMesh/CD_BoxLoops.H`, `Source/AmrMesh/CD_BoxLoopsImplem.H`
- The `Box` overload is the vectorization target. Strides are compile-time non-type template
  parameters; unit-stride is `BoxLoops::loop<D_DECL(1,1,1)>(box, kernel)`.
- Other overloads (`IntVectSet`, `DenseIntVectSet`, `VoFIterator`, `FaceIterator`,
  `Vector<T>`) are inherently irregular and are **not** vectorization targets — only audit
  them when normalizing raw iterator loops.
- Audit tool: `python3 CheckDocs.py` (compares against `main`).

### Working notes / conventions

- **Work one file at a time.** Do not start a new file until the current one is finished and
  checked off.
- **Only tick a file off once its loops have been *verified* to vectorize** (e.g. via
  compiler vectorization reports / inspection of the generated code), not merely rewritten.
  If a loop cannot vectorize, the file is "done" only once the reason is documented (see
  below) and the user has been consulted where required.
- **The user must be notified whenever the behavior of a loop is expected to change.**
  Changing the *contents* of a loop is normally **not acceptable** — vectorization work
  should preserve behavior. Consult the user when or if a loop is rewritten in a way that
  could alter what it computes.
- **If a loop CAN'T be vectorized**, it may be acceptable to write the loop directly,
  inline, *without* the `BoxLoops::loop` wrapper. In that case there must be **extensive
  documentation above the call site** explaining why the `BoxLoops::loop` wrapper was not
  used and why the loop cannot vectorize.
- Per item, the count in parentheses is the total `BoxLoops::loop` occurrences in that file
  (all overloads), not only the `Box` overload — triage each call to see if it is a `Box`
  loop before acting.
- Prefer rewrites that keep dimension independence via `D_DECL` / `D_TERM`.
- Do not change numerical behavior. Vectorization rewrites must be bit-for-bit equivalent
  (modulo FP reassociation that the existing `CD_PRAGMA_SIMD` already permits).

---

## Task 1 — Vectorize `BoxLoops::loop(Box, Functor)` functors

Files sorted by occurrence count (all overloads). Triage each call for the `Box` overload.

### Source/Utilities
- [ ] `Source/Utilities/CD_DataOps.cpp` (72)
- [ ] `Source/Utilities/CD_VofUtils.cpp` (1)

### Source/AmrMesh
- [ ] `Source/AmrMesh/CD_EBCoarseToFineInterp.cpp` (18)
- [ ] `Source/AmrMesh/CD_EBCoarAve.cpp` (18)
- [ ] `Source/AmrMesh/CD_EBGradient.cpp` (9)
- [ ] `Source/AmrMesh/CD_Realm.cpp` (7)
- [ ] `Source/AmrMesh/CD_PetscGrid.cpp` (7)
- [ ] `Source/AmrMesh/CD_EBReflux.cpp` (5)
- [ ] `Source/AmrMesh/CD_EBLeastSquaresMultigridInterpolator.cpp` (5)
- [ ] `Source/AmrMesh/CD_EBGhostCellInterpolator.cpp` (5)
- [ ] `Source/AmrMesh/CD_EBCentroidInterpolation.cpp` (4)
- [ ] `Source/AmrMesh/CD_CellCentroidInterpolationImplem.H` (4)
- [ ] `Source/AmrMesh/CD_EBMGRestrict.cpp` (2)
- [ ] `Source/AmrMesh/CD_EBMGProlong.cpp` (2)
- [ ] `Source/AmrMesh/CD_PhaseRealm.cpp` (1)
- [ ] `Source/AmrMesh/CD_EBNonConservativeDivergence.cpp` (1)
- [ ] `Source/AmrMesh/CD_EBFluxRedistribution.cpp` (1)

### Source/Elliptic
- [ ] `Source/Elliptic/CD_EBHelmholtzOp.cpp` (29)
- [ ] `Source/Elliptic/CD_MFHelmholtzRobinEBBC.cpp` (2)
- [ ] `Source/Elliptic/CD_MFHelmholtzOp.cpp` (2)
- [ ] `Source/Elliptic/CD_MFHelmholtzJumpBC.cpp` (2)
- [ ] `Source/Elliptic/CD_MFHelmholtzEBBC.cpp` (2)
- [ ] `Source/Elliptic/CD_MFHelmholtzDirichletEBBC.cpp` (2)
- [ ] `Source/Elliptic/CD_EBHelmholtzRobinEBBC.cpp` (2)
- [ ] `Source/Elliptic/CD_EBHelmholtzNeumannEBBC.cpp` (2)
- [ ] `Source/Elliptic/CD_EBHelmholtzDirichletEBBC.cpp` (2)
- [ ] `Source/Elliptic/CD_MFHelmholtzNeumannEBBC.cpp` (1)
- [ ] `Source/Elliptic/CD_EBHelmholtzRobinDomainBC.cpp` (1)
- [ ] `Source/Elliptic/CD_EBHelmholtzNeumannDomainBC.cpp` (1)
- [ ] `Source/Elliptic/CD_EBHelmholtzDomainBC.cpp` (1)
- [ ] `Source/Elliptic/CD_EBHelmholtzDirichletDomainBC.cpp` (1)
- [ ] `Source/Elliptic/CD_EddingtonSP1.cpp` — **NOTE: deleted on this branch (`git status` shows `D`)**; confirm before touching.

### Source/ConvectionDiffusionReaction
- [ ] `Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp` (30)
- [ ] `Source/ConvectionDiffusionReaction/CD_CdrCTU.cpp` (10)
- [ ] `Source/ConvectionDiffusionReaction/CD_CdrGodunov.cpp` (2)

### Source/Electrostatics
- [ ] `Source/Electrostatics/CD_FieldSolver.cpp` (12)
- [ ] `Source/Electrostatics/CD_MFHelmholtzElectrostaticEBBC.cpp` (2)

### Source/Particle
- [ ] `Source/Particle/CD_EBCoarseFineParticleMesh.cpp` (12)
- [ ] `Source/Particle/CD_EBParticleMeshImplem.H` (4)
- [ ] `Source/Particle/CD_ParticleContainerImplem.H` (3)
- [ ] `Source/Particle/CD_EBAMRParticleMesh.cpp` (2)

### Source — other modules
- [ ] `Source/ItoDiffusion/CD_ItoSolver.cpp` (8)
- [ ] `Source/RadiativeTransfer/CD_McPhoto.cpp` (7)
- [ ] `Source/RadiativeTransfer/CD_EddingtonSP1.cpp` (5)
- [ ] `Source/Driver/CD_Driver.cpp` (6)
- [ ] `Source/MeshODESolver/CD_MeshODESolverImplem.H` (4)
- [ ] `Source/SurfaceODESolver/CD_SurfaceODESolverImplem.H` (3)

### Physics
- [ ] `Physics/DischargeInception/CD_DischargeInceptionStepperImplem.H` (34)
- [ ] `Physics/ItoKMC/CD_ItoKMCStepperImplem.H` (22)
- [ ] `Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp` (19)
- [ ] `Physics/ItoKMC/TimeSteppers/ItoKMCBackgroundEvaluator/CD_ItoKMCBackgroundEvaluatorImplem.H` (5)
- [ ] `Physics/DischargeInception/CD_DischargeInceptionTagger.cpp` (4)
- [ ] `Physics/CdrPlasma/CD_CdrPlasmaTagger.cpp` (4)
- [ ] `Physics/ItoKMC/CD_ItoKMCTaggerImplem.H` (2)
- [ ] `Physics/ItoKMC/CD_ItoKMCFieldTaggerImplem.H` (2)
- [ ] `Physics/CdrPlasma/CD_CdrPlasmaFieldTagger.cpp` (2)
- [ ] `Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp` (1)
- [ ] `Physics/AdvectionDiffusion/CD_AdvectionDiffusionTagger.cpp` (1)

---

## Task 2 — Normalize raw `BoxIterator` loops

- [ ] `Source/Utilities/CD_DataOps.cpp` (3)
- [ ] `Source/Utilities/CD_DischargeIO.cpp` (2)
- [ ] `Source/AmrMesh/CD_EBCoarAve.cpp` (3)
- [ ] `Source/AmrMesh/CD_TiledMeshRefine.cpp` (2)
- [ ] `Source/AmrMesh/CD_CoarseInterpQuadCF.cpp` (2)
- [ ] `Source/AmrMesh/CD_Realm.cpp` (1)
- [ ] `Source/AmrMesh/CD_EBReflux.cpp` (1)
- [ ] `Source/AmrMesh/CD_EBMGRestrict.cpp` (1)
- [ ] `Source/AmrMesh/CD_EBMGProlong.cpp` (1)
- [ ] `Source/AmrMesh/CD_EBCoarseToFineInterp.cpp` (1)
- [ ] `Source/Particle/CD_EBCoarseFineParticleMesh.cpp` (2)
- [ ] `Source/Geometry/CD_ScanShopImplem.H` (2)

---

## Task 3 — Documentation audit (pre-merge)

- [ ] Run `python3 CheckDocs.py` and collect all `.rst` literalincludes referencing files
      changed on this branch.
- [ ] For each hit, diff the referenced code block against `main` and confirm the
      literalinclude line ranges / anchors still bracket the intended code. Fix drift.
- [ ] Verify Doxygen builds clean: `pre-commit run doxygen-check --all-files`.
