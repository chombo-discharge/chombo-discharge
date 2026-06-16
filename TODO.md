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

### Verification harness (how to actually check vectorization)

`CD_PRAGMA_SIMD` is **empty in debug builds** — vectorization only happens with
`OPT=HIGH DEBUG=FALSE` (which defines `NDEBUG`; GCC then uses `#pragma GCC ivdep`). The
active config builds with `g++ -O3 -march=native`.

To get a per-call-site verdict for one translation unit, compile it standalone with GCC's
optimization-record JSON (basic `-fopt-info` loses the inlining context because the
`BoxLoops` inner loop is `ALWAYS_INLINE`d and reported at `CD_BoxLoopsImplem.H:31`):

```bash
cd $DISCHARGE_HOME/Source
# 1. capture the real flags from the Chombo build database
make discharge-lib DIM=2 OPT=HIGH MPI=TRUE DEBUG=FALSE -pn 2>/dev/null \
  | grep -E '^CPPFLAGS :=' | head -1 | sed 's/^CPPFLAGS := //' > /tmp/cppflags.txt
make discharge-lib DIM=2 OPT=HIGH MPI=TRUE DEBUG=FALSE -pn 2>/dev/null \
  | grep -E '^CXXFLAGS :=' | head -1 | sed 's/^CXXFLAGS := //' > /tmp/cxxflags.txt
# 2. compile the single TU. NOTE: -DCH_LANG_CC is required (part of Chombo's recipe).
mpic++ $(cat /tmp/cxxflags.txt) $(cat /tmp/cppflags.txt) -DCH_LANG_CC \
  -fsave-optimization-record -c Utilities/CD_DataOps.cpp -o /tmp/probe.o
# 3. parse /tmp/probe.cpp.opt-record.json.gz: top level is [meta, passes, records];
#    each record has kind (success/failure/note), message[], location, inlining_chain[].
#    Filter location.file == CD_BoxLoopsImplem.H, attribute via innermost 'DataOps::' frame
#    in inlining_chain. message 'loop vectorized' = success.
```

Caveats: per-*function* attribution is reliable; exact per-*line* attribution is fuzzy when a
function holds several Box kernels. `not vectorized: multiple nested loops` /
`...consecutive inner loops` are about the *outer* j-loop and are expected/benign; the
meaningful blocker is `not vectorized: control flow in loop` (data-dependent branch in the
kernel) or a non-inlinable call (e.g. `std::function`).

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
- **DANGER — never pass a `std::function` as a `BoxLoops::loop` kernel.** Type-erasing the
  kernel into a `std::function<void(const IntVect&)>` (often done to select between variants in
  a `switch`) forces an **indirect, non-inlinable call on every cell**, which both kills
  performance and blocks vectorization of the entire kernel body. Watch for the idiom
  `std::function<void(...)> k; switch(...) { k = [&]{...}; }  BoxLoops::loop(box, k);`.
  Instead use a **direct lambda**. The fix when the *whole kernel* was a `std::function`: keep one
  direct lambda and, if a per-variant helper (e.g. a slope limiter) must be selected at runtime, pass
  *that helper* as a **function pointer** (captureless lambdas and `static` functions convert to
  `Real(*)(...)`). Fixed (whole-kernel std::function) instances: `CD_EBGhostCellInterpolator.cpp`,
  `CD_EBCoarseToFineInterp.cpp`.
  - **Distinct from a `std::function` used as a genuine *API callback parameter*** — e.g.
    `DataOps::setValue`'s user function — which is unavoidable and just documented as non-vectorizable.
  - **DO NOT, conversely, "fix" a plain `switch` on a loop-INVARIANT enum that calls INLINE helpers**
    (e.g. the centroid-interpolation limiter switch). Verified empirically: GCC unswitches the
    loop-invariant switch and *inlines* the tiny helpers (0 calls), whereas converting to a function
    pointer forces an indirect, non-inlined call per cell — a slight **regression**. Reverted on
    `CD_EBCentroidInterpolation.cpp` / `CD_CellCentroidInterpolationImplem.H`. Rule of thumb: a
    `std::function` *kernel* is always bad; a `switch` selecting an *inline* helper is fine.

---

## Task 1 — Vectorize `BoxLoops::loop(Box, Functor)` functors

Files sorted by occurrence count (all overloads). Triage each call for the `Box` overload.

### Source/Utilities
- [x] `Source/Utilities/CD_DataOps.cpp` (76) — DONE. Done: 3 IVS-path plain-Box loops
      converted to `loop<D_DECL(1,1,1)>` (averageCellVelocityToFaceVelocity L129,
      averageCellToFace L339, filterSmooth L748) — all verified vectorizing. Already
      vectorizing: averageCellToFace, averageCellVelocityToFaceVelocity, filterSmooth, floor,
      roof, invert, setInvalidValue. Remaining non-vectorizing kernels are pending a user
      decision on acceptable behavioral changes (see report). `squareRoot` needs
      `-fno-math-errno` (build flag, no code change; documented in Installation.rst).
      `compute`/`setValue` documented as non-vectorizable (std::function callback).
      `setCoveredValue` rewritten with a precomputed covered-cell mask → vectorizes.
      Finding: a ternary directly on `ebisbox.isRegular/isCovered(iv)` does NOT vectorize
      (out-of-line EBGraph.cpp call); a precomputed contiguous mask does. kappaSum/norm/
      getMaxMinNorm need the mask AND the postponed reduction fix (cat 4) to vectorize.
      Mask infrastructure added: PhaseRealm/Realm/AmrMesh now build+expose per-phase
      m_regularCells / m_coveredCells / m_irregularCells (EBAMRCellData, ghosted) via
      get{Regular,Covered,Irregular}Cells(realm, phase). kappaSum + norm now take the
      regular-cell mask as an argument (callers updated); setCoveredValue kept on its local
      covered-mask (already vectorizes; Realm-mask threading through its 21 callers deferred).
      Full discharge-lib build passes.
      Added a 4th mask m_notCoveredCells (1 in regular+irregular, 0 in covered) to
      PhaseRealm/Realm/AmrMesh (getNotCoveredCells). vectorLength + vectorLength2 rewritten
      branchless (lhs = notCovered * f(rhs)); all 9 call sites pass the mask. Both VECTORIZE.
      discharge-lib + ItoKMC test build pass.
      Point 4 (local accumulators): kappaSum now VECTORIZES. norm got the accumulator too but
      still does not vectorize — p>0 is blocked by std::pow (runtime int exponent), p==0
      (inf-norm max-reduction) needs -ffinite-math-only which is unsafe for a norm.
      norm, getMaxMin, getMaxMinNorm documented as non-vectorizable (min/max reductions need
      -ffinite-math-only; norm p>0 also blocked by std::pow runtime exponent).
      Point 5 (hoist component/direction loop outside the cell loop): averageFaceToCell,
      dotProduct, max all now VECTORIZE. discharge-lib builds.
      setCoveredValue: local mask removed; now takes the Realm covered-cell mask as an argument
      (m_amr->getCoveredCells(realm, phase)); all 20 cell call sites updated (the EBFluxFAB
      overload is unchanged). Uses a ternary/select (not multiply) to avoid 0*inf on covered
      cells. Vectorizes. Verified by discharge-lib + 7 physics test builds (CdrPlasma, AdvDiff,
      BrownianWalker, ItoKMC, DischargeInception, MeshODE, TracerParticles).
      -fno-math-errno added to all GNU cxxoptflags (top-level Lib/Local + betzy; fram/Intel get
      it via -Ofast). With it, squareRoot now VECTORIZES (verified). norm still does not
      (pow runtime exponent; p==0 max needs unsafe -ffinite-math-only) — could special-case
      p=1/p=2 later if wanted. CD_DataOps.cpp is COMPLETE: every kernel vectorizes or is
      documented non-vectorizable.
      VERIFY-BEFORE-MERGE: the two tagger plot-output callers (CdrPlasmaTagger, ItoKMCTagger)
      use getCoveredCells(a_outputRealm, m_phase) — confirm a_output's realm/phase matches.
- [x] `Source/Utilities/CD_VofUtils.cpp` (1) — DONE. The single loop (getAllVofsInRadius) gathers
      cut-cell VoFs via the out-of-line ebisbox.getVoFs(iv) and appends to a growing Vector;
      documented non-vectorizable (irregular/allocating/order-dependent). No raw BoxIterator loops.

### Source/AmrMesh
- [x] `Source/AmrMesh/CD_EBCoarseToFineInterp.cpp` (18) — DONE. 7 template Box loops, all
      inherently non-vectorizable and documented: interpolatePWC (coarse->fine scatter to strided
      fine cells), interpolateConservativePWC / slopeExtrapRegular / regularConstantTerm (coarse
      gather via coarsen(fineIV) — non-contiguous broadcast + floor-div branches), interpolate-
      ConservativeSlope interiorKernel (slope-limiter control flow), checkConservation x2 (debug-only,
      under #ifndef NDEBUG). Other 11 are VoF-iterator loops (not targets). The 1 raw BoxIterator is
      a nested inner refinement-box loop inside the PWC scatter kernel — left as-is (intrinsic).
      EFFICIENCY FIX: interpolateConservativeSlope's interiorKernel was a std::function selected by a
      3-case switch (per-cell indirect call); replaced with a function-pointer limiter + direct lambda
      (made superbee captureless via a nested minmod). Behavior-identical; CdrPlasma test runs clean.
      AUDIT (std::function-as-BoxLoops-kernel across all branch-modified files): only two instances
      existed — this one and CD_EBGhostCellInterpolator.cpp's regularSlopeKernel — both now fixed.
      Other std::function uses are legitimate API callbacks (DataOps compute/setValue user functions;
      ItoParticle/P& particle modifiers).
- [x] `Source/AmrMesh/CD_EBCoarAve.cpp` (18) — DONE. 6 template Box loops (arithmetic/harmonic/
      conservative averaging, cell + face), all inherently non-vectorizable: fine->coarse gather-
      reduction over the refRat^D fine cells (strided gather + inner refinement loop). Documented.
      The 3 raw BoxIterator are the cell-averaging inner refinement loops (intrinsic; face variants
      use explicit i/j/k loops) — left as-is. Multi-cut/double-counting check (per user): NOT
      applicable — irregular kernels overwrite (no double count) and use geometry-aware coarsening
      stencils that differ from the regular gather-sum, so they must run over all cut cells.
      Other 12 occurrences are VoF/Face-iterator loops (not targets).
- [x] `Source/AmrMesh/CD_EBGradient.cpp` (9) — DONE. 4 template Box loops. computeLevelGradient
      already vectorized. computeNormalDerivative NOW vectorizes — fix: hoisted the per-direction
      shift `const IntVect shift = BASISV(dir)` out of the kernel so the offset is loop-invariant
      (runtime BASISV(dir) inline was "more than one data ref"); bit-identical. defineIteratorsEBCF
      and defineStencilsEBCF are one-time regrid setup loops (data-dependent branch + IntVectSet/
      DenseIntVectSet ops) — documented non-vectorizable. Other 5 are irregular iterator loops.
      Lib builds.
- [x] `Source/AmrMesh/CD_Realm.cpp` (7) — DONE. All 7 template loops are one-time (regrid) mask/
      valid-cell setup loops, documented non-vectorizable: defineOuterHaloMask/defineOuterCFMask/
      defineValidCells (conditional bool-mask writes), defineInnerHaloMask (box-grow scatter +
      coarsen gather), defineInnerCFMask (per-cell box-grow + neighbor-overlap counting). The 1 raw
      BoxIterator (568) is the inner box-grow scatter inside flagGrownRegion (intrinsic) — left as-is.
- [~] `Source/AmrMesh/CD_PetscGrid.cpp` (7) — SKIPPED (Petsc-related, per user request).
- [x] `Source/AmrMesh/CD_EBReflux.cpp` (5) — DONE. 2 template Box loops, both non-vectorizable and
      documented: coarsenFluxesCF regular kernel (fine->coarse flux gather-reduction over fine faces,
      nested i/j/k sum) and defineRegionsCF findIrregCells (one-time setup: branch + isIrregular query
      + IntVectSet insert). The 1 raw BoxIterator (146, defineRegionsCF) is the same setup pattern
      (isRegular query + DenseIntVectSet insert) — documented, left as-is. Irregular face kernel uses
      a geometry stencil (all cut faces). Other 3 occurrences are Face/IVS-iterator loops.
- [x] `Source/AmrMesh/CD_EBLeastSquaresMultigridInterpolator.cpp` (5) — DONE. coarseFineInterpH
      vectorizes. interpOnFine NOW vectorizes — hoisted `const IntVect shift = iHiLo*BASISV(dir)` out
      of the kernel (runtime BASISV(dir) inline was "more than one data ref"; same fix as
      computeNormalDerivative); bit-identical. Earlier ivdep/dependency worry RESOLVED: interpBox is
      only the first ghost layer (1 cell thick in dir), so the stencil reads valid interior cells the
      loop never writes -> no loop-carried dependency, in-place is a safe shortcut (no transient
      buffer). applyDerivs: converted plain-Box loop(interpBox,...) (was implicit IntVectSet path) to
      template form; non-vectorizable (coarsen gather + out-of-line CoarseInterpQuadCF derivative
      calls) — documented. Lib builds. Other 2 occurrences are VoF-iterator loops.
- [x] `Source/AmrMesh/CD_EBGhostCellInterpolator.cpp` (5) — DONE. 3 template Box loops, all
      non-vectorizable and documented (coarse->fine PWL ghost interpolation, same family as
      EBCoarseToFineInterp): regSetFineToCoar + addRegularSlopeContribution are coarsen(fineIV)
      gathers (non-contiguous broadcast + floor-div); regularSlopeKernel (3 limiter variants) has
      per-cell domain-boundary checks + slope-limiter control flow. Other 2 are VoF-iterator loops.
      EFFICIENCY FIX: regularSlopeKernel was wrapped in a std::function (per-cell indirect,
      non-inlinable call) selected by a 3-case switch; replaced with a function-pointer limiter +
      a single direct lambda (behavior-identical, removes per-cell indirection, collapses the 3
      ~25-line switch arms into 1). Multi-cut check: NOT safe — irregular slopes use EB connectivity
      (getFaces/one-sided) differing from the regular Cartesian slopes even for singly-cut cells.
- [x] `Source/AmrMesh/CD_EBCentroidInterpolation.cpp` (4) — DONE. All 4 loops are VoFIterator
      (cut-cell) loops applying centroid VoFStencils — inherently irregular, not vectorization
      targets. Kernels are direct lambdas (no std::function). NO CHANGE (matches main). slopeKernel
      has a redundant inner switch on m_interpolationType (loop-invariant) to pick the limiter; I tried
      replacing it with a function-pointer limiter but reverted — empirically the switch form INLINES
      the (inline) limiters (GCC unswitches the loop-invariant switch → 0 calls) whereas a function
      pointer forces an indirect call per cell (no inline), i.e. a slight regression. So the original
      is optimal; left unchanged. No bugs (one-sided slope handling + face-count guards correct).
      KEY LESSON: a switch on a loop-INVARIANT enum calling INLINE helpers is fine (compiler
      unswitches + inlines); do NOT "fix" it with a function pointer (kills inlining). This differs
      from the std::function-as-whole-kernel case, which is always worth fixing.
- [x] `Source/AmrMesh/CD_CellCentroidInterpolationImplem.H` (4) — DONE, NO CHANGE (matches main). Same
      structure/conclusion as sibling EBCentroidInterpolation: VoFIterator cut-cell loops (not targets),
      direct-lambda kernels; the redundant inner limiter-switch is optimal as-is (inline limiters,
      compiler unswitches). No bugs.
- [x] `Source/AmrMesh/CD_EBMGRestrict.cpp` (2) — DONE. 1 template loop (restrictResidual regular
      kernel): fine->coarse gather-reduction over refRat^D fine cells (strided gather + inner
      refinement BoxIterator) — non-vectorizable, documented. The 1 raw BoxIterator is that inner
      refinement loop (intrinsic). Other loop is a VoFStencil cut-cell loop. Direct-lambda kernels (no
      std::function, no redundant dispatch); no bugs.
- [x] `Source/AmrMesh/CD_EBMGProlong.cpp` (2) — DONE. 1 template loop (prolongResidual regular kernel):
      coarse->fine scatter to refRat^D fine cells (strided writes + inner refinement BoxIterator),
      gated per fine cell by out-of-line ebisBoxFine.isIrregular(ivFine) — non-vectorizable, documented.
      The 1 raw BoxIterator is the inner refinement loop (intrinsic). Other loop is a VoFStencil cut-cell
      loop. Direct-lambda kernels; no bugs. (Note: MG-level grids, so the Realm masks don't apply, and
      the scatter wouldn't vectorize regardless.)
- [x] `Source/AmrMesh/CD_PhaseRealm.cpp` (1) — DONE. The 1 template loop is defineLevelSet's kernel:
      `fab(iv,comp) = m_baseif->value(pos)` — a virtual call on the polymorphic implicit function per
      cell (non-inlinable → non-vectorizable), one-time regrid setup. Documented. (This file also holds
      the new mask infrastructure; defineMasks uses setVal/setCoveredCellVal/IVSIterator, not a Box
      loop.) Direct-lambda kernel, no std::function, no bugs.
- [x] `Source/AmrMesh/CD_EBNonConservativeDivergence.cpp` (1) — DONE. The 1 loop is a VoFIterator
      cut-cell loop applying a VoFStencil (non-conservative divergence is inherently a cut-cell op) —
      not a vectorization target. Direct-lambda kernel, no std::function, no bugs.
- [x] `Source/AmrMesh/CD_EBFluxRedistribution.cpp` (1) — DONE. The 1 template loop is a one-time
      (regrid) setup loop with a data-dependent conditional bool-mask write (`if (mask(iv)>0)
      validCells(iv)=false`), same family as the CD_Realm valid-cell masks — documented
      non-vectorizable. Direct-lambda kernel, no std::function, no bugs.
- [x] `Source/AmrMesh/CD_AmrMesh.H` (1) — DONE. The single `BoxLoops::loop` occurrence is only a
      reference in a doc comment (line 1770), not an actual loop. Nothing to do.

### Source/Elliptic
- [x] `Source/Elliptic/CD_EBHelmholtzOp.cpp` (29)
      - BUG FIX: dotProduct fully-regular-box branch was `else if (isCovered)` (no-op) so regular
        boxes were dropped from the Krylov inner product. Now `else if (isRegular)`.
      - Vectorization wins via shift-hoist (runtime BASISV -> loop-invariant const IntVect shift):
        computeFaceCenteredFlux, refluxFreeAMROperator, computeRelaxationCoefficient (bcoef kernel).
        All verified vectorizing via opt-record.
      - Already vectorizing (compile-time strides): applyOpRegular, gauSaiRedBlackKernel (manual
        stride-2 i_start), gauSaiMultiColorKernel, computeRelaxationCoefficient inversion kernel.
      - Inherent non-vectorizers (documented in source):
        * dotProduct/norm regular kernels -- FP sum / max reduction; GCC won't reassociate without
          -fassociative-math/-ffast-math (project does not set these). Empirically verified 0 packed
          insns even after removing the isRegular guard, so the guard-removal rewrite is perf-neutral
          (NOT done, per the EBCentroid lesson).
        * applyDomainFlux Lo/Hi -- div-by-zero guard, GCC won't speculate division into masked lanes
          (thin boundary slab anyway). A mask-multiply rewrite DOES vectorize but is below the noise
          floor and weakens the defensive guard (0*inf=NaN), so not done.
        * computeFaceCentroidFlux, computeDiagWeight, applyOpIrregular, makeAggStencil, defineStencils,
          all *irregular* kernels -- sparse VoF/Face index lists + indirect stencil gather; inherently
          scalar (iterate an explicit cut-cell list).
      - PERF: precomputed domain-boundary boxes (m_domainBndryBoxes, struct DomainBndryBoxes) built
        once in defineStencils, so applyDomainFlux/fillDomainFlux no longer call EBArith::loHi every
        V-cycle. Interior patches (the common case) early-out on a single touchesDomain bool. Pure
        refactor (identical boxes); RodSphere GMG converges identically (3.28e4 -> 3.2e-7, rate ~11).
      - CONSIDERED, not done: persisting the faceFlux scratch FArrayBox in applyDomainFlux/fillDomainFlux
        to avoid per-visit malloc (Chombo BArena = plain malloc, no pool). Boundary-patches-only and the
        getFaceFlux fill (which must run regardless) dominates the alloc, so below the noise floor.
- [x] `Source/Elliptic/CD_MFHelmholtzRobinEBBC.cpp` (2)
      - Both BoxLoops are sparse VoFIterator sweeps over cut-cells (defineSinglePhase stencil build,
        applyEBFluxSinglePhase flux) with out-of-line per-cell work (LeastSquares stencils, std::function
        coefficients). Inherently scalar -- nothing to vectorize.
      - BUG FIX: in defineSinglePhase the stencil-search while-loops never set `foundStencil` from the
        returned stencil (foundStencil stayed false), so isStencilValidCF never ran and the `else` dead-cell
        branch always fired -> the Robin gradPhi stencil (the A*phi/B coupling) was cleared for EVERY cut-cell,
        silently degenerating the Robin EB BC. Added `foundStencil = (fluxStencil.size() > 0);` after each
        getInterpolationStencil call (matches the canonical Dirichlet/Neumann pattern). Same bug fixed in
        the sibling [[CD_EBHelmholtzRobinEBBC.cpp]]. Verified: builds clean; Eddington (larsen/Robin EB BC)
        runs 100 steps with no NaN/divergence.
- [x] `Source/Elliptic/CD_MFHelmholtzOp.cpp` (2)
      - Both BoxLoops are in dotProduct (the rest of the class dispatches to the per-phase EBHelmholtzOps).
        regularKernel = FP sum-reduction + isRegular guard; irregularKernel = sparse VoF reduction.
        Inherently non-vectorizable (same as EBHelmholtzOp::dotProduct) -- documented in source.
      - NO bug here (unlike EBHelmholtzOp::dotProduct): the regular kernel runs for every non-covered
        box via `if (!isCovered)`, so fully-regular boxes are correctly summed; the per-cell isRegular
        guard prevents double-counting the cut/multi-valued cells handled by the irregular kernel.
- [x] `Source/Elliptic/CD_MFHelmholtzJumpBC.cpp` (2)
      - Both BoxLoops are one-time, sparse VoFIterator stencil builders (defineStencils gradient
        stencils, buildAverageStencils AggStencil gather). Inherently scalar, not hot -- documented.
      - defineStencils correctly binds `foundStencil = getLeastSquaresBoundaryGradStencil(...)` (returns
        bool), so it does NOT have the [[CD_MFHelmholtzRobinEBBC.cpp]] bug. Verified.
      - HOT path matchBC (called every relax color via MFHelmholtzOp) is already perf-conscious: the
        heavy work is two AggStencil applies; m_denom is precomputed and folded into the stencil weights;
        m_avgVoFs precomputed. The remaining IVSIterator loop is sparse with delicate ordered
        read-after-write dependencies (the phase-1 loop reads what the phase-0 loop wrote) and size-1
        inner loops -- micro-opts there are marginal and risky, so left as-is.
- [x] `Source/Elliptic/CD_MFHelmholtzEBBC.cpp` (2)
      - Both BoxLoops are sparse VoFIterator sweeps over multi-phase cut-cells; inherently scalar.
        defineMultiPhase (setup, not hot); applyEBFluxMultiPhase (HOT -- called every applyOp).
      - getLeastSquaresBoundaryGradStencil correctly returns/binds the bool (no Robin-style bug).
      - PERF: in the hot applyEBFluxMultiPhase kernel, hoisted the loop-invariant
        m_jumpBC->getBndryPhi(m_phase, a_dit) and m_boundaryWeights[a_dit] out of the per-vof lambda.
        getBndryPhi is an out-of-line call (different TU) the compiler cannot lift itself, so it was
        re-invoked per cut-cell. Behavior-preserving (same objects bound once). Verified: clean build;
        MechShaft (multifluid electrostatics) regression converges with no NaN/divergence.
- [x] `Source/Elliptic/CD_MFHelmholtzDirichletEBBC.cpp` (2)
      - Both BoxLoops are sparse VoFIterator sweeps; inherently scalar. defineSinglePhase (setup,
        correct foundStencil binding -- no Robin bug); applyEBFluxSinglePhase (HOT, called every applyOp).
      - NO perf change (unlike the EBBC sibling): the hot kernel's only loop-invariant is
        m_boundaryWeights[a_dit], an INLINE LayoutData index the compiler already hoists -- there is no
        out-of-line per-vof container call (no getBndryPhi equivalent) to lift. Verified via asm that the
        per-vof cost is dominated by out-of-line BaseIVFAB/EBCellFAB operator() and the std::function BC.
        Hoisting would be perf-neutral, so left as-is per the no-perf-neutral rule. Documented.
- [x] `Source/Elliptic/CD_EBHelmholtzRobinEBBC.cpp` (2)
      - Same two sparse VoFIterator kernels as the MF sibling -- inherently scalar, nothing to vectorize.
      - BUG FIX: same foundStencil bug as [[CD_MFHelmholtzRobinEBBC.cpp]] (stencil always cleared);
        fixed identically. See that entry for details.
- [x] `Source/Elliptic/CD_EBHelmholtzNeumannEBBC.cpp` (2)
      - Both BoxLoops are sparse VoFIterator sweeps; inherently scalar. define (setup, clears stencils);
        applyEBFlux (HOT, called every applyOp when inhomogeneous).
      - PERF: in the hot applyEBFlux kernel, hoisted `const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit]`
        out of the per-vof lambda. getEBISL() returns an EBISLayout BY VALUE (RefCountedPtr copy w/ atomic
        refcount) and EBISLayout::operator[] is out-of-line, so this was doing an atomic layout copy + an
        out-of-line call per cut-cell. Loop-invariant; the EBISBox ref stays valid (layout owned by m_eblg,
        and the existing in-kernel code already relied on that lifetime). Behavior-preserving.
        Verified: clean build; RodSphere GMG converges identically (3.2e-7, rate ~11). No dedicated
        inhomogeneous-EB-Neumann regression exists, but the hoist is referentially transparent.
- [x] `Source/Elliptic/CD_EBHelmholtzDirichletEBBC.cpp` (2)
      - Both BoxLoops are sparse VoFIterator sweeps; inherently scalar. define (setup, correct
        foundStencil binding -- no Robin bug); applyEBFlux (HOT, called every applyOp when inhomogeneous).
      - NO perf change (same as the MF Dirichlet sibling): the hot kernel's only loop-invariant is
        m_boundaryWeights[a_dit], an inline LayoutData index the compiler hoists; areaFrac is folded into
        those precomputed weights so there is no getEBISL()/out-of-line lookup to lift (unlike the Neumann
        EBBC). Hoisting would be perf-neutral, so left as-is. Documented.
- [x] `Source/Elliptic/CD_MFHelmholtzNeumannEBBC.cpp` (1)
      - Single BoxLoop is the hot applyEBFluxSinglePhase (sparse VoFIterator, inherently scalar);
        defineSinglePhase has no stencils for Neumann.
      - PERF: same hoist as [[CD_EBHelmholtzNeumannEBBC.cpp]] -- moved the loop-invariant
        `m_eblg.getEBISL()[a_dit]` (EBISLayout-by-value copy w/ atomic refcount + out-of-line operator[])
        out of the per-vof kernel. Behavior-preserving; clean build. (MechShaft multifluid regression
        already validates the MF EB-flux machinery; this hoist is referentially transparent.)
- [x] `Source/Elliptic/CD_EBHelmholtzRobinDomainBC.cpp` (1)
      - The one BoxLoop (getFaceFlux BaseFab version) runs over a regular domain-face slab but does NOT
        vectorize: per-cell out-of-line ebisbox.isCovered(ivNear) query (gates the extrapolation) + the
        std::function coefficients in the function-BC case + division. getEBISL() is already hoisted out
        of the kernel. The remaining loop-invariant (isign*BASISV(a_dir)) hoist would not enable
        vectorization, so it's perf-neutral and skipped (thin boundary slab, like applyDomainFlux).
        Documented. (The VolIndex overload is a per-vof irregular helper -- no loop.)
- [x] `Source/Elliptic/CD_EBHelmholtzNeumannDomainBC.cpp` (1)
      - One BoxLoop (getFaceFlux BaseFab version), runs ONLY in the function-BC case: per-cell
        getBoundaryPosition() + m_functionDphiDn(pos). Non-vectorizable (out-of-line position +
        std::function indirect call), no hoist available (a_dit unused -> no getEBISL). The constant
        case is a single setVal (optimal). Documented. No code change (std::function family pattern
        left as-is per user).
- [ ] `Source/Elliptic/CD_EBHelmholtzDomainBC.cpp` (1)
- [ ] `Source/Elliptic/CD_EBHelmholtzDirichletDomainBC.cpp` (1)

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
