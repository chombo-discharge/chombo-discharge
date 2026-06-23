# SoA particle migration: Dev/ -> Source/Particle/

Tracking the move of the reviewed SoA staging code into production `Source/Particle/`.
One file at a time, on review. Each migrated file: switch `"CD_X.H"` includes to the
`<CD_X.H>` Source convention, drop it from the `Docs/doxygen.conf` EXCLUDE (so it is
doxygen-checked), and keep REUSE/clang-format/doxygen green.

## Library headers (dependency order)

- [x] **CD_ParticleLoops.H** -> `Source/Particle/` (reviewed). Note: references `CD_ParticleSoA.H`,
      which is not migrated yet; nothing in `Source/` includes ParticleLoops, so it does not yet
      compile from the library -- it is exercised by the Dev demo/benchmark (which carry `-IDev` for
      the not-yet-migrated headers). Compile-from-Source coverage arrives once its dependencies move.
- [x] **CD_ParticleSoA.H** + **CD_ParticleSoAImplem.H** -> `Source/Particle/` (reviewed; doxygen-checked).
      This resolves CD_ParticleLoops.H's forward dependency (both now in Source). The container's
      `ParticleContainerSoA` bulk paths (distributeFromPool, transferParticles) were simplified to use
      the new `catenate()` (O(1) arena-swap when the destination is empty). Audit added
      deepCopy/deepCopyTo/append(bulk)/catenate/swap/shrinkToFit/cellRange + the two-overload append.
- [x] **CD_EBParticleMeshSoA.H** -> `Source/Particle/` (reviewed; now doxygen-checked, passes as-is).
      `CD_EBAMRParticleMeshSoA.H`'s quoted include switched to `<...>`.
- [x] **CD_ParticleContainerSoA.H** + **CD_ParticleContainerSoAImplem.H** -> `Source/Particle/`
      (reviewed; now doxygen-checked, passes as-is). Header's Implem include + CD_EBAMRParticleMeshSoA.H's
      include switched to `<...>`.
- [x] **CD_EBAMRParticleMeshSoA.H** -> `Source/Particle/` (reviewed; now doxygen-checked, passes
      as-is). All its `<...>` dependencies were already in Source. **The SoA library stack is now fully
      migrated to Source/Particle/; only test/benchmark scaffolding remains in Dev/.**

## Support / fixtures (decide on migration vs retire)

- [ ] **CD_DemoParticle.H** (demo payload) -- likely becomes a test fixture, not a library header.
- [ ] Dev tests (`TestParticleContainer*`, `TestEBAMRParticleMesh`, `TestRemapDepositParity`) ->
      `Exec/Tests/...` once their headers are in Source.
- [ ] `Dev/Benchmark`, `Dev/main.cpp` -- demo/benchmark scaffolding; migrate or retire last.

## Consumer rewiring: TracerParticle/TracerParticleSolver -> SoA, then DischargeInception

Decisions (locked): (1) per-consumer NAMED payloads (drop the generic `TracerParticle<M,N>`);
(2) add AmrMesh SoA overloads (don't bypass AmrMesh); (3) retire the AoS `TracerParticle<M,N>` +
AoS `TracerParticleSolver` once all consumers move (no other users -- only TracerParticles,
Physics/TracerParticle + CoaxialCable exec, DischargeInception); (4) OK to break HDF5 checkpoint
compatibility on this branch.

- [x] **Phase 0 -- library prerequisites.**
      - `interpolateWeight` added to `EBParticleMeshSoA` + `EBAMRParticleMeshSoA` (mirror of
        `depositWeight`; scatters into the container-owned weight column). Validated bit-for-bit vs
        the AoS `interpolate<...,&weight>` in `TestRemapDepositParity` (Check D), NGP+CIC, 2D/3D, np1+np2.
      - AmrMesh SoA overloads: `allocate`/`remapToNewGrids`/`depositWeight`/`depositParticles<Members>`/
        `interpolateWeight`/`interpolateParticles<Members>` + `getParticleMeshSoA`. `PhaseRealm` now holds
        `m_particleMeshSoA` (defined alongside `m_particleMesh`); `Realm`/`AmrMesh` expose it. discharge-lib
        compiles; the templated overloads get full runtime coverage via the tracer solver (Phase 2/3).
- [x] **Phase 0.5 -- SoA particle HDF5 checkpoint I/O.** DischargeIO::writeCheckParticlesToHDF /
      readCheckParticlesFromHDF overloads for LayoutData<ParticleSoA> (mirror the AoS ParticleData<P>
      path using the leaf h5* primitives; records are position+weight+payload-h5-subset, all doubles).
- [x] **Phase 1 -- named payload** (TracerParticlePayload = velocity vx/vy/vz + RK scratch xk/k1/k2/k3,
      per-component columns; defined in CD_TracerParticleStepper.H with its ParticleTraits).
- [x] **Phase 2 -- TracerParticleSolver -> SoA** (in place): ParticleContainerSoA + AmrMesh SoA overloads;
      deposit->depositWeight; interpolateWeight; interpolateVelocities->interpolate<&P::vx,...>; computeDt
      SoA loop; checkpoint via the new SoA DischargeIO; getParticles returns ParticleContainerSoA. Also
      added removeCoveredParticlesIF(SoA); ParticleManagement::drawRandomParticles(ParticleSoA) fills a free
      buffer (List<P> analogue) and ParticleContainerSoA::addParticlesDestructive(ParticleSoA) routes it.
- [x] **Phase 3 -- TracerParticleStepper -> SoA** (Euler/RK2/RK4 loops -> SoA column access; seeding via
      drawRandomParticles(SoA)). CoaxialCable exec instantiates TracerParticleStepper<TracerParticlePayload>;
      builds + runs (OPT 2D, plot output advancing).
- [ ] **Phase 4 -- DischargeInceptionStepper -> SoA** (seed/integration loops, gradAlpha interpolation,
      rewind/reset).
- [ ] **Retire** AoS `TracerParticle<M,N>` + AoS `TracerParticleSolver` after Phases 2-4.

## Conventions for each migration

1. `git mv Dev/CD_X.H Source/Particle/CD_X.H` (preserve history).
2. Change internal `#include "CD_Y.H"` -> `#include <CD_Y.H>` (Source style).
3. Remove `./Dev/CD_X.H` from `Docs/doxygen.conf` EXCLUDE.
4. Verify: `pre-commit run clang-format/reuse/doxygen-check`, and rebuild a consumer.
