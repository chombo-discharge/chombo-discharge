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

## Conventions for each migration

1. `git mv Dev/CD_X.H Source/Particle/CD_X.H` (preserve history).
2. Change internal `#include "CD_Y.H"` -> `#include <CD_Y.H>` (Source style).
3. Remove `./Dev/CD_X.H` from `Docs/doxygen.conf` EXCLUDE.
4. Verify: `pre-commit run clang-format/reuse/doxygen-check`, and rebuild a consumer.
