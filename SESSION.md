# SESSION.md — `ito_cdr_deposition` / PR #723

Handoff notes for picking this up cold. Delete this file before merging; it is deliberately its
own commit so `git rebase -i` or a single `git revert` drops it.

## State

| | |
|---|---|
| Branch | `ito_cdr_deposition` |
| PR | [chombo-discharge/chombo-discharge#723](https://github.com/chombo-discharge/chombo-discharge/pull/723), open, ready for review |
| Base | upstream `main` @ `03f64735` (the squash-merge of #722) |
| Content | one commit, 13 files, +709 / -176 |
| Worktree | `.claude/worktrees/ito_cdr_deposition` |

`@claude review` was requested **before** the rebase, so its comments refer to a diff that then had
six commits in it. Consider re-triggering.

### Building here

The worktree has its own `json` and `EBGeometry` submodules but shares the main checkout's
already-built Chombo:

```bash
export DISCHARGE_HOME=/home/robertm/Projects/chombo-discharge/.claude/worktrees/ito_cdr_deposition
export CHOMBO_HOME=/home/robertm/Projects/chombo-discharge/Submodules/Chombo-3.3/lib
```

Do not point `CHOMBO_HOME` at this worktree — its `Submodules/Chombo-3.3` is empty and building it
from scratch is a long detour for a dependency this branch does not touch.

## What the change does

The reaction network handed its two kinds of product to the mesh in different ways. Ito products
became particles at a sub-cell position and reached the mesh through the solver's deposition kernel
(CIC by default). CDR products were added straight into the producing cell — NGP by construction. An
ionization event creating an Ito electron and a CDR ion therefore spread the electron's charge over
a deposition stencil and left the ion's charge in a single cell: a charge separation on the cell
scale, proportional to the ionization rate, so largest at the streamer head.

Photoionization showed both branches side by side. `ItoKMCPhysics::reconcilePhotoionization` appends
the Ito product and the CDR product at the *same* absorption position
(`CD_ItoKMCPhysicsImplem.H`, the `SpeciesType::Ito` / `SpeciesType::CDR` arms), after which the Ito
one was deposited with the solver's kernel and the CDR one with a hardcoded NGP.

New key `<stepper>.cdr_products`:

| value | behaviour |
|---|---|
| `mesh` | Production goes straight into the producing cell. The old behaviour, bit-for-bit. |
| `particle` | Production is emitted as particles and deposited with the Ito solvers' kernel. **Default.** |

### Map of the change

| symbol | file | note |
|---|---|---|
| `CdrProductInjection` (enum) | `CD_ItoKMCStepper.H` | `Mesh` / `Particle`. The only new type in the PR. |
| `m_cdrProductInjection` | `CD_ItoKMCStepper.H` | Set by `parseCdrProducts()`. |
| `m_cdrProducts` | `CD_ItoKMCStepper.H` | **Renamed** from `m_cdrPhotoiProducts`; now carries all CDR production, not just photoionization. |
| `m_particleCdrProduction` | `CD_ItoKMCStepper.H` | Particle-realm copy of the per-cell production. The only new data holder. |
| `depositCdrProducts()` | `CD_ItoKMCStepperImplem.H` | Deposits `m_cdrProducts`, adds to the per-cell change. |
| `reconcileCdrDensities()` | `CD_ItoKMCStepperImplem.H` | Signature changed: takes one *change* field instead of `(new, old)`. |
| `reconcileCdrParticles()` | `CD_ItoKMCPhysicsImplem.H` | Turns per-cell production into particles. |
| `drawNewParticlePosition()` | `CD_ItoKMCPhysicsImplem.H` | The `ParticlePlacement` switch, now shared by the Ito and CDR branches instead of duplicated. |
| `sampleParentPosition()` | `CD_ItoKMCPhysicsImplem.H` | Two overloads; the `Vector` one pools parents across Ito species for CDR products. |

## Decisions already made — do not relitigate

The user settled these explicitly. Reopening them needs a new conversation, not a judgement call.

1. **Particle-mediated, not a mesh smoothing stencil.** A stencil would be cheaper and noise-free
   but would have to re-implement the EB and coarse-fine handling `AmrMesh::depositWeight` already
   provides.
2. **Production only becomes particles; removal stays in-cell.** Removal is proportional to the
   cell-averaged density, which carries no sub-cell information, so there is nothing to be
   consistent about. The network's net change is split by sign for this.
3. **Independent placement** under the same `ParticlePlacement` policy — explicitly *not*
   co-locating each CDR product with its Ito partner.
4. **The quadrature budget is `max_new_particles`**, the same one the Ito products use. A separate
   `max_new_cdr_particles` knob was proposed and rejected **twice** (first "inherit from", then
   "must BE the same as"). Do not reintroduce one.

## The invariant that is easy to get wrong

The CDR per-cell-particle count has a normalization pinned at *both* ends:

- `computeReactiveCdrParticlesPerCell` reads a cut cell as `floor(kappa * phi * V)` — the physical
  mass actually inside the cut cell.
- `reconcileCdrDensities` writes it back as `phi += change / V`, **not** `change / (kappa*V)`
  ("Increment, but don't divide by kappa"). A cut-cell reaction is therefore deliberately damped by
  kappa, and `redistribute_cdr` exists to push the missing `(1-kappa)` share out.

So for the deposition, the cut-cell selector is a **normalization**, not a spreading rule, and is
deliberately *not* inherited from the Ito solver. It is fixed at `IrregularDeposition::NGP`:

- `Mirror` (the ItoSolver default) makes a cut cell hold `n`, not `kappa*n` — a different quantity
  from the one the removal branch and the redistribution are written in.
- `Native` lets a CIC cloud leak into the covered region, where it is zeroed. Silent mass loss.
- `NGP` keeps `kappa*n` **and** loses nothing.

The deposition *kernel* (NGP/CIC/TSC) and the *coarse-fine* strategy are inherited from the Ito
solvers; only the cut-cell selector is constrained. Under an NGP kernel, `NGP` and `Native` are the
same thing, which is what makes `cdr_products = mesh` exact.

## Verified

All on `Exec/Examples/ItoKMC/AirBasic` — the mixed example (`e:ito` plus `O2+:cdr`, `O2-:cdr`,
`N2+:cdr`). Most other ItoKMC examples, `PartialDischarge` included, are **all-Ito** and exercise
none of this.

- **`cdr_products = mesh` is bit-for-bit identical to upstream `main`.** `h5diff -c`, zero
  differences, 200 steps on 8 ranks, `Random.seed=12345`, default configuration (CIC/mirror
  deposition, downstream placement, photoionization). Re-established against the post-rebase base,
  not inherited from the earlier measurement.
- `cdr_products = particle` with `ItoSolver.deposition = ngp` plus centroid placement also
  reproduces the old behaviour bit-for-bit — the particle path collapses onto the mesh path in the
  NGP limit. (Measured on the pre-rebase base.)
- `Exec/Tests/ItoKMC/JSON` passes in 2D and 3D, no NaNs or aborts.
- `pre-commit` clean; `python3 -m sphinx -W --keep-going -b html source build/html` clean.

Reproduce the headline check:

```bash
cd $DISCHARGE_HOME/Exec/Examples/ItoKMC/AirBasic
make -j$(nproc) DIM=2 MPI=TRUE OPT=HIGH DEBUG=FALSE USE_HDF=TRUE
mpirun -np 8 ./main2d.*.ex example.inputs \
  Random.seed=12345 Driver.max_steps=200 Driver.plot_interval=200 \
  Driver.output_dt=-1 Driver.restart=0 Driver.output_directory=<dir> \
  ItoKMCGodunovStepper.cdr_products=mesh
```

then `h5diff -c` against the same run built from upstream `main`.

## Not verified — the honest gap

**The physics improvement has not been demonstrated.** A fixed-step-count A/B cannot resolve it:

| run | sim time reached at step 200 |
|---|---|
| old scheme, seed 12345 | `1.079e-08` |
| new scheme, seed 12345 | `9.930e-09` (−8.6% vs old) |
| new scheme, seed 999 | `9.318e-09` (−6.2% vs the same scheme) |

The scheme difference is the same size as pure seed-to-seed scatter. Any stochastic ItoKMC model
diverges chaotically once the RNG streams differ, and they will differ, because creating CDR product
particles consumes draws. Sizing the real effect needs an ensemble over seeds or a purpose-built
deterministic diagnostic. Do not quote the 8.6% as an effect — it is noise.

## Open item

The quadrature is stochastic, so its sampling error falls off as `1/sqrt(N)`. A **deterministic**
alternative is strictly better and was never built: under CIC, placing `2^D` particles at the corners
of a box of half-width `dx/4` about the cell centre reproduces the *expected* uniform-placement
stencil exactly, with zero noise, at 4 particles in 2D / 8 in 3D, and degenerates correctly to NGP.
It would slot in as a third `cdr_products` mode. Worth proposing if the noise turns out to matter.

## Two claims of mine that were wrong

Recorded so they are not resurrected from the PR history:

1. **"This fixes a pre-existing kappa bug in the photoionization path."** It does not. Under the old
   NGP kernel, `Native` and `NGP` are identical, and `NGP + Native + volumeScale` yields exactly the
   PPC unit the CDR path expects. The old path was internally consistent. What is real is the
   Ito-vs-CDR mismatch in cut cells (the Ito partner deposits with `Mirror`) — an instance of the
   inconsistency this PR fixes, not an arithmetic error. The kappa damping is the documented
   convention, not a defect.
2. **"Inheriting `max_new_particles` (16–32) will be prohibitively expensive."** It is not. AirBasic,
   8 ranks, 200 steps: 62.9 s with the inherited budget of 16, against 59.2 s for the pre-change
   binary and 69.3 s for an earlier N=1 run. The N=1 vs N=16 difference is not resolvable against
   run-to-run scatter.

## Gotchas

- `Driver.stop_time=1E99` is safe on AirBasic and PartialDischarge (`max_dt` caps `dt`) but destroys
  the `Exec/Tests/ItoKMC/JSON` run — `dt` blows up to 1e99 via the Hardcap once the plasma decays.
- Lowering `AmrMesh.max_amr_depth` on AirBasic under-resolves the discharge and the Poisson solve
  fails to converge. It looks like a code bug and is not.
- `Random.seed = -1` in the example inputs is clock-seeded. Always override it for any comparison.
- `Driver.output_dt` must be set `-1` or plotting goes time-based and ignores `plot_interval`,
  giving non-comparable file sets between runs.
- The `format-input-files` pre-commit hook rewrites column alignment in `.options`/`.inputs`. Expect
  it to modify files on the first commit attempt after adding or removing a key.
- `backup/ito_cdr_deposition-prerebase` @ `367a9edb` is a local ref holding the pre-rebase shape
  (six commits). Safe to delete.
