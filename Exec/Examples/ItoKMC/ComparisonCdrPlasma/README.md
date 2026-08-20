## Exec/Examples/ItoKMC/ComparisonCdrPlasma

This example runs the *same* 2D streamer discharge problem with two different plasma models:

* A fluid (drift-diffusion) model, `Physics/CdrPlasma`, using `CdrCTU` for transport and `McPhoto` for radiative transfer.
* A particle model, `Physics/ItoKMC`, using `ItoSolver` for the electrons, `McPhoto` for radiative transfer, and a
  kinetic Monte Carlo description of the plasma chemistry.

The comparison is arranged so that the **electron transport model is the only thing that differs**. The three ion
species are CDR solvers in both chemistry files, and are immobile, so they evolve as identical mesh densities on both
sides. Only the electron is carried as Ito particles by the Ito-KMC model and as a CDR density by the CdrPlasma model.
The geometry, the mesh, the applied voltage, the transport data, the plasma chemistry, and the initial seed are shared.

The case is a rod-needle electrode at 35 kV in a grounded 8 cm x 8 cm box of atmospheric air, and is derived from
`Exec/Examples/ItoKMC/AirBasic`.

Files:

* Model switch, mesh, solvers, and time stepping: `example.inputs`
* Plasma chemistry for the fluid model: `cdr_chemistry.json`
* Plasma chemistry for the particle model: `ito_chemistry.json`
* Transport data (BOLSIG+) for both models: `bolsig_air.dat`

# Selecting the model

The model is selected with a single switch in `example.inputs`:

```
Comparison.model = ito_kmc      ## or 'cdr_plasma'
```

Everything else in the input script is shared. The `ItoSolver` options are only read by the Ito-KMC model and the
`CdrCTU` options are only read by the CdrPlasma model; the remaining sections apply to both.

# What is held fixed between the models

* **Transport data.** Both chemistry files read the electron mobility, diffusion coefficient, mean energy, and all
  reaction rates from the same `bolsig_air.dat`, using the same table headers, the same truncation limits, and the same
  number of resampling points.
* **Gas.** Ideal gas at 300 K and 1 atm with 20% O2 and 80% N2. Note that `CdrPlasmaJSON` takes the pressure in
  atmospheres and `ItoKMCJSON` takes it in Pascal, so the two files say `1` and `101325` for the same gas.
* **Ions.** O2+, O2- and N2+ are CDR species in both files, and all three are immobile and non-diffusive. They are
  therefore the same mesh densities evolving under the same chemistry in both models, and contribute nothing to the
  conductivity that enters the semi-implicit field solve.
* **Reactions.** Both models use N2 and O2 ionization, dissociative and three-body attachment to O2-, photon generation
  from N2 with efficiency 0.06 and a quenching pressure of 4000 Pa, and photoionization of O2. There is no secondary
  emission in either model.
* **Townsend coefficients.** Both chemistry files read alpha/N and eta/N directly from the tabulated BOLSIG+
  coefficients in `bolsig_air.dat`. `ItoKMCJSON` could instead derive them from the reaction list (`"type" : "auto"`),
  which `CdrPlasmaJSON` cannot do, so the tabulated form is used in both. These coefficients only enter the cell
  tagger, not the transport or the reactions.
* **No gradient correction.** The `soloviev` field is absent from every reaction in `cdr_chemistry.json`, i.e. the
  density-gradient correction to the LFA rates is off, so the rates are the same functions of E in both models.
* **Refinement criteria.** Both cell taggers refine on `(alpha - eta)*dx` with the same thresholds, so the two runs
  build their grids by the same rule.
* **Initial seed.** Both models are seeded with 10^5 computational particles of weight 10^4 (10^9 physical electrons
  and as many O2+), drawn from a Gaussian distribution of width 200 um centred 1 mm below the needle tip. Both files
  use `ParticleManagement::drawGaussianParticles` through their respective JSON interfaces, and with `Random.seed = 0`
  the two models therefore start from the same particles. The O2+ seed is deposited on the mesh in both models; the
  electron seed is deposited in the CdrPlasma model and kept as particles in the Ito-KMC model.

# What necessarily differs

* **The electron representation**, which is the point of the example. The Ito-KMC model advances electrons as particles
  with an Euler-Maruyama update and resolves the chemistry with a kinetic Monte Carlo solver, so its results carry
  statistical noise and depend on the target particles-per-cell, the superparticle merging algorithm, and the placement
  of secondary particles. None of these have a fluid analogue.
* **The form of the diffusion operator.** `ItoKMCPhysics::isotropicDiffusion` draws `sqrt(2*D*dt)*N(0,1)` and the
  stepper applies `X += dt*v + hop`, which is an Ito-interpretation Euler-Maruyama update and therefore transports the
  flux `v*n - grad(D*n)`. `CdrSolver` integrates the flux `v*n - D*grad(n)`. The two coincide only where `D` is
  constant. Since this example reads the tabulated `D(E/N)`, the models integrate slightly different equations wherever
  `grad(D)` is non-zero, which in atmospheric air is a per-cent-level effect concentrated at the streamer head. Setting
  `D` to a constant in both chemistry files removes the ambiguity if you would rather compare the models without it.
* **The time integrators** are different classes (`CdrPlasmaGodunovStepper` and `ItoKMCGodunovStepper`) and therefore
  have their own option sections. Both are set to a semi-implicit field coupling with a maximum time step of 100 ps,
  but they select their time steps by different rules: `CdrPlasmaGodunovStepper` takes the smaller of the separate
  advective and diffusive limits, while `ItoKMCGodunovStepper` uses the combined advection-diffusion limit and applies
  additional restrictions of its own. Expect the two runs to need a different number of steps to reach the same time.

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2```

This example is 2D only.

# Running the example

```./main2d.*.ex example.inputs```

or with MPI:

```mpirun -np <num_proc> main2d.*.ex example.inputs```

To run the other model without editing the input script:

```./main2d.*.ex example.inputs Comparison.model=cdr_plasma```

# Output

Output is given to HDF5 files in the plt folder. `Driver` writes one plot file every 1 ns and the simulation stops at
5 ns, so both models produce plot files at the same six output times.
