## Exec/Examples/ItoKMC/ComparisonCdrPlasma

This example runs the *same* 2D streamer discharge problem with two different plasma models:

* A fluid (drift-diffusion) model, `Physics/CdrPlasma`, using `CdrCTU` for transport and `McPhoto` for radiative transfer.
* A particle model, `Physics/ItoKMC`, using `ItoSolver` for transport, `McPhoto` for radiative transfer, and a kinetic
  Monte Carlo description of the plasma chemistry.

The geometry, the mesh, the applied voltage, the transport data, the plasma chemistry, and the initial particles are the
same for both models, so the two runs differ only in how the plasma itself is represented.
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
* **Reactions.** Both models use N2 and O2 ionization, dissociative and three-body attachment to O2-, photon generation
  from N2 with efficiency 0.06 and a quenching pressure of 4000 Pa, and photoionization of O2. There is no secondary
  emission in either model.
* **No gradient correction.** The `soloviev` field is absent from every reaction in `cdr_chemistry.json`, i.e. the
  density-gradient correction to the LFA rates is off, so the rates are the same functions of E in both models.
* **Initial particles.** Both models are seeded with 10^5 computational particles of weight 10^4 (10^9 physical
  electrons and as many O2+), drawn from a Gaussian distribution of width 200 um centred 1 mm below the needle tip.
  Both files use `ParticleManagement::drawGaussianParticles` through their respective JSON interfaces, and with
  `Random.seed = 0` the two models therefore start from the same particles. The Ito solver keeps them as particles; the
  CDR solver deposits them on the mesh.

# What necessarily differs

* `ItoKMCJSON` derives the Townsend coefficients from the reaction list (`"type": "auto"`), which `CdrPlasmaJSON` cannot
  do, so `cdr_chemistry.json` uses the tabulated BOLSIG+ alpha/N and eta/N instead. These coefficients only enter the
  cell tagger, not the transport or the reactions.
* The time integrators are different classes (`CdrPlasmaGodunovStepper` and `ItoKMCGodunovStepper`) and therefore have
  their own option sections. Both are set to a semi-implicit field coupling with a maximum time step of 100 ps.

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
