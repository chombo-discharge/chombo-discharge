## Exec/Examples/ItoKMC/ComparisonCdrPlasma

This example runs the *same* 2D streamer discharge problem with two different plasma models, selected by a single
switch, so that the results can be compared directly:

* A fluid model, `Physics/CdrPlasma`, using `CdrCTU` for transport and `McPhoto` for radiative transfer.
* A particle model, `Physics/ItoKMC`, using `ItoSolver` for the electrons, `McPhoto` for radiative transfer, and a
  kinetic Monte Carlo description of the plasma chemistry.

The electrons are the only species treated differently: they are Ito particles in the Ito-KMC model and a CDR density
in the CdrPlasma model. The three ion species are immobile CDR densities in both. The geometry, mesh, applied voltage,
transport data, chemistry, and initial seed are shared.

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

# Notes on the comparison

Both chemistry files read the electron mobility, diffusion coefficient, mean energy, Townsend coefficients, and all
reaction rates from the same `bolsig_air.dat`, and both cell taggers refine on the same criterion with the same
thresholds. Note that `CdrPlasmaJSON` takes the gas pressure in atmospheres while `ItoKMCJSON` takes it in Pascal, so
the two files say `1` and `101325` for the same gas.

The two models will not agree exactly. The Ito-KMC results carry statistical noise and depend on the target
particles-per-cell and the superparticle merging algorithm; the two time steppers also select their time steps by
different rules, so the runs take a different number of steps to reach the same time.

Both models do now transport the same diffusive flux. The particle update carries a `grad(D)` drift correction
(`ItoSolver.diffusion_grad_drift`, on by default) without which an Ito-interpretation update transports
`v*n - grad(D*n)` rather than the `v*n - D*grad(n)` the CDR solvers integrate -- a difference that vanishes only
where `D` is constant, and that made the particle swarm lag the fluid solution by a few percent at the head with
the tabulated `D(E/N)` used here. See
[#706](https://github.com/chombo-discharge/chombo-discharge/issues/706).

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
