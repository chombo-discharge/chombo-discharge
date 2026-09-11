## Exec/Examples/KineticMonteCarlo/StiffTwoGroup

This example runs the Kinetic Monte Carlo solver on a stiff two-group electron chemistry and compares
the mean population against the exact solution. It is the reproducer for
[issue #728](https://github.com/chombo-discharge/chombo-discharge/issues/728), where the hybrid
tau-leaping algorithms lost reactions when one reaction turned its reactant over faster than the step.

The model is two electron groups, a slow branching reaction feeding a fast one that relaxes back:

```
e  -> el + el + M+     nu_i = 7.468e9  1/s
el -> e                nu_r = 1.741e11 1/s     (23x faster)
```

All rates are constants, so no transport data is needed and the field is irrelevant. Starting from
`initial_particles` electrons in `e`, the total `e + el` grows at the dominant eigenvalue of
`[-nu_i, nu_r; 2 nu_i, -nu_r]` = `6.8987e9 1/s`. Every reaction is first order, so the moment
hierarchy closes and the mean populations are known exactly (**993.89** per initial electron at
`t = 1 ns`). The program reads the two rates from `chemistry.json` and integrates the moment
equations itself, so the reference always matches the chemistry that was run.

# Compilation

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 main```

# Running the example

```mpirun -np <num_proc> main2d.*ex example.inputs```

Every rank prints the same table of means over all realizations on all ranks to its `pout.*` file,
one row per step handed to `advanceKMC`, with the columns

```
time  <e>  <el>  <e+el>  stderr(<e+el>)  e_moments  el_moments  total_moments  deviation(%)
```

where the `_moments` columns are the solution of the moment equations and `deviation` is the relative
deviation of `<e+el>` from it. The last line repeats the deviation at `stop_time`. With 12 ranks,
`num_runs = 20000` and one initial electron the sampling error on the mean is about 0.2 %. To plot the
KMC solution against the moment solution, e.g.

```
grep -E '^[0-9]' pout.0 > kmc.dat
gnuplot -p -e "set log y; plot 'kmc.dat' u 1:4 w p t 'KMC', '' u 1:8 w l t 'moments'"
```

# What it shows

The chemistry is only mildly stiff (a rate ratio of 23), but with `max_dt = 1e-11` the relaxation
reaction turns the `el` group over `nu_r * dt = 1.74` times per step. A tau-leap that bounds only the
net change in each population does not see this, because production and loss of `el` nearly cancel,
and the Poisson-sampled firings then routinely exceed the population they draw on. The solver limits
the leap by the gross consumption of every reactant so that the fast reaction is resolved instead.

The algorithm and the solver parameters can be varied in the input script or on the command line,
e.g. `ItoKMCJSON.algorithm=ssa` or `initial_particles=1000`. `ssa` is exact and serves as a
reference; the second-order `hybrid_midpoint` and `hybrid_prc` agree with it to within the sampling
error, while the first-order Euler propagators carry a truncation error that is linear in the leap
and hence in `prop_eps`.
