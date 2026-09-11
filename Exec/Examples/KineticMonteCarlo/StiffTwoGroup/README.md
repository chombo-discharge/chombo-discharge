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
a single `e`, the total `e + el` grows at the dominant eigenvalue of `[-nu_i, nu_r; 2 nu_i, -nu_r]`
= `6.8987e9 1/s`. Every reaction is first order, so the moment hierarchy closes and the mean is
known exactly: **993.89** at `t = 1 ns`. This value is passed to the program as `exact_mean` and must
be updated if `stop_time` or the rates are changed.

# Compilation

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 main```

# Running the example

```mpirun -np <num_proc> main2d.*ex example.inputs```

The program prints the mean total electron number over all realizations on all ranks, the exact mean,
and the relative deviation, to `pout.0`. With 12 ranks and `num_runs = 20000` the sampling error on
the mean is about 0.2 %.

# What it shows

The chemistry is only mildly stiff (a rate ratio of 23), but with `max_dt = 1e-11` the relaxation
reaction turns the `el` group over `nu_r * dt = 1.74` times per step. A tau-leap that bounds only the
net change in each population does not see this, because production and loss of `el` nearly cancel,
and the Poisson-sampled firings then routinely exceed the population they draw on. The solver limits
the leap by the gross consumption of every reactant so that the fast reaction is resolved instead.

The algorithm and the solver parameters can be varied in the input script. `ssa` is exact and
serves as a reference; the `hybrid_*` algorithms should agree with it to within the sampling error
for any `max_dt`.
