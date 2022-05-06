## Convergence/Poisson

This example computes convergence rates for MFHelmholtzOp, a Helmholtz operator with support for discontinuous material coefficients.
The example uses a "coaxial cable" geometry with an embedded dielectric.

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs```

# Output

The output shows the Linf/L1/L2 errors which are obtained with the help of an exact solution.
The typical output is

Convergence rates: 
|Norm    |Coar err|       |Fine err|       |Conv. rate|
|Inf     |0.00461732      |0.00121248      |1.92909|
|1       |0.00317894      |0.00074213      |2.0988|
|2       |0.00330765      |0.000776896     |2.09001|

|Inf     |0.00121248      |0.000347123     |1.80444|
|1       |0.00074213      |0.000178246     |2.0578|
|2       |0.000776896     |0.000187113     |2.05381|

|Inf     |0.000347123     |9.41944e-05     |1.88173|
|1       |0.000178246     |4.47828e-05     |1.99285|
|2       |0.000187113     |4.7206e-05      |1.98687|

|Inf     |9.41944e-05     |3.82267e-05     |1.30106|
|1       |4.47828e-05     |1.12161e-05     |1.99738|
|2       |4.7206e-05      |1.2025e-05      |1.97293|
