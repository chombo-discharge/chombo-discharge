# Regression
This folder can run various regression or unit tests.

## Main usage
Use various flags for controlling how the applications are run and whether or not benchmark files for comparison will be generated.
Configurable options are

* ```--compile``` Force application to compile first.
* ```--silent``` Turn off compiler and run-time messages.
* ```--clean``` Do a clean compilation.
* ```--benchmark``` Generate benchmark files.
* ```--no_exec``` Compile, but do not run tests.
* ```--compare``` Run, and then compare with benchmark files.
* ```-cores``` Run with specified number of cores. If using MPI+OpenMP this will use 1 MPI rank. 
* ```-suites``` Tests suites to run. E.g. ```-suite AdvectionDiffusion Electrostatics```.
* ```-tests``` Specific tests to run. E.g. ```-tests AdvectionDiffusion/Godunov```.

By default, the compilation is done using the configured options from Make.defs.local.

```shell
python tests.py --compile --clean --silent --no_exec -cores 12
```

will compile all the tests using 12 cores, but not run them.
Benchmark files can then be generated with

```shell
python tests.py -cores 12 --benchmark
```

After adding some new code to chombo-discharge, one can compare against the benchmark files by

```shell
python3 tests.py --compile --silent --compare -cores 12
```

## Running a specific suite
To run a specific suite of tests do
```shell
python3 tests.py -suite AdvectionDiffusion
```

## Adjusting compilation parameters
The user can override compilation parameters that are defined in Make.defs.local as follows:

* ```-mpi``` Use MPI. E.g., -mpi=true will enable MPI. 
* ```-hdf``` Use HDF5. Using -hdf=true will enable HDF5. 
* ```-openmp``` Use OpenMP. Using -openmp=true will enable OpenmP.
* ```-petsc``` Use OpenMP. Using -openmp=true will enable OpenmP. 
* ```-dim```Test dimensionality. E.g. ```-dim 2``` does all 2d tests. If ```-dim``` is not 2 or 3, both 2d and 3d are run. 

E.g.,

```shell
python3 tests.py --compile --clean --silent --benchmark -dim=2 -mpi=true -hdf=true -openmp=true -petsc=true -cores 12
```

will do a clean compile of all tests and generate benchmark files (using MPI, HDF5, PETSC, and OpenMP).
