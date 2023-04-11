# Regression
This folder can run various regression or unit tests.

Applications have been set up using various modules. To see how to run the tests do

```shell
python3 tests.py --help
```

## Running a specific suite
To run a specific suite of tests do
```shell
python3 tests.py -suite AdvectionDiffusion
```

## Other options
Use various flags for controlling how the applications are run and whether or not benchmark files for comparison will be generated.

Configrable options are

* ```--compile``` Force application to compile first.
* ```--silent``` Turn off compiler and run-time messages.
* ```--clean``` Do a clean compilation.
* ```--benchmark``` Generate benchmark files.
* ```--no_exec``` Compile, but do not run tests.
* ```--compare``` Run, and then compare with benchmark files.
* ```-dim```Test dimensionality. E.g. ```-dim 2``` does all 2d tests. If ```-dim``` is not 2 or 3, both 2d and 3d are run. 
* ```-mpi``` (true/false) Use MPI or not
* ```-hdf``` (true/false) Use HDF5 or not
* ```-cores``` Run with specified number of cores.
* ```-suites``` Tests suites to run. E.g. ```-suite AdvectionDiffusion Electrostatics```.
* ```-tests``` Specific tests to run. E.g. ```-tests AdvectionDiffusion/Godunov```.

E.g.,

```shell
python3 tests.py --compile --clean --silent --benchmark -dim=2 -mpi=true -hdf=true -cores 12
```

will do a clean compile of all tests and generate benchmark files (using MPI and HDF5).
If one wants to compare previous output with benchmark files, one may later do

```shell
python3 tests.py --compile --silent --compare -dim=2 -mpi=true -hdf=true -cores 12
```
