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
E.g.,

```shell
python3 tests.py --compile --clean --silent --no_compare
```

will do a clean compile of all tests and generate benchmark files.
Configrable options are

* ```--compile``` Force application to compile first.
* ```--silent``` Turn off compiler and run-time messages.
* ```--clean``` Do a clean compilation.
* ```--benchmark``` Generate benchmark files.
* ```--no_exec``` Compile, but do not run tests.
* ```--no_compare``` Run, but do not compare with benchmark files.
* ```--parallel``` Compile with MPI.
* ```-cores``` Run with specifies number of cores.
* ```-suites``` Tests suites to run. E.g. ```-suite AdvectionDiffusion Electrostatics```.
* ```-tests``` Specific tests to run. E.g. ```-tests AdvectionDiffusion/Godunov```.
* ```-dim```Test dimensionality. E.g. ```-dim 2``` does all 2d tests. If ```-dim``` is not 2 or 3, both 2d and 3d are run. 

## Comparing benchmark files
To compare the benchmark files with a new run, do
```shell
python3 tests.py --compile --clean --silent
```
