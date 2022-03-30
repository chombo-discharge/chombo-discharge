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
python3 tests.py --compile --clean --silent --benchmark
```

will do a clean compile of all tests and generate benchmark files.

To compare the benchmark files with a new run, do
```shell
python3 tests.py --compile --clean --silent
```
