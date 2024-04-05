import os
import glob
import argparse
import sys
import configparser
import subprocess
import time
from subprocess import DEVNULL

# This script requires Python3.5 to work properly
MIN_PYTHON = (3,5)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)
    
# --------------------------------------------------
# Set up arguments that can be passed into this
# script
# --------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('--compile',              help="Compile executables.",           action='store_true')
parser.add_argument('--silent',               help="Turn off unnecessary output.",   action='store_true')
parser.add_argument('--clean',                help="Do a clean compile.",            action='store_true')
parser.add_argument('--benchmark',            help="Generate benchmark files only.", action='store_true')
parser.add_argument('--no_exec',              help="Do not run executables.",        action='store_true')
parser.add_argument('--compare',              help="Turn off HDF5 comparisons",      action='store_true')
parser.add_argument('--mpi',                  help="Use MPI or not",                 action='store_true')
parser.add_argument('--hdf',                  help="Use HDF5 or not",                action='store_true')
parser.add_argument('--openmp',               help="Use OpenMP or not",              action='store_true')
parser.add_argument('-dim',                   help="Test dimensionality", type=int, default=-1, required=False)
parser.add_argument('-cores',                 help="Number of cores to use", type=int, default=2)
parser.add_argument('-exec_mpi',              help="MPI run command.", type=str, default="mpirun")
parser.add_argument('-suites',                help="Test suite (e.g. 'geometry' or 'field')", nargs='+', default="all")
parser.add_argument('-tests',                 help="Individual tests in test suite.", nargs='+', required=False)


args = parser.parse_args()

# --------------
# Base directory.
# --------------
baseDir = os.getcwd()

# ---------
# Exit code
# ---------
ret_code = 0


# -------------------------------------------------------
# Set up test suite. If 'suite' is all, do all test files
# -------------------------------------------------------
test_files = []
if args.suites == "all":
    for file in glob.glob("*.ini"):
        test_files.append(file)
else:
    for s in args.suites:
        file = str(s) + ".ini"
        test_files.append(file)

# --------------------------------
# Print which suites we're running
# --------------------------------
print("Running test suites: " + str(test_files))


# ---------------
# Parse ini files
# ---------------
config = configparser.ConfigParser()
config.read(test_files)

# --------------------------------------------------
# Moron check for running the test suite
# --------------------------------------------------
def pre_check(silent):
    """ Check that chombo-discharge has been appropriately set up with an environment variable. 
        Print some error messages and what to do if we can't run anything. """
    print("Running " + __file__ + "...")
    discharge_home = os.environ.get("DISCHARGE_HOME")
    if not discharge_home:
        print("Error: Cannot run regression tests because the DISCHARGE_HOME environment variable has not been set.")
        print("Please set DISCHARGE_HOME, for example: '> export  DISCHARGE_HOME=<directory>'")
        print("Aborting regtest suite")
        exit()

    print("CWD            = " + os.getcwd())
    print("DISCHARGE_HOME = " + discharge_home)

# --------------------------------------------------
# Function that compiles a test
# --------------------------------------------------
def compile_test(silent, build_procs, dim, mpi, omp, hdf, clean, main):
    """ Set up and run a compilation of the target test. """

    makeCommand = "make "
    if silent:
        makeCommand += "-s "
    makeCommand += "-j" + str(build_procs) + " "
    makeCommand += "DIM=" + str(dim) + " "
    makeCommand += "MPI=" + str(mpi).upper() + " "
    makeCommand += "OPENMPCC=" + str(omp).upper() + " "    
    makeCommand += "USE_HDF=" + str(hdf).upper() + " "
    if omp:
        makeCommand += "USE_MT=FALSE "
    if clean:
        makeCommand += "clean "
    makeCommand += str(config[str(test)]['exec'])        

    print("\t Compiling with   = '" + str(makeCommand) + "'")

    if args.silent:
        exit_code = subprocess.call(makeCommand, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    else:
        exit_code = subprocess.call(makeCommand, shell=True)
    return exit_code

# --------------------------------------------------
# Do a sanity check before trying tests.
# --------------------------------------------------
pre_check(args.silent)


# --------------------------------------------------
# Run all tests
# --------------------------------------------------
for test in config.sections():

    # --------------------------------------------------
    # Check that the configuration parser has the
    # appropriate keys
    # --------------------------------------------------
    do_test = True
    if not config.has_option(str(test), 'directory'):
        print(tests_file + " does not contain option [" + str(test) + "][directory]. Skipping this test")
        do_test = False
    if not config.has_option(str(test), 'exec'):
        print(tests_file + " does not contain option [" + str(test) + "][exec]. Skipping this test")
        do_test = False
    if not config.has_option(str(test), 'input'):
        do_test = False
        print(tests_file + " does not contain option [" + str(test) + "][input]. Skipping this test")
    if not config.has_option(str(test), 'output'):
        do_test = False
        print(tests_file + " does not contain option [" + str(test) + "][output]. Skipping this test")
    if not config.has_option(str(test), 'benchmark'):
        do_test = False
        print(tests_file + " does not contain option [" + str(test) + "][benchmark]. Skipping this test")
    if not config.has_option(str(test), 'nsteps'):
        do_test = False
        print(tests_file + " does not contain option [" + str(test) + "][nsteps]. Skipping this test")
    if not config.has_option(str(test), 'plot_interval'):
        do_test = False
        print(tests_file + " does not contain option [" + str(test) + "][plot_interval]. Skipping this test")
    if not config.has_option(str(test), 'restart'):
        do_test = False
        print(tests_file + " does not contain option [" + str(test) + "][restart]. Skipping this test")
    if not config.has_option(str(test), 'dim'):
        do_test = False
        print(tests_file + " does not contain option [" + str(test) + "][dim]. Skipping this test")        

    # --------------------------------------------------
    # If we're not running all tests, check that the test string equals args.test
    # --------------------------------------------------
    if not args.tests:
        do_test = True
        dim     = int(config[str(test)]['dim'])

        if dim != args.dim and (args.dim==2 or args.dim==3):
            do_test = False
    else:
        do_test = False
        if str(test) in args.tests:
            do_test = True
        
    # --------------------------------------------------
    # If moron check passed, try to run the test
    # --------------------------------------------------
    if do_test:
        # --------------------------------------------------
        # Get test suite parameters from .ini file and
        # convert them to the types that they represent. 
        # --------------------------------------------------
        directory  = str(config[str(test)]['directory'])
        dim        = int(config[str(test)]['dim'])    
        executable = str(config[str(test)]['exec']) + str(dim) + "d.*.ex"
        input      = str(config[str(test)]['input'])
        nplot      = int(config[str(test)]['plot_interval'])
        nsteps     = int(config[str(test)]['nsteps'])
        restart    = int(config[str(test)]['restart'])
        if args.benchmark:
            output     = str(config[str(test)]['benchmark'])
        else:
            output     = str(config[str(test)]['output'])

        start = time.time()        
        
        # --------------------------------------------------
        # Print some information about the regression test 
        # being run. 
        # --------------------------------------------------
        print("\nRunning regression test '" + str(test) + "'...")
        if args.benchmark:
            print("\t Running benchmark!")
            
        print("\t Directory is     = " + directory)
        print("\t Input file is    = " + input)
        print("\t Output files are = " + str(output) + ".stepXXXXXXX." + str(dim) + "d.hdf5")

        # --------------------------------------------------
        # Now change to test directory
        # --------------------------------------------------
        os.chdir(baseDir + "/" + directory) 

        # --------------------------------------------------
        # Check if executable exists. Recompile test if
        # user has called for it
        # --------------------------------------------------
        run_suite = True
        if args.compile:
            compile_code = compile_test(silent=args.silent,
                                        build_procs=args.cores,
                                        dim=dim,
                                        mpi=args.mpi,
                                        omp=args.openmp,
                                        hdf=args.hdf,
                                        clean=args.clean,
                                        main = str(config[str(test)]['exec']))
            if compile_code != 0:
                print("\t Compilation of test '" + str(test) + "' failed. Aborting this regression test.")
                run_suite = False
                ret_code = 1

        # --------------------------------------------------
        # Run the test suite
        # --------------------------------------------------
        if run_suite:
            if args.no_exec:
                print(f'\t Test completed   = {time.time() - start:.2f}s')                
            else:
                # --------------------------------------------------
                # Set up the run command
                # --------------------------------------------------
                if str(args.mpi).upper() == "TRUE":
                    if str(args.openmp).upper() == "TRUE":
                        runCommand = "OMP_NUM_THREADS=" + str(args.cores) + " "
                        runCommand += "OMP_PLACES=cores "
                        runCommand += args.exec_mpi + " -np 1" + " ./" + executable + " " + input                        
                    else:
                        runCommand = args.exec_mpi   + " -np " + str(args.cores) + " ./" + executable + " " + input
                else:
                    runCommand = "./" + executable + " " + input                

                runCommand = runCommand + " Driver.output_names="  + str(output)
                runCommand = runCommand + " Driver.plot_interval=" + str(nplot)
                runCommand = runCommand + " Driver.checkpoint_interval=" + str(nplot)
                runCommand = runCommand + " Driver.max_steps="     + str(nsteps)
                    
                if args.benchmark:
                    runCommand = runCommand + " Driver.restart=0"
                else:
                    runCommand = runCommand + " Driver.restart="     + str(restart)

                print("\t Executing with   = '" + str(runCommand) + "'")

                # --------------------------------------------------
                # Run the executable and print the exit code
                # --------------------------------------------------
                #        exit_code = os.system(runCommand)
                if args.silent:
                    exit_code = subprocess.call(runCommand, shell=True, stdout=DEVNULL, stderr=DEVNULL)
                else:
                    exit_code = subprocess.call(runCommand, shell=True)


                if exit_code != 0:
                    print("\t Test run failed with exit code = " + str(exit_code))
                    ret_code = 2
                else:
                    print(f'\t Test completed   = {time.time() - start:.2f}s')                    
                    # --------------------------------------------------
                    # Do file comparison if the test ran successfully
                    # --------------------------------------------------
                    if args.benchmark:
                        print("\t Regression test '" + str(test) + "' has generated benchmark files.")
                    elif not args.benchmark and args.compare and str(args.hdf).upper == "TRUE":
                        # --------------------------------------------------
                        # Loop through all files that were generated and
                        # compare them with h5diff. Print an error message
                        # if files don't match. 
                        # --------------------------------------------------
                        for i in range (0, nsteps+nplot, nplot):
                        
                            # --------------------------------------------------
                            # Get the two files that will be compared
                            # --------------------------------------------------
                            regFile =  "plt/" + str(config[str(test)]['output'])
                            benFile =  "plt/" + str(config[str(test)]['benchmark'])
                            
                            regFile = regFile + (".step{0:07}.".format(i)) + str(dim) + "d.hdf5"
                            benFile = benFile + (".step{0:07}.".format(i)) + str(dim) + "d.hdf5"

                            if not os.path.exists(benFile):
                                print("\t Benchmark file(s) not found, generate them with --benchmark")
                            else:
                                print("\t Comparing files " + regFile +  " and " + str(benFile))
                                
                                # --------------------------------------------------
                                # Run h5diff and compare the two files. Print a
                                # petite message if they match, and a huge-ass 
                                # warning if they don't. 
                                # --------------------------------------------------
                                compare_command = "h5diff " + regFile + " " + benFile
                                if args.silent:
                                    compare_code = subprocess.call(compare_command, shell=True, stdout=DEVNULL, stderr=DEVNULL)
                                else:
                                    compare_code = subprocess.call(compare_command, shell=True)
                                    
                                    if compare_code != 0:
                                        print("\t FILES '" + regFile +  "' AND '" + benFile + "' DO NOT MATCH - REGRESSION TEST FAILED")
                                    else:
                                        print("\t Benchmark test succeded for files " + regFile +  " and " + str(benFile))

exit(ret_code)
