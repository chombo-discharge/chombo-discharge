import os
import argparse
import sys
from glob import glob
sys.path.append('./python')

regression_rules = "regression_rules.py" 

# Get arguments from input script
parser = argparse.ArgumentParser()
parser.add_argument('-all', '--all',          help="Run all tests", action='store_true', )
parser.add_argument('-compile', '--compile',  help="Compile executables", action='store_false')
parser.add_argument('-procs',                 help="Number of processors to use for compiling and running", type=int, default=1)
parser.add_argument('--benchmark',            help="Generate benchmark files only", action='store_true')
parser.add_argument('-tests',                 help="Run one or more regression tests", nargs='+', required=False)
parser.add_argument('--silent',               help="Turn off unnecessary output", action='store_true')
parser.add_argument('--clean',                help="Do a clean compile", action='store_true')

args = parser.parse_args()

def sanityCheck():
    """ Check that PLASMAC has been set as appropriate"""
    print("Running " + __file__ + "...")
    plasmac_home = os.environ.get("PLASMAC_HOME")
    if not plasmac_home:
        print("Error: Cannot run regression tests because the PLASMAC_HOME environment variable has not been set.")
        print("Please set PLASMAC_HOME, for example: '> export  PLASMAC_HOME=<directory>'")
        print("Aborting regtest suite")
        exit()
    else:
        print("CWD          = " + os.getcwd())
        print("PLASMAC_HOME = " + plasmac_home)

def setupTests(args):
    """ Check that the regression tests have a regression rules file """
    sanityCheck();
    RegTestDirectory = os.getcwd()
    if args.all:
        tests = [ f.name for f in os.scandir(RegTestDirectory) if f.is_dir()]
    else:
        tests = args.tests
    return tests

def compile(args):
    if args.compile:
        makeCommand = "make "
        if args.silent:
            makeCommand += "-s "
        makeCommand += "-j" + str(args.procs) + " "
        if args.clean:
            makeCommand += "clean "
        makeCommand += "main"
        os.system(makeCommand)


baseDir = os.getcwd()
regressionTests = setupTests(args)
for test in regressionTests:
    os.chdir(baseDir + "/" + str(test)) # Change to test directory
    compile(args)                       # Compile if called for
    os.system("mpirun -np 6 main2d*.ex regression2d.inputs")

