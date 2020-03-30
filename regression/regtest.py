import os
import argparse
import sys
import configparser
sys.path.append('./python')

tests_file       = "tests.ini"
regression_rules = "regression_rules.py"

# Arguments for this script
parser = argparse.ArgumentParser()
parser.add_argument('-all', '--all',          help="Run all tests", action='store_true', )
parser.add_argument('-compile', '--compile',  help="Compile executables", action='store_true')
parser.add_argument('-procs',                 help="Number of processors to use for compiling and running", type=int, default=1)
parser.add_argument('--benchmark',            help="Generate benchmark files only", action='store_true')
parser.add_argument('-tests',                 help="Run one or more regression tests", nargs='+', required=False)
parser.add_argument('--silent',               help="Turn off unnecessary output", action='store_true')
parser.add_argument('--clean',                help="Do a clean compile", action='store_true')
parser.add_argument('-dim',                   help="Regression suite dimensions", type=int, default=2)
parser.add_argument('-run',                   help="MPI run command", type=str, default="mpirun")

# Read arguments and configuration files
args = parser.parse_args()
config = configparser.ConfigParser()
config.read(tests_file)
baseDir = os.getcwd()

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


def compile(args):
    if args.compile:
        makeCommand = "make "
        if args.silent:
            makeCommand += "-s "
        makeCommand += "-j" + str(args.procs) + " "
        makeCommand += "DIM=" + str(args.dim) + " "
        if args.clean:
            makeCommand += "clean "
        makeCommand += "main"

        print("\t Compiling with = '" + str(makeCommand) + "'\n")
        os.system(makeCommand)

# Do a sanity check before trying tests. 
sanityCheck()


# Run all tests
for test in config.sections():

    directory  = str(config[str(test)]['directory'])
    executable = str(config[str(test)]['exec']) + str(args.dim) + "d.*.ex"
    input      = str(config[str(test)]['input']) + str(args.dim) + "d.inputs"
    nplot      = int(config[str(test)]['plot_interval'])
    nsteps     = int(config[str(test)]['nsteps'])
    if args.benchmark:
        output     = str(config[str(test)]['benchmark'])
    else:
        output     = str(config[str(test)]['output'])

    
    print("\n")
    print("Running regression test '" + str(test) + "' with DIM=" + str(args.dim))
    if not args.silent:
        if args.benchmark:
            print("\t Running benchmark!")
            print("\t Directory is  = " + directory)
            print("\t Input file is = " + input)
            print("\t Output files are = " + str(output) + ".stepXXXXXXX." + str(args.dim) + "d.hdf5")

    # Change to test directory and compile if necessary
    os.chdir(baseDir + "/" + directory) # Change to test directory
    compile(args)                       # Compile if called for

    # Run the executable
    runCommand = args.run + " -np " + str(args.procs) + " " + executable + " " + input
    runCommand = runCommand + " driver.output_names=" + str(output)
    runCommand = runCommand + " driver.plot_interval=" + str(nplot)
    runCommand = runCommand + " driver.max_steps=" + str(nsteps)
    print("\t Executing with '" + str(runCommand) + "'")
    os.system(runCommand)

    # Output at the end
    if args.benchmark:
        print("Regression test '" + str(test) + "' has generated benchmark files \n")
    else:
        print("\t Comparing files for test '" + str(test) + "'\n")

        for i in range (0, nsteps+nplot, nplot):
            regFile =  "plt/" + str(config[str(test)]['output'])
            benFile =  "plt/" + str(config[str(test)]['benchmark'])

            regFile = regFile + (".step{0:07}.".format(i)) + str(args.dim) + "d.hdf5"
            benFile = benFile + (".step{0:07}.".format(i)) + str(args.dim) + "d.hdf5"
            
            print("\t Comparing files " + regFile +  " and " + str(benFile))
            print(os.getcwd())
            os.system("cmp -l " + regFile + " " + benFile)
