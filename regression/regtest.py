import os
import argparse
import sys
from glob import glob
sys.path.append('./python')


# Get arguments from input script
parser = argparse.ArgumentParser()
parser.add_argument('-all', '--all',          help="Run all tests", action='store_true', )
parser.add_argument('-compile', '--compile',  help="Compile executables", action='store_true')
parser.add_argument('-procs',                 help="Number of processors to use for compiling and running", type=int, default=1)
parser.add_argument('--benchmark',            help="Generate benchmark files only", action='store_true', )
parser.add_argument('-tests',                 help="Run one or more regression tests", nargs='+', required=False)

args = parser.parse_args()



def plasmac_sanity():
    """ Check that PLASMAC has been set as appropriate"""
    print("Running " + __file__ + "...")
    
    plasmac_home = os.environ.get("PLASMAC_HOME")
    if not plasmac_home:
        print("Error: Cannot run regression tests because the PLASMAC_HOME environment variable has not been set.")
        print("Please set PLASMAC_HOME, for example:")
        print(">export  PLASMAC_HOME=<directory>")
        print("Aborting regtest suite")
        exit()
    else:
        print("CWD          = " + os.getcwd())
        print("PLASMAC_HOME = " + plasmac_home)


# Check that PLASMAC_HOME variable has been set
plasmac_sanity()

# Set up tests to run


# Set up test suite
# print "Input arguments = "
# print args

# print os.walk(os.getcwd())

# print glob(os.getcwd()+"/*/")

folder=os.getcwd()
subfolders = [ f.path for f in os.scandir(folder) if f.is_dir() ]

print (subfolders)
