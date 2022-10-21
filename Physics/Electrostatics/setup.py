#!/usr/bin/python
import argparse
import os
import sys
sys.path.append('./python')
import app_main
import app_options
import app_inc

# Get arguments from input script
parser = argparse.ArgumentParser();
parser.add_argument('-discharge_home',  type=str,  help="Source code base directory (default: %(default)s)", default=os.environ.get('DISCHARGE_HOME', os.getcwd()))
parser.add_argument('-base_dir',        type=str,  help="Base directory of application (default: %(default)s)", default="MyApplications")
parser.add_argument('-app_name',        type=str,  help="Application name (default: %(default)s)", default="MyElectrostaticProblem")
parser.add_argument('-field_solver',    type=str,  help="FieldSolver implementation (default: %(default)s)", default="FieldSolverMultigrid")
parser.add_argument('-geometry',        type=str,  help="Geometry class (default: %(default)s)", default="RegularGeometry")

args = parser.parse_args()

# Check if DISCHARGE_HOME and CHOMBO_HOME has been set
discharge_home = args.discharge_home
if not args.discharge_home:
    print("Error: Cannot set application because the DISCHARGE_HOME environment variable has not been set.")
    print("       Please set DISCHARGE_HOME, for example:")
    print("       >export  DISCHARGE_HOME=<directory>")
else:
    print("DISCHARGE_HOME is " + args.discharge_home)
    print('Setting up problem in directory ' + args.discharge_home + "/" + args.base_dir + "/" + args.app_name)

    app_main.write_template(args)    # Write main file
    app_options.write_template(args) # Write options file

    # Copy the makefile
    os.system('cp ./python/GNUmakefile ' + args.discharge_home + "/" + args.base_dir + "/" + args.app_name + "/GNUmakefile")            

    print('Problem setup successful')
