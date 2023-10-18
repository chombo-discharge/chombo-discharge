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
parser.add_argument('-discharge_home', type=str,  default=os.environ.get('DISCHARGE_HOME', os.getcwd()), help = "Source code base directory (default: %(default)s)")
parser.add_argument('-base_dir',       type=str,  default="MyApplications",                              help="Base directory of mini-app (default: %(default)s")
parser.add_argument('-field_solver',   type=str,  default="FieldSolverMultigrid",                        help="Poisson solver implementation (default: %(default)s)")
parser.add_argument('-cdr_solver',     type=str,  default="CdrGodunov",                                  help="CDR solver implementation (default: %(default)s)")
parser.add_argument('-rte_solver',     type=str,  default="EddingtonSP1",                                help="RTE solver implementation (default: %(default)s)")
parser.add_argument('-physics',        type=str,  default="CdrPlasmaJSON",                               help="Plasma model (default: %(default)s)")
parser.add_argument('-geometry',       type=str,  default="RegularGeometry",                             help="Geometry class (default: %(default)s)")
parser.add_argument('-time_stepper',   type=str,  default="CdrPlasmaGodunovStepper",                     help="Time stepping method (default: %(default)s)")
parser.add_argument('-cell_tagger',    type=str,  default="CdrPlasmaStreamerTagger",                     help="Cell tagging method (default: %(default)s)")
parser.add_argument('-app_name',       type=str,  default="MyApplication",                               help="Mini app name. An error message is issued if the name already exists")

args = parser.parse_args()

# Check if DISCHARGE_HOME and CHOMBO_HOME has been set
discharge_home = args.discharge_home
if not args.discharge_home:
    print("Error: Cannot set application because the DISCHARGE_HOME environment variable has not been set.")
    print("       Please set DISCHARGE_HOME, for example:")
    print("       >export  DISCHARGE_HOME=<directory>")
else:
    print("DISCHARGE_HOME is " + args.discharge_home)
    print('Setting up problem in directory ' + args.discharge_home + '/' + args.base_dir + '/' + args.app_name)

    app_main.write_template(args)    # Write main file
    app_options.write_template(args) # Write options file
    app_inc.copy_dependencies(args)  # Copy depencies

    # # Copy the makefile
    os.system('cp ./python/GNUmakefile ' + args.discharge_home + '/' + args.base_dir + "/" + args.app_name + "/GNUmakefile")        

    print('Problem setup successful')
        
