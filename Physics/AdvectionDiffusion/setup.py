#!/usr/bin/python
import argparse
import os
import sys
sys.path.append('./python')
import app_main
import app_make
import app_options
import app_inc

# Get arguments from input script
parser = argparse.ArgumentParser();
parser.add_argument('-dim',             type=int,  help='Dimension', default=2)
parser.add_argument('-procs',           type=int,  help='Processors to use when building executable', default=1)
parser.add_argument('-use_mpi',         type=bool, help='MPI enabled (default true)', default=True)
parser.add_argument('-build',           type=bool, help='Build executable at end', default=False)
parser.add_argument('-silent',          type=bool, help='Silent build of executable', default=False)
parser.add_argument('-discharge_home',  type=str,  help="Source code base directory", default=os.environ.get('DISCHARGE_HOME'))
parser.add_argument('-base_dir',        type=str,  help="Base directory of mini-app", default="./mini_apps")
parser.add_argument('-app_name',        type=str,  help="Mini app name. An error message is issued if the name already exists")
parser.add_argument('-filename',        type=str,  help="File name of main file", default="main")
parser.add_argument('-cdrsolver',       type=str,  help="CDR solver", default="CdrGodunov")
parser.add_argument('-geometry',        type=str,  help="Geometry class", default="RegularGeometry")
parser.add_argument('-stepper',         type=str,  help="Time stepper class", default="AdvectionDiffusionStepper")

args = parser.parse_args()

# Check if DISCHARGE_HOME and CHOMBO_HOME has been set
discharge_home = args.discharge_home
if not args.discharge_home:
    print "Error: Cannot set application because the DISCHARGE_HOME environment variable has not been set."
    print "       Please set DISCHARGE_HOME, for example:"
    print "       >export  DISCHARGE_HOME=<directory>"
else:
    print "DISCHARGE_HOME is " + args.discharge_home
    print 'Setting up problem in directory ' + args.discharge_home + "/" + args.base_dir + "/" + args.app_name

    app_main.write_template(args)    # Write main file
    app_make.write_template(args)    # Write makefile
    app_options.write_template(args) # Write options file

    # Build executable if called for it
    if args.build:
        os.chdir(args.discharge_home + "/" + args.base_dir + "/" + args.app_name)
        os.system('pwd')
        if args.silent:
            os.system('make -s -j ' + str(args.procs) + ' ' + args.filename)
        else:
            os.system('make -j ' + str(args.procs) + ' ' + args.filename)
            print 'Created and built your mini app - it resides in ' + args.discharge_home + "/" + args.base_dir + "/" + args.app_name
    else:
        print 'Problem setup successful'
