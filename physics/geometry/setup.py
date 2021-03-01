#!/usr/bin/python
import argparse
import os
import sys
sys.path.append('./python')
import app_main
import app_make
import app_options

# Get arguments from input script
parser = argparse.ArgumentParser();
parser.add_argument('-dim',             type=int,  help='Dimension', default=2)
parser.add_argument('-procs',           type=int,  help='Processors to use when building executable', default=1)
parser.add_argument('-use_mpi',         type=bool, help='MPI enabled (default true)', default=True)
parser.add_argument('-build',           type=bool, help='Build executable at end', default=False)
parser.add_argument('-silent',          type=bool, help='Silent build of executable', default=False)
parser.add_argument('-plasmac_home',    type=str,  help="Source code base directory", default=os.environ.get('PLASMAC_HOME', os.getcwd()))
parser.add_argument('-base_dir',        type=str,  help="Base directory of mini-app", default="./mini_apps")
parser.add_argument('-app_name',        type=str,  help="Mini app name. An error message is issued if the name already exists")
parser.add_argument('-filename',        type=str,  help="File name of main file", default="main")
parser.add_argument('-geometry',        type=str,  help="Geometry class", default="regular_geometry")

args = parser.parse_args()

app_main.write_template(args)    # Write main file
app_make.write_template(args)    # Write makefile
app_options.write_template(args) # Write options file

# Check if PLASMAC_HOME and CHOMBO_HOME has been set
plasmac_home = args.plasmac_home
if not args.plasmac_home:
    print "Error: Cannot set application because the PLASMAC_HOME environment variable has not been set."
    print "       Please set PLASMAC_HOME, for example:"
    print "       >export  PLASMAC_HOME=<directory>"
else:
    print "PLASMAC_HOME is " + args.plasmac_home
    print 'Setting up problem in directory ' + args.plasmac_home + "/" + args.base_dir + "/" + args.app_name

    app_main.write_template(args)    # Write main file
    app_make.write_template(args)    # Write makefile
    app_options.write_template(args) # Write options file

    # Build executable if called for it
    if args.build:
        os.chdir(args.plasmac_home + "/" + args.base_dir + "/" + args.app_name)
        os.system('pwd')
        if args.silent:
            os.system('make -s -j ' + str(args.procs) + ' ' + args.filename)
        else:
            os.system('make -j ' + str(args.procs) + ' ' + args.filename)
            print 'Created and built your mini app - it resides in ' + args.plasmac_home + "/" + args.base_dir + "/" + args.app_name
    else:
        print 'Problem setup successful'
