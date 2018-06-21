#!/usr/bin/python
import argparse
import os
import sys
sys.path.append('./app_builder')
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
parser.add_argument('-chombo_home',     type=str,  help="Chombo source code base directory", default=os.environ.get('CHOMBO_HOME', os.getcwd()))
parser.add_argument('-streamer_home',   type=str,  help="Source code base directory", default=os.environ.get('STREAMER_HOME', os.getcwd()))
parser.add_argument('-base_dir',        type=str,  help="Base directory of mini-app", default="./mini_apps")
parser.add_argument('-app_name',        type=str,  help="Mini app name. An error message is issued if the name already exists")
parser.add_argument('-filename',        type=str,  help="File name of main file", default="main")
parser.add_argument('-plasma_kinetics', type=str,  help="Plasma kinetics class", default="")
parser.add_argument('-geometry',        type=str,  help="Geometry class", default="")
parser.add_argument('-time_stepper',    type=str,  help="Time stepping method", default="rk2")
parser.add_argument('-cell_tagger',     type=str,  help="Cell tagging method", default="field_tagger")
args = parser.parse_args()


app_main.write_template(args)
app_make.write_template(args)
app_options.write_template(args)

if args.build:
    os.chdir(args.streamer_home + "/" + args.base_dir + "/" + args.app_name)
    os.system('pwd')
    if args.silent:
        os.system('make -s -j ' + str(args.procs) + ' ' + args.filename)
    else:
        os.system('make -j ' + str(args.procs) + ' ' + args.filename)
    print 'Created and built your mini app - it resides in ' + args.streamer_home + "/" + args.base_dir + "/" + args.app_name
else:
    print 'Created (but did not build) your mini app - it resides in ' + args.streamer_home + "/" + args.base_dir + "/" + args.app_name
