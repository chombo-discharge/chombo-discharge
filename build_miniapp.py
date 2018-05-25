import argparse
import os

# Get arguments from input script
parser = argparse.ArgumentParser();
parser.add_argument('-dim',             type=int,  help='Dimension', default=2)
parser.add_argument('-procs',           type=int,  help='Processors to use when building executable', default=1)
parser.add_argument('-use_mpi',         type=bool, help='MPI enabled (default true)', default=True)
parser.add_argument('-build',           type=bool, help='Build executable at end', default=False)
parser.add_argument('-silent',          type=bool, help='Silent build of executable', default=False)
parser.add_argument('-chombo_home',     type=str,  help="Chombo source code base directory", default=os.getcwd())
parser.add_argument('-streamer_home',   type=str,  help="Source code base directory", default=os.getcwd())
parser.add_argument('-base_dir',        type=str,  help="Base directory of mini-app", default="./mini_apps")
parser.add_argument('-app_name',        type=str,  help="Mini app name. An error message is issued if the name already exists")
parser.add_argument('-filename',        type=str,  help="File name of main file", default="main")
parser.add_argument('-plasma_kinetics', type=str,  help="Plasma kinetics class", default="")
parser.add_argument('-geometry',        type=str,  help="Geometry class", default="")
parser.add_argument('-time_stepper',    type=str,  help="Time stepping method", default="")
parser.add_argument('-cell_tagger',     type=str,  help="Cell tagging method", default="")
args = parser.parse_args()


# Make sure that every class can be found where they should
geofile = args.streamer_home + "/geometries_prebuilt" + "/" + args.geometry + ".H"
tsfile  = args.streamer_home + "/time_steppers" + "/" + args.time_stepper + ".H"
kinfile = args.streamer_home + "/plasma_models" + "/" + args.plasma_kinetics + "/" + args.plasma_kinetics + ".H"
tagfile = args.streamer_home + "/cell_taggers" +  "/" + args.cell_tagger + ".H"
if not os.path.exists(geofile):
    print 'Could not find ' + geofile
if not os.path.exists(tsfile):
    print 'Could not find ' + tsfile
if not os.path.exists(kinfile):
    print 'Could not find ' + kinfile
if not os.path.exists(tagfile):
    print 'Could not find ' + tagfile

# Create app directory if it does not exist
app_dir = args.base_dir + "/" + args.app_name
if not os.path.exists(app_dir):
    os.makedirs(app_dir)

# Write main file. This should be a separate routine. 
main_filename = app_dir + "/" + args.filename + ".cpp"
mainf = open(main_filename, "w")
mainf.write('#include "plasma_engine.H"\n')
mainf.write('#include "' + args.plasma_kinetics + '.H"\n')
mainf.write('#include "' + args.geometry + '.H"\n')
mainf.write('#include "' + args.time_stepper + '.H"\n')
mainf.write('#include "' + args.cell_tagger + '.H"\n')
mainf.write('#include "ParmParse.H"\n')
mainf.write("\n")

mainf.write("// This is the potential curve (constant in this case). Modify it if you want to.\n")
mainf.write("Real g_potential;\n");
mainf.write("Real potential_curve(const Real a_time){\n")
mainf.write("  return g_potential;\n")
mainf.write("}\n")

mainf.write("\n")
mainf.write("int main(int argc, char* argv[]){\n")

mainf.write("\n")
if args.use_mpi:
    mainf.write("#ifdef CH_MPI\n")
    mainf.write("  MPI_Init(&argc, &argv);\n")
    mainf.write("#endif\n")
    
mainf.write("\n")
mainf.write("  // Build class options from input script and command line options\n")
mainf.write("  char* input_file = argv[1];\n")
mainf.write("  ParmParse pp(argc-2, argv+2, NULL, input_file);")
mainf.write("\n")

mainf.write("\n")
mainf.write("  { // Get potential from input script \n")
mainf.write('    ParmParse pp("' + args.app_name + '");\n')
mainf.write('    pp.get("potential", g_potential);\n')
mainf.write("  }\n")

mainf.write("\n")
mainf.write("  // Set up everything \n")
mainf.write("  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics> (new " + args.plasma_kinetics + "());\n")
mainf.write("  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new " + args.geometry + "());\n")
mainf.write("  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper> (new " + args.time_stepper + "());\n")
mainf.write("  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new " + args.cell_tagger + "());\n")
mainf.write("  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());\n")
mainf.write("  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());\n")
mainf.write("  RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom, compgeom, plaskin, timestepper, amr, tagger));\n")
mainf.write("\n")
mainf.write("  // Run the plasma engine\n")
mainf.write("  engine->set_potential(potential_curve);\n");
mainf.write("  engine->setup_and_run();\n");
mainf.write("\n")

if args.use_mpi:
    mainf.write("#ifdef CH_MPI\n")
    mainf.write("  CH_TIMER_REPORT();\n")
    mainf.write("  MPI_Finalize();\n")
    mainf.write("#endif\n")
mainf.write("}\n")
mainf.close()

# Write makefile. This should be a separate routine.
make_filename = app_dir + "/GNUmakefile"
makef = open(make_filename, "w")
makef.write("# Chombo and chombo-streamer directories \n")
makef.write("CHOMBO_HOME   := " + args.chombo_home + "/lib\n")
makef.write("STREAMER_HOME := " + args.streamer_home + "\n")

# Make rules
makef.write("\n")
makef.write("include $(CHOMBO_HOME)/mk/Make.defs")
makef.write("\n")

makef.write("\n")
makef.write("USE_EB=TRUE\n")
makef.write("USE_MF=TRUE\n")
makef.write("DIM=" + str(args.dim) + "\n")
makef.write("\n")

makef.write("\n")
makef.write("# Base file containing int main()\n")
makef.write("ebase := " + args.filename)
makef.write("\n")

makef.write("\n")
makef.write("LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \\\n")
makef.write("\tAMRTimeDependent BaseTools BoxTools Workshop\n")
makef.write("\n")
 	
makef.write("\n")
makef.write('# Target\n')
makef.write('all: all-test')
makef.write('\n')

# #
makef.write('\n')
makef.write('base_dir = .\n')
makef.write('src_dirs = $(STREAMER_HOME)/src \\\n')
makef.write('\t$(STREAMER_HOME)/src/amr_mesh \\\n')
makef.write('\t$(STREAMER_HOME)/src/cdr_solver \\\n')
makef.write('\t$(STREAMER_HOME)/src/elliptic \\\n')
makef.write('\t$(STREAMER_HOME)/src/geometry \\\n')
makef.write('\t$(STREAMER_HOME)/src/global \\\n')
makef.write('\t$(STREAMER_HOME)/src/plasma_solver \\\n')
makef.write('\t$(STREAMER_HOME)/src/poisson_solver \\\n')
makef.write('\t$(STREAMER_HOME)/src/rte_solver \\\n')
makef.write('\t$(STREAMER_HOME)/src/sigma_solver \\\n')
makef.write('\t$(STREAMER_HOME)/geometries_prebuilt \\\n')
makef.write('\t$(STREAMER_HOME)/cell_taggers \\\n')
makef.write('\t$(STREAMER_HOME)/time_steppers \\\n')
makef.write('\t$(STREAMER_HOME)/plasma_models/' + args.plasma_kinetics + '\n')
makef.write('\n')

# Define rules to build everything

makef.write('\n')
makef.write('include $(CHOMBO_HOME)/mk/Make.example\n')
makef.write('\n')

makef.close()


# Write an options file. This should be a separate routine
options_filename = app_dir + "/template.inputs"
optf = open(options_filename, 'w')

# Write plasma kinetics options
optf.write("# ====================================================================================================\n")
optf.write('# POTENTIAL CURVE\n')
optf.write("# ====================================================================================================\n")
optf.write(args.app_name + ".potential = 1\n")
optf.write('\n')
options_files = [args.streamer_home + "/src/geometry/physical_domain.options", \
                 args.streamer_home + "/geometries_prebuilt/" + args.geometry + ".options", \
                 args.streamer_home + "/src/amr_mesh/amr_mesh.options", \
                 args.streamer_home + "/src/plasma_solver/plasma_engine.options", \
                 args.streamer_home + "/src/plasma_solver/time_stepper.options", \
                 args.streamer_home + "/src/poisson_solver/poisson_solver.options", \
                 args.streamer_home + "/src/poisson_solver/poisson_multifluid_gmg.options", \
                 args.streamer_home + "/src/plasma_solver/cdr_layout.options", \
                 args.streamer_home + "/src/plasma_solver/rte_layout.options", \
                 args.streamer_home + "/src/rte_solver/eddington_sp1.options", \
                 args.streamer_home + "/plasma_models/" + args.plasma_kinetics + "/" + args.plasma_kinetics + ".options"]
for opt in options_files:
    if os.path.exists(opt):
        with open(opt) as f:
            lines = f.readlines()
            optf.writelines(lines)
            optf.write('\n\n')
    else:
        print 'Could not find options file ' + opt
optf.close()
