import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.discharge_home + "/Geometries" + "/" + args.geometry + "/" + args.geometry + ".H"
                    
    # Create app directory if it does not exist
    app_dir = args.discharge_home + "/" + args.base_dir + "/" + args.app_name
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
        # Write main file. This should be a separate routine. 
    main_filename = app_dir + "/" + args.filename + ".cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include "CD_Driver.H"\n')
    mainf.write('#include "' + args.RtSolver + '.H"\n')
    mainf.write('#include "' + args.geometry + '.H"\n')
    mainf.write('#include <CD_RtPhysicsStepper.H>\n')
    mainf.write('#include "ParmParse.H"\n')
    mainf.write("\n")

    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::RadiativeTransfer;\n\n")
    mainf.write("int main(int argc, char* argv[]){\n")

    mainf.write("\n")
    if args.use_mpi:
        mainf.write("#ifdef CH_MPI\n")
        mainf.write("  MPI_Init(&argc, &argv);\n")
        mainf.write("#endif\n")
    
    mainf.write("\n")
    mainf.write("  // Build class options from input script and command line options\n")
    mainf.write("  const std::string input_file = argv[1];\n")
    mainf.write("  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());")
    mainf.write("\n")

    mainf.write("\n")
    mainf.write("  // Set geometry and AMR \n")
    mainf.write("  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());\n")
    mainf.write("  RefCountedPtr<GeoCoarsener> geocoarsen        = RefCountedPtr<GeoCoarsener> (new GeoCoarsener());\n")
    mainf.write("  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (NULL);\n")

    mainf.write("\n")
    mainf.write("  // Set up basic Poisson, potential = 1 \n")
    mainf.write("  auto timestepper = RefCountedPtr<RtPhysicsStepper<" + args.RtSolver + "> >\n")
    mainf.write("     (new RtPhysicsStepper<" + args.RtSolver + ">());\n")
    mainf.write("\n")



    
    mainf.write("  // Set up the Driver and run it\n")
    mainf.write("  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger, geocoarsen));\n")
    mainf.write("  engine->setupAndRun(input_file);\n");
    mainf.write("\n")

    if args.use_mpi:
        mainf.write("#ifdef CH_MPI\n")
        mainf.write("  CH_TIMER_REPORT();\n")
        mainf.write("  MPI_Finalize();\n")
        mainf.write("#endif\n")
        mainf.write("}\n")
        
    mainf.close()
