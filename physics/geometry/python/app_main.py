import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.plasmac_home + "/geometries" + "/" + args.geometry + "/" + args.geometry + ".H"
                    
    # Create app directory if it does not exist
    app_dir = args.plasmac_home + "/" + args.base_dir + "/" + args.app_name
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
        # Write main file. This should be a separate routine. 
    main_filename = app_dir + "/" + args.filename + ".cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include "driver.H"\n')
    mainf.write('#include "' + args.geometry + '.H"\n')
    mainf.write('#include "geometry_stepper.H"\n')
    mainf.write('#include "ParmParse.H"\n')
    mainf.write("\n")

    mainf.write("using namespace physics::geometry;\n\n")
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
    mainf.write("  // Set geometry and AMR \n")
    mainf.write("  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new " + args.geometry + "());\n")
    mainf.write("  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());\n")
    mainf.write("  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());\n")
    mainf.write("  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (NULL);\n")

    mainf.write("\n")
    mainf.write("  // Set up basic geometry stepper = 1 \n")
    mainf.write("  auto timestepper = RefCountedPtr<geometry_stepper> (new geometry_stepper());\n")
    mainf.write("\n")
    
    mainf.write("  // Set up the driver and run it\n")
    mainf.write("  RefCountedPtr<driver> engine = RefCountedPtr<driver> (new driver(compgeom, timestepper, amr, tagger, geocoarsen));\n")
    mainf.write("  engine->setup_and_run();\n");
    mainf.write("\n")

    if args.use_mpi:
        mainf.write("#ifdef CH_MPI\n")
        mainf.write("  CH_TIMER_REPORT();\n")
        mainf.write("  MPI_Finalize();\n")
        mainf.write("#endif\n")
        mainf.write("}\n")
        
    mainf.close()
