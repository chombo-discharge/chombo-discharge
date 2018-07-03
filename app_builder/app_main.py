import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    app_dir = args.base_dir + "/" + args.app_name
    geofile = args.streamer_home + "/geometries_prebuilt" + "/" + args.geometry + "/" + args.geometry + ".H"
    tsfile  = args.streamer_home + "/time_steppers" + "/" + args.time_stepper + "/" + args.time_stepper + ".H"
    kinfile = args.streamer_home + "/plasma_models" + "/" + args.plasma_kinetics + "/" + args.plasma_kinetics + ".H"
    tagfile = args.streamer_home + "/cell_taggers" +  "/" + args.cell_tagger + "/" + args.cell_tagger + ".H"
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
    mainf.write('#include "geo_coarsener.H"\n')
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
    mainf.write("  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());\n")
    mainf.write("  RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom, compgeom, plaskin, timestepper, amr, tagger, geocoarsen));\n")
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
