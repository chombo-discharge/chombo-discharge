import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.discharge_home + "/geometries" + "/" + args.geometry + "/" + args.geometry + ".H"
    tsfile  = args.discharge_home + "/physics/ito_plasma/timesteppers" + "/" + args.TimeStepper + "/" + args.TimeStepper + ".H"
    kinfile = args.discharge_home + "/physics/ito_plasma/plasma_models" + "/" + args.physics + "/" + args.physics + ".H"
    tagfile = args.discharge_home + "/physics/ito_plasma/CellTaggers" +  "/" + args.CellTagger + "/" + args.CellTagger + ".H"
    if not os.path.exists(geofile):
        print 'Could not find ' + geofile
    if not os.path.exists(tsfile):
        print 'Could not find ' + tsfile
    if not os.path.exists(kinfile):
        print 'Could not find ' + kinfile
    if not os.path.exists(tagfile) and args.CellTagger != "none":
        print 'Could not find ' + tagfile
                    
    # Create app directory if it does not exist
    app_dir = args.discharge_home + "/" + args.base_dir + "/" + args.app_name
    print(app_dir)
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
        # Write main file. This should be a separate routine. 
    main_filename = app_dir + "/" + args.filename + ".cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include "CD_Driver.H"\n')
    mainf.write('#include "geo_coarsener.H"\n')
    mainf.write('#include "CD_FieldSolverFactory.H"\n')
    mainf.write('#include "CD_' + args.field_solver + '.H"\n')
    mainf.write('#include "ito_layout.H"\n')
    mainf.write('#include "' + args.ito_solver + '.H"\n')
    mainf.write('#include "rte_layout.H"\n')
    mainf.write('#include "mc_photo.H"\n')
    mainf.write('#include "' + args.physics + '.H"\n')
    mainf.write('#include "' + args.geometry + '.H"\n')
    mainf.write('#include "' + args.TimeStepper + '.H"\n')
    if not args.CellTagger == "none":
        mainf.write('#include "' + args.CellTagger + '.H"\n')
    mainf.write('#include "ParmParse.H"\n')
    mainf.write("\n")
    
    mainf.write("// This is the potential curve (constant in this case). Modify it if you want to.\n")
    mainf.write("Real g_potential;\n");
    mainf.write("Real potential_curve(const Real a_time){\n")
    mainf.write("  return g_potential;\n")
    mainf.write("}\n")

    mainf.write("\n")
    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace physics::ito_plasma;\n\n")
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
    mainf.write("  // Get potential from input script \n")
    mainf.write("  std::string basename; \n");
    mainf.write("  {\n")
    mainf.write('     ParmParse pp("' + args.app_name + '");\n')
    mainf.write('     pp.get("potential", g_potential);\n')
    mainf.write('     pp.get("basename",  basename);\n')
    mainf.write('     setPoutBaseName(basename);\n')
    mainf.write("  }\n")

    mainf.write("\n")
    mainf.write("  // Set geometry and AMR \n")
    mainf.write("  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new " + args.geometry + "());\n")
    mainf.write("  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());\n")
    mainf.write("  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());\n")

    mainf.write("\n")
    mainf.write("  // Set up physics \n")
    mainf.write("  RefCountedPtr<ito_plasma_physics> physics      = RefCountedPtr<ito_plasma_physics> (new " + args.physics + "());\n")
    mainf.write("  RefCountedPtr<ito_plasma_stepper> timestepper  = RefCountedPtr<ito_plasma_stepper> (new " + args.TimeStepper + "(physics));\n")
    if args.CellTagger != "none":
        mainf.write("  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (new " + args.CellTagger + "(physics, timestepper, amr, compgeom));\n")
    else:
        mainf.write("  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (NULL);\n")

    mainf.write("\n")

    mainf.write("  // Create solver factories\n")
    mainf.write("  auto poi_fact = new FieldSolverFactory<" + args.field_solver + ">();\n")
    mainf.write("  auto ito_fact = new ito_factory<ito_solver, " + args.ito_solver + ">();\n")
    mainf.write("  auto rte_fact = new rte_factory<mc_photo, mc_photo>();\n")
    mainf.write("\n")
    
    mainf.write("  // Instantiate solvers\n")
    mainf.write("  auto poi = poi_fact->newSolver();\n");
    mainf.write("  auto cdr = ito_fact->new_layout(physics->get_ito_species());\n");
    mainf.write("  auto rte = rte_fact->new_layout(physics->get_rte_species());\n");
    mainf.write("\n")

    mainf.write("  // Send solvers to TimeStepper \n")
    mainf.write("  timestepper->set_poisson(poi);\n");
    mainf.write("  timestepper->set_ito(cdr);\n");
    mainf.write("  timestepper->set_rte(rte);\n");
    mainf.write("\n")

    mainf.write("  // Set potential \n")
    mainf.write("timestepper->set_potential(potential_curve);\n")
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
