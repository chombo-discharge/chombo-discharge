import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.discharge_home + "/Geometries" + "/" + args.geometry + "/CD_" + args.geometry + ".H"
    tsfile  = args.discharge_home + "/Physics/ItoPlasma/TimeSteppers" + "/" + args.time_stepper + "/CD_" + args.time_stepper + ".H"
    kinfile = args.discharge_home + "/Physics/ItoPlasma/PlasmaModels" + "/" + args.physics + "/CD_" + args.physics + ".H"
    tagfile = args.discharge_home + "/Physics/ItoPlasma/CellTaggers" +  "/" + args.cell_tagger + "/CD_" + args.cell_tagger + ".H"
    if not os.path.exists(geofile):
        print('Could not find ' + geofile)
    if not os.path.exists(tsfile):
        print('Could not find ' + tsfile)
    if not os.path.exists(kinfile):
        print('Could not find ' + kinfile)
    if not os.path.exists(tagfile) and args.cell_tagger != "none":
        print('Could not find ' + tagfile)
                    
    # Create app directory if it does not exist
    app_dir = args.discharge_home + "/" + args.base_dir + "/" + args.app_name
    print(app_dir)
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
        # Write main file. This should be a separate routine. 
    main_filename = app_dir + "/program.cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include <CD_Driver.H>\n')
    mainf.write('#include <CD_GeoCoarsener.H>\n')
    mainf.write('#include <CD_FieldSolverFactory.H>\n')
    mainf.write('#include <CD_' + args.field_solver + '.H>\n')
    mainf.write('#include <CD_ItoLayout.H>\n')
    mainf.write('#include <CD_' + args.ito_solver + '.H>\n')
    mainf.write('#include <CD_RtLayout.H>\n')
    mainf.write('#include <CD_McPhoto.H>\n')
    mainf.write('#include <CD_' + args.physics + '.H>\n')
    mainf.write('#include <CD_' + args.geometry + '.H>\n')
    mainf.write('#include <CD_' + args.time_stepper + '.H>\n')
    if not args.cell_tagger == "none":
        mainf.write('#include <CD_' + args.cell_tagger + '.H>\n')
    mainf.write('#include "ParmParse.H"\n')
    mainf.write("\n")
    
    mainf.write("// This is the potential curve (constant in this case). Modify it if you want to.\n")
    mainf.write("Real g_potential;\n");
    mainf.write("Real potential_curve(const Real a_time){\n")
    mainf.write("  return g_potential;\n")
    mainf.write("}\n")

    mainf.write("\n")
    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::ItoPlasma;\n\n")
    mainf.write("int main(int argc, char* argv[]){\n")

    mainf.write("\n")
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
    mainf.write("  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());\n")
    mainf.write("  RefCountedPtr<GeoCoarsener> geocoarsen        = RefCountedPtr<GeoCoarsener> (new GeoCoarsener());\n")

    mainf.write("\n")
    mainf.write("  // Set up physics \n")
    mainf.write("  RefCountedPtr<ItoPlasmaPhysics> physics      = RefCountedPtr<ItoPlasmaPhysics> (new " + args.physics + "());\n")
    mainf.write("  RefCountedPtr<ItoPlasmaStepper> timestepper  = RefCountedPtr<ItoPlasmaStepper> (new " + args.time_stepper + "(physics));\n")
    if args.cell_tagger != "none":
        mainf.write("  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (new " + args.cell_tagger + "(physics, timestepper, amr, compgeom));\n")
    else:
        mainf.write("  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (NULL);\n")

    mainf.write("\n")

    mainf.write("  // Create solver factories\n")
    mainf.write("  auto poi_fact = new FieldSolverFactory<" + args.field_solver + ">();\n")
    mainf.write("  auto ito_fact = new ItoFactory<ItoSolver, " + args.ito_solver + ">();\n")
    mainf.write("  auto rte_fact = new RtFactory<McPhoto, McPhoto>();\n")
    mainf.write("\n")
    
    mainf.write("  // Instantiate solvers\n")
    mainf.write("  auto poi = poi_fact->newSolver();\n");
    mainf.write("  auto cdr = ito_fact->newLayout(physics->getItoSpecies());\n");
    mainf.write("  auto rte = rte_fact->newLayout(physics->getRtSpecies());\n");
    mainf.write("\n")

    mainf.write("  // Send solvers to TimeStepper \n")
    mainf.write("  timestepper->setFieldSolver(poi);\n");
    mainf.write("  timestepper->setIto(cdr);\n");
    mainf.write("  timestepper->setRadiativeTransferSolvers(rte);\n");
    mainf.write("\n")

    mainf.write("  // Set potential \n")
    mainf.write("timestepper->setVoltage(potential_curve);\n")
    mainf.write("\n")
    
    mainf.write("  // Set up the Driver and run it\n")
    mainf.write("  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger, geocoarsen));\n")
    mainf.write("  engine->setupAndRun(input_file);\n");
    mainf.write("\n")

    mainf.write("#ifdef CH_MPI\n")
    mainf.write("  CH_TIMER_REPORT();\n")
    mainf.write("  MPI_Finalize();\n")
    mainf.write("#endif\n")
    mainf.write("}\n")
        
    mainf.close()
