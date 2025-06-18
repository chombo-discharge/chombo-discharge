import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.discharge_home + "/Geometries" + "/" + args.geometry + "/CD_" + args.geometry + ".H"
    tsfile  = args.discharge_home + "/Physics/CdrPlasma/Timesteppers" + "/" + args.time_stepper + "/CD_" + args.time_stepper + ".H"
    kinfile = args.discharge_home + "/Physics/CdrPlasma/PlasmaModels" + "/" + args.physics + "/CD_" + args.physics + ".H"
    tagfile = args.discharge_home + "/Physics/CdrPlasma/CellTaggers" +  "/" + args.cell_tagger + "/CD_" + args.cell_tagger + ".H"
    if not os.path.exists(geofile):
        print('Could not find ' + geofile)
    if not os.path.exists(tsfile):
        print('Could not find ' + tsfile)
    if not os.path.exists(kinfile):
        print('Could not find ' + kinfile)
    if not os.path.exists(tagfile) and args.cell_tagger != "none":
        print('Could not find ' + tagfile)
                    
    # Create app directory if it does not exist
    app_dir = args.discharge_home + '/' + args.base_dir + "/" + args.app_name
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
        # Write main file. This should be a separate routine. 
    main_filename = app_dir + "/program.cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include "CD_Driver.H"\n')
    mainf.write('#include <CD_FieldSolverFactory.H>\n')
    mainf.write('#include <CD_' + args.field_solver + '.H>\n')
    mainf.write('#include <CD_CdrLayoutImplem.H>\n')
    mainf.write('#include <CD_' + args.cdr_solver + '.H>\n')
    mainf.write('#include <CD_RtLayoutImplem.H>\n')
    mainf.write('#include <CD_' + args.rte_solver + '.H>\n')
    mainf.write('#include <CD_' + args.physics + '.H>\n')
    mainf.write('#include <CD_' + args.geometry + '.H>\n')
    mainf.write('#include <CD_' + args.time_stepper + '.H>\n')
    if not args.cell_tagger == "none":
        mainf.write('#include <CD_' + args.cell_tagger + '.H>\n')
    mainf.write('#include <ParmParse.H>\n')
    mainf.write("\n")
    
    mainf.write("// This is the voltage curve (constant in this case). Modify it if you want to.\n")
    mainf.write("Real g_voltage;\n");
    mainf.write("Real voltageCurve(const Real a_time){\n")
    mainf.write("  return g_voltage;\n")
    mainf.write("}\n")

    mainf.write("\n")
    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::CdrPlasma;\n\n")
    mainf.write("int main(int argc, char* argv[]){\n")

    mainf.write("\n")
    mainf.write("#ifdef CH_MPI\n")
    mainf.write("  MPI_Init(&argc, &argv);\n")
    mainf.write("#endif\n")
    
    mainf.write("\n")
    mainf.write("  // Build class options from input script and command line options\n")
    mainf.write("  const std::string input_file = argv[1];\n")
    mainf.write("  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());\n")
    mainf.write("  Random::seed();\n")
    mainf.write("\n")

    mainf.write("\n")
    mainf.write("  // Get voltage from input script \n")
    mainf.write("  std::string basename; \n");
    mainf.write("  {\n")
    mainf.write('    ParmParse pp("' + args.app_name + '");\n')
    mainf.write('    pp.get("voltage", g_voltage);\n')
    mainf.write('    pp.get("basename",  basename);\n')
    mainf.write('    setPoutBaseName(basename);\n')
    mainf.write("  }\n")

    mainf.write("\n")
    mainf.write("  // Set geometry and AMR \n")
    mainf.write("  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());\n")

    mainf.write("\n")
    mainf.write("  // Set up physics \n")
    mainf.write("  RefCountedPtr<CdrPlasmaPhysics> physics      = RefCountedPtr<CdrPlasmaPhysics> (new " + args.physics + "());\n")
    mainf.write("  RefCountedPtr<CdrPlasmaStepper> timestepper  = RefCountedPtr<CdrPlasmaStepper> (new " + args.time_stepper + "(physics));\n")
    if args.cell_tagger != "none":
        mainf.write("  RefCountedPtr<CellTagger> tagger             = RefCountedPtr<CellTagger> (new " + args.cell_tagger + "(physics, timestepper, amr, compgeom));\n")
    else:
        mainf.write("  RefCountedPtr<CellTagger> tagger             = RefCountedPtr<CellTagger> (nullptr);\n")

    mainf.write("\n")

    mainf.write("  // Create solver factories\n")
    mainf.write("  auto poi_fact = new FieldSolverFactory<" + args.field_solver + ">();\n")
    mainf.write("  auto cdr_fact = new CdrFactory<CdrSolver, " + args.cdr_solver + ">();\n")
    mainf.write("  auto rte_fact = new RtFactory<RtSolver, " + args.rte_solver + ">();\n")
    mainf.write("\n")
    
    mainf.write("  // Instantiate solvers\n")
    mainf.write("  auto poi = poi_fact->newSolver();\n");
    mainf.write("  auto cdr = cdr_fact->newLayout(physics->getCdrSpecies());\n");
    mainf.write("  auto rte = rte_fact->newLayout(physics->getRtSpecies());\n");
    mainf.write("\n")

    mainf.write("  // Send solvers to TimeStepper \n")
    mainf.write("  timestepper->setFieldSolver(poi);\n");
    mainf.write("  timestepper->setCdrSolvers(cdr);\n");
    mainf.write("  timestepper->setRadiativeTransferSolvers(rte);\n");
    mainf.write("\n")

    mainf.write("  // Set voltage \n")
    mainf.write("  timestepper->setVoltage(voltageCurve);\n")
    mainf.write("\n")
    
    mainf.write("  // Set up the Driver and run it\n")
    mainf.write("  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger));\n")
    mainf.write("  engine->setupAndRun(input_file);\n");
    mainf.write("\n")

    mainf.write("  // Clean up memory\n")
    mainf.write("  delete poi_fact;\n")
    mainf.write("  delete cdr_fact;\n")
    mainf.write("  delete rte_fact;\n")
    mainf.write("\n")    

    mainf.write("#ifdef CH_MPI\n")
    mainf.write("  CH_TIMER_REPORT();\n")
    mainf.write("  MPI_Finalize();\n")
    mainf.write("#endif\n")
    mainf.write("}\n")
        
    mainf.close()
