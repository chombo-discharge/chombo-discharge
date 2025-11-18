import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.discharge_home + "/Geometries" + "/" + args.geometry + "/CD_" + args.geometry + ".H"
    tsfile  = args.discharge_home + "/Physics/ItoKMC/TimeSteppers" + "/" + args.time_stepper + "/CD_" + args.time_stepper + ".H"
    kinfile = args.discharge_home + "/Physics/ItoKMC/PlasmaModels" + "/" + args.physics + "/CD_" + args.physics + ".H"
    tagfile = args.discharge_home + "/Physics/ItoKMC/CellTaggers" +  "/" + args.cell_tagger + "/CD_" + args.cell_tagger + ".H"
    if not os.path.exists(geofile):
        print('Could not find ' + geofile)
    if not os.path.exists(tsfile):
        print('Could not find ' + tsfile)
    if not os.path.exists(kinfile):
        print('Could not find ' + kinfile)
    if not os.path.exists(tagfile) and args.cell_tagger != "none":
        print('Could not find ' + tagfile)
                    
    # Create app directory if it does not exist
    app_dir = args.base_dir + "/" + args.app_name
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
        # Write main file. This should be a separate routine. 
    main_filename = app_dir + "/main.cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include <CD_Driver.H>\n')
    mainf.write('#include <CD_' + args.physics + '.H>\n')
    mainf.write('#include <CD_' + args.geometry + '.H>\n')
    mainf.write('#include <CD_' + args.time_stepper + '.H>\n')
    if not args.cell_tagger == "none":
        mainf.write('#include <CD_' + args.cell_tagger + '.H>\n')
    mainf.write("\n")
    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::ItoKMC;\n")
    mainf.write("\n")    
    mainf.write("using I = " + args.ito_solver + ";\n")
    mainf.write("using C = " + args.cdr_solver + ";\n")
    mainf.write("using R = " + args.rte_solver + ";\n")
    mainf.write("using F = " + args.field_solver + ";\n")
    mainf.write("\n")        
    mainf.write("int main(int argc, char* argv[]){\n")
    mainf.write("  ChomboDischarge::initialize(argc, argv);\n")
    mainf.write("\n")    
    mainf.write("  Random::seed();\n");    
    mainf.write("\n")
    mainf.write("  // Define a voltage curve\n")
    mainf.write("  Real U0 = 1.0;\n")
    mainf.write("  ParmParse pp;\n")
    mainf.write('  pp.get("voltage", U0);\n')
    mainf.write("  auto U = [&](const Real t) -> Real { return U0;};\n")    
    mainf.write("\n")
    mainf.write("  auto compgeom    = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  auto amr         = RefCountedPtr<AmrMesh> (new AmrMesh());\n")
    mainf.write("  auto physics     = RefCountedPtr<ItoKMCPhysics> (new " + args.physics + "());\n")
    mainf.write("  auto timestepper = RefCountedPtr<ItoKMCStepper<I, C, R, F>> (new " + args.time_stepper + "<>(physics));\n")
    if args.cell_tagger != "none":
        mainf.write("  auto tagger      = RefCountedPtr<CellTagger> (new " + args.cell_tagger + "<ItoKMCStepper<I, C, R, F>>(physics, timestepper, amr));\n")
    else:
        mainf.write("  auto tagger      = RefCountedPtr<CellTagger> (nullptr);\n")
    mainf.write("  auto driver      = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger));\n")
    mainf.write("\n")
    mainf.write("  timestepper->setVoltage(U);\n")
    mainf.write("\n")
    mainf.write("  driver->setupAndRun();\n");
    mainf.write("\n")
    mainf.write("  ChomboDischarge::finalize();\n")
    mainf.write("}\n")
        
    mainf.close()
