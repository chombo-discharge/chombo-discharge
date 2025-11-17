import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.discharge_home + "/Geometries" + "/" + args.geometry + "/" + args.geometry + ".H"
                    
    # Create app directory if it does not exist
    app_dir = args.base_dir + "/" + args.app_name
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
    # Write the main file.
    main_filename = app_dir + "/main.cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include <CD_Driver.H>\n')
    mainf.write('#include <CD_' + args.ItoSolver + '.H>\n')
    mainf.write('#include <CD_' + args.geometry + '.H>\n')
    mainf.write('#include <CD_' + args.stepper + '.H>\n')
    mainf.write('#include <CD_BrownianWalkerTagger.H>\n')
    mainf.write("\n")
    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::BrownianWalker;\n\n")
    mainf.write("int main(int argc, char* argv[]){\n")
    mainf.write("  ChomboDischarge::initialize(argc, argv);\n")
    mainf.write("\n")    
    mainf.write("  auto compgeom    = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  auto amr         = RefCountedPtr<AmrMesh> (new AmrMesh());\n")
    mainf.write("  auto solver      = RefCountedPtr<ItoSolver> (new " + args.ItoSolver + "());\n")
    mainf.write("  auto timestepper = RefCountedPtr<TimeStepper> (new " + args.stepper + "(solver));\n")
    mainf.write("  auto tagger      = RefCountedPtr<CellTagger> (new BrownianWalkerTagger(solver, amr));\n")
    mainf.write("  auto driver      = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger));\n")
    mainf.write("\n")    
    mainf.write("  driver->setupAndRun();\n");
    mainf.write("\n")
    mainf.write("  ChomboDischarge::finalize();\n")
    mainf.write("}\n")    
        
    mainf.close()
