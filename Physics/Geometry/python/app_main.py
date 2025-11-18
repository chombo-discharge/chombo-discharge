import os
import sys

def write_template(args):
    # Make sure that every class can be found where they should
    geofile = args.discharge_home + "/Geometries" + "/" + args.geometry + "/" + args.geometry + ".H"
                    
    # Create app directory if it does not exist
    app_dir = args.base_dir + "/" + args.app_name
    if not os.path.exists(app_dir):
        os.makedirs(app_dir)
                        
        # Write main file. This should be a separate routine. 
    main_filename = app_dir + "/main.cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include <CD_Driver.H>\n')
    mainf.write('#include "CD_' + args.geometry + '.H"\n')
    mainf.write('#include <CD_GeometryStepper.H>\n')
    mainf.write("\n")
    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::Geometry;\n\n")
    mainf.write("int main(int argc, char* argv[]){\n")
    mainf.write("  ChomboDischarge::initialize(argc, argv);\n")        
    mainf.write("\n")
    mainf.write("  auto compgeom    = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  auto amr         = RefCountedPtr<AmrMesh> (new AmrMesh());\n")
    mainf.write("  auto tagger      = RefCountedPtr<CellTagger> (nullptr);\n")
    mainf.write("  auto timestepper = RefCountedPtr<GeometryStepper> (new GeometryStepper());\n")
    mainf.write("  auto engine      = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger));\n")
    mainf.write("\n")    
    mainf.write("  engine->setupAndRun();\n");
    mainf.write("\n")
    mainf.write("  ChomboDischarge::finalize();\n")
    mainf.write("}\n")    
        
    mainf.close()
