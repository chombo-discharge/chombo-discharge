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
    main_filename = app_dir + "/program.cpp"
    mainf = open(main_filename, "w")
    mainf.write('#include <CD_Driver.H>\n')
    mainf.write('#include <CD_' + args.geometry + '.H>\n')
    mainf.write('#include <CD_StreamerInceptionStepper.H>\n')
    mainf.write('#include <CD_StreamerInceptionTagger.H>\n')    
    mainf.write('#include <ParmParse.H>\n')
    mainf.write("\n")

    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::StreamerInception;\n\n")
    mainf.write("int main(int argc, char* argv[]){\n")

    mainf.write("\n")
    
    mainf.write("#ifdef CH_MPI\n")
    mainf.write("  MPI_Init(&argc, &argv);\n")
    mainf.write("#endif\n")
    
    mainf.write("\n")
    
    mainf.write("  // Read the input file into the ParmParse table\n")
    mainf.write("  const std::string input_file = argv[1];\n")
    mainf.write("  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());\n")
    
    mainf.write("\n")
    
    mainf.write("  // Set geometry and AMR \n")
    mainf.write("  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());\n")

    mainf.write("\n")
    
    mainf.write("  // Define transport data\n")
    mainf.write("  auto alpha         = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto eta           = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto alphaEff      = [&](const Real& E, const RealVect& x) -> Real { return alpha(E,x) - eta(E,x); };\n")
    mainf.write("  auto bgRate        = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto detachRate    = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto fieldEmission = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto ionMobility   = [&](const Real& E) -> Real { return 0.0; };\n")
    mainf.write("  auto ionDiffusion  = [&](const Real& E) -> Real { return 0.0; };\n")    
    mainf.write("  auto ionDensity    = [&](const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto voltageCurve  = [&](const Real& time) -> Real { return 1.0;};\n")
    
    mainf.write("\n")

    mainf.write("  // Set up time stepper \n")
    mainf.write("  auto timestepper = RefCountedPtr<StreamerInceptionStepper<>> (new StreamerInceptionStepper<>());\n")
    mainf.write("  auto celltagger  = RefCountedPtr<StreamerInceptionTagger> (new StreamerInceptionTagger(amr, timestepper->getElectricField(), alphaEff));\n");
    
    mainf.write("\n")

    mainf.write("  // Set transport data\n")
    mainf.write("  timestepper->setAlpha(alpha);\n")
    mainf.write("  timestepper->setEta(eta);\n")
    mainf.write("  timestepper->setBackgroundRate(bgRate);\n")
    mainf.write("  timestepper->setDetachmentRate(detachRate);\n")
    mainf.write("  timestepper->setFieldEmission(fieldEmission);\n")
    mainf.write("  timestepper->setIonMobility(ionMobility);\n")
    mainf.write("  timestepper->setIonDiffusion(ionDiffusion);\n")    
    mainf.write("  timestepper->setIonDensity(ionDensity);\n")
    mainf.write("  timestepper->setVoltageCurve(voltageCurve);\n")

    mainf.write("\n");
    
    mainf.write("  // Set up the Driver and run it\n")
    mainf.write("  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, celltagger));\n")
    mainf.write("  engine->setupAndRun(input_file);\n");
    
    mainf.write("\n")

    mainf.write("#ifdef CH_MPI\n")
    mainf.write("  CH_TIMER_REPORT();\n")
    mainf.write("  MPI_Finalize();\n")
    mainf.write("#endif\n")
    mainf.write("}\n")
        
    mainf.close()
