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
    mainf.write('#include <CD_' + args.geometry + '.H>\n')
    mainf.write('#include <CD_DischargeInceptionStepper.H>\n')
    mainf.write('#include <CD_DischargeInceptionTagger.H>\n')    
    mainf.write("\n")
    mainf.write("using namespace ChomboDischarge;\n")
    mainf.write("using namespace Physics::DischargeInception;\n\n")
    mainf.write("int main(int argc, char* argv[]){\n")
    mainf.write("  ChomboDischarge::initialize(argc, argv);\n")        
    mainf.write("\n")    
    mainf.write("  auto alpha         = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto eta           = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto alphaEff      = [&](const Real& E, const RealVect& x) -> Real { return alpha(E,x) - eta(E,x); };\n")
    mainf.write("  auto bgRate        = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto detachRate    = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto fieldEmission = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto secondCoeff   = [&](const Real& E, const RealVect& x) -> Real { return 0.0; };\n")    
    mainf.write("  auto ionMobility   = [&](const Real& E) -> Real { return 0.0; };\n")
    mainf.write("  auto ionDiffusion  = [&](const Real& E) -> Real { return 0.0; };\n")    
    mainf.write("  auto ionDensity    = [&](const RealVect& x) -> Real { return 0.0; };\n")
    mainf.write("  auto voltageCurve  = [&](const Real& time) -> Real { return 1.0;};\n")
    mainf.write("\n")
    mainf.write("  auto compgeom    = RefCountedPtr<ComputationalGeometry> (new " + args.geometry + "());\n")
    mainf.write("  auto amr         = RefCountedPtr<AmrMesh> (new AmrMesh());\n")
    mainf.write("  auto timestepper = RefCountedPtr<DischargeInceptionStepper<>> (new DischargeInceptionStepper<>());\n")
    mainf.write("  auto celltagger  = RefCountedPtr<DischargeInceptionTagger> (new DischargeInceptionTagger(amr, timestepper->getElectricField(), alphaEff));\n");
    mainf.write("  auto driver      = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, celltagger));\n")
    mainf.write("\n")    
    mainf.write("  timestepper->setAlpha(alpha);\n")
    mainf.write("  timestepper->setEta(eta);\n")
    mainf.write("  timestepper->setBackgroundRate(bgRate);\n")
    mainf.write("  timestepper->setDetachmentRate(detachRate);\n")
    mainf.write("  timestepper->setFieldEmission(fieldEmission);\n")
    mainf.write("  timestepper->setSecondaryEmission(secondCoeff);\n")    
    mainf.write("  timestepper->setIonMobility(ionMobility);\n")
    mainf.write("  timestepper->setIonDiffusion(ionDiffusion);\n")    
    mainf.write("  timestepper->setIonDensity(ionDensity);\n")
    mainf.write("  timestepper->setVoltageCurve(voltageCurve);\n")
    mainf.write("\n");
    mainf.write("  driver->setupAndRun();\n");
    mainf.write("\n")
    mainf.write("  ChomboDischarge::finalize();\n")
    mainf.write("}\n")
    mainf.close()
