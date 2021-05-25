import os
import sys

# Write an options file. This should be a separate routine
def write_template(args):
    app_dir = args.discharge_home + "/" + args.base_dir + "/" + args.app_name
    options_filename = app_dir + "/template.inputs"
    optf = open(options_filename, 'w')
    
    # Write plasma kinetics options
    optf.write("# ====================================================================================================\n")
    optf.write('# POTENTIAL CURVE\n')
    optf.write("# ====================================================================================================\n")
    optf.write(args.app_name + ".potential = 1\n")
    optf.write(args.app_name + ".basename  = pout\n")
    optf.write('\n')

    options_files = [args.discharge_home + "/Source/AmrMesh/AmrMesh.options", \
                     args.discharge_home + "/Source/Driver/Driver.options", \
                     args.discharge_home + "//Source/Electrostatics/CD_" + args.field_solver + ".options",\
                     args.discharge_home + "/Source/ItoDiffusion/" + args.ItoSolver + ".options",\
                     args.discharge_home + "/Source/RadiativeTransfer/McPhoto.options",\
                     args.discharge_home + "/Source/Geometry/GeoCoarsener.options", \
                     args.discharge_home + "/Geometries/" + args.geometry + "/" + args.geometry + ".options", \
                     args.discharge_home + "/Physics/ItoPlasma/timesteppers/" + args.TimeStepper + "/" + args.TimeStepper + ".options", \
                     args.discharge_home + "/Physics/ItoPlasma/PlasmaModels/" + args.physics + "/" + args.physics + ".options"]

    if not args.CellTagger == "none":
        options_files.append(args.discharge_home + "/Physics/ItoPlasma/CellTaggers/" + args.CellTagger + "/" + args.CellTagger + ".options")
        
    for opt in options_files:
        if os.path.exists(opt):
            f = open(opt, 'r')
            lines = f.readlines()
            optf.writelines(lines)
            optf.write('\n\n')
            f.close()
        else:
            print 'Could not find options file (this _may_ be normal behavior) ' + opt
    optf.close()
