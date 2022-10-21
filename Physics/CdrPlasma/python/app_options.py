import os
import sys

# Write an options file. This should be a separate routine
def write_template(args):
    app_dir = args.discharge_home + '/' + args.base_dir + "/" + args.app_name
    options_filename = app_dir + "/template.inputs"
    optf = open(options_filename, 'w')
    
    # Write plasma kinetics options
    optf.write("# ====================================================================================================\n")
    optf.write('# Voltage curve\n')
    optf.write("# ====================================================================================================\n")
    optf.write(args.app_name + ".voltage   = 1\n")
    optf.write(args.app_name + ".basename  = pout\n")
    optf.write('\n')

    options_files = [args.discharge_home + "/Source/AmrMesh/CD_AmrMesh.options", \
                     args.discharge_home + "/Source/Driver/CD_Driver.options", \
                     args.discharge_home + "//Source/Electrostatics/CD_" + args.field_solver + ".options",\
                     args.discharge_home + "//Source/ConvectionDiffusionReaction/CD_" + args.cdr_solver + ".options",\
                     args.discharge_home + "/Source/RadiativeTransfer/CD_" + args.rte_solver + ".options",\
                     args.discharge_home + "/Source/SurfaceODESolver/CD_SurfaceODESolver.options",\
                     args.discharge_home + "/Source/Geometry/CD_GeoCoarsener.options", \
                     args.discharge_home + "/Geometries/" + args.geometry + "/CD_" + args.geometry + ".options", \
                     args.discharge_home + "/Physics/CdrPlasma/Timesteppers/" + args.time_stepper + "/CD_" + args.time_stepper + ".options", \
                     args.discharge_home + "/Physics/CdrPlasma/PlasmaModels/" + args.physics + "/CD_" + args.physics + ".options"]

    if not args.cell_tagger == "none":
        options_files.append(args.discharge_home + "/Physics/CdrPlasma/CellTaggers/" + args.cell_tagger + "/CD_" + args.cell_tagger + ".options")
        
    for opt in options_files:
        if os.path.exists(opt):
            f = open(opt, 'r')
            lines = f.readlines()
            optf.writelines(lines)
            optf.write('\n\n')
            f.close()
        else:
            print('Could not find options file (this _may_ be normal behavior) ' + opt)
    optf.close()
