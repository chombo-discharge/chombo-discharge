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

    options_files = [args.discharge_home + "/src/AmrMesh/AmrMesh.options", \
                     args.discharge_home + "/src/Driver/Driver.options", \
                     args.discharge_home + "/src/FieldSolver/CD_" + args.field_solver + ".options",\
                     args.discharge_home + "/src/CdrSolver/" + args.CdrSolver + ".options",\
                     args.discharge_home + "/src/rte_solver/" + args.rte_solver + ".options",\
                     args.discharge_home + "/src/geometry/geo_coarsener.options", \
                     args.discharge_home + "/geometries/" + args.geometry + "/" + args.geometry + ".options", \
                     args.discharge_home + "/physics/cdr_plasma/timesteppers/" + args.TimeStepper + "/" + args.TimeStepper + ".options", \
                     args.discharge_home + "/physics/cdr_plasma/plasma_models/" + args.physics + "/" + args.physics + ".options"]

    if not args.CellTagger == "none":
        options_files.append(args.discharge_home + "/physics/cdr_plasma/CellTaggers/" + args.CellTagger + "/" + args.CellTagger + ".options")
        
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
