import os
import sys

# Write an options file. This should be a separate routine
def write_template(args):
    app_dir = args.base_dir + "/" + args.app_name
    options_filename = app_dir + "/template.inputs"
    optf = open(options_filename, 'w')
    
    # Write plasma kinetics options
    optf.write("# ====================================================================================================\n")
    optf.write('# POTENTIAL CURVE\n')
    optf.write("# ====================================================================================================\n")
    optf.write(args.app_name + ".potential = 1\n")
    optf.write(args.app_name + ".basename  = pout\n")
    optf.write('\n')
    options_files = [args.streamer_home + "/src/amr_mesh/amr_mesh.options", \
                     args.streamer_home + "/src/plasma_solver/driver.options", \
                     args.streamer_home + "/src/poisson_solver/" + args.poisson_solver + ".options",\
                     args.streamer_home + "/src/cdr_solver/" + args.cdr_solver + ".options",\
                     args.streamer_home + "/src/rte_solver/" + args.rte_solver + ".options",\
                     args.streamer_home + "/src/plasma_solver/geo_coarsener.options", \
                     args.streamer_home + "/geometries_prebuilt/" + args.geometry + "/" + args.geometry + ".options", \
                     args.streamer_home + "/time_steppers/" + args.time_stepper + "/" + args.time_stepper + ".options", \
                     args.streamer_home + "/plasma_models/" + args.plasma_kinetics + "/" + args.plasma_kinetics + ".options"]

    if not args.cell_tagger == "none":
        options_files.append(args.streamer_home + "/cell_taggers/" + args.cell_tagger + "/" + args.cell_tagger + ".options")
        
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
