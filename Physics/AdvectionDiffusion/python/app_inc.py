import os
import sys
from shutil import copyfile

# Write an options file. This should be a separate routine
def copy_dependencies(args):
    app_dir = args.base_dir + "/" + args.app_name
    
    kin_home = args.discharge_home + "/Physics/CdrPlasma/PlasmaModels" + "/" + args.physics
    kin_name = args.discharge_home + "/Physics/CdrPlasma/PlasmaModels" + "/" + args.physics + "/" + args.physics

    if os.path.exists(kin_name + ".inc"): # Read file and copy all listed dependencies
        f = open(kin_name + ".inc", 'r')
        incf = f.readline();
        while incf:
            if os.path.exists(kin_home + "/" + incf):
                copyfile(kin_home + "/" + incf, app_dir + "/" + incf)
            incf = f.readline()
        f.close()

        
