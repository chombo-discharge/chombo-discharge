import os
import sys

def write_template(args):
    app_dir = args.discharge_home + "/" + args.base_dir + "/" + args.app_name
    make_filename = app_dir + "/GNUmakefile"
    makef = open(make_filename, "w")
    
    # Make rules
    makef.write("include $(CHOMBO_HOME)/mk/Make.defs")
    makef.write("\n")

    makef.write("\n")
    makef.write("USE_EB=TRUE\n")
    makef.write("USE_MF=TRUE\n")
    makef.write("DIM=" + str(args.dim) + "\n")
    makef.write("\n")

    makef.write("\n")
    makef.write("# Base file containing int main()\n")
    makef.write("ebase := " + args.filename)
    makef.write("\n")

    makef.write("\n")
    makef.write("LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \\\n")
    makef.write("\tAMRTimeDependent BaseTools BoxTools Workshop ParticleTools\n")
    makef.write("\n")
 	
    makef.write("\n")
    makef.write('# Target\n')
    makef.write('all: all-test')
    makef.write('\n')

    makef.write('\n')
    makef.write('base_dir = .\n')
    makef.write('src_dirs = $(DISCHARGE_HOME)/Source \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/AmrMesh \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/CdrSolver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/CellTagger \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/Elliptic \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/geometry \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/global \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/Driver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/Particle \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/FieldSolver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/RadiativeTransfer \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Source/SigmaSolver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/Geometries/' + args.geometry + '\\\n')
    makef.write('\t$(DISCHARGE_HOME)/Physics/cdr_plasma \\\n')
    if not args.CellTagger == "none":
        makef.write('\t$(DISCHARGE_HOME)/Physics/cdr_plasma/CellTaggers/'  + args.CellTagger + '\\\n')
    makef.write('\t$(DISCHARGE_HOME)/Physics/cdr_plasma/timesteppers/' + args.TimeStepper + '\\\n')
    makef.write('\t$(DISCHARGE_HOME)/Physics/cdr_plasma/plasma_models/' + args.physics + '\n')
    makef.write('\n')

    makef.write('\n')
    makef.write('include $(CHOMBO_HOME)/mk/Make.example\n')
    makef.write('\n')

    makef.close()
