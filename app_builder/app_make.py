import os
import sys

def write_template(args):
    app_dir = args.base_dir + "/" + args.app_name
    make_filename = app_dir + "/GNUmakefile"
    makef = open(make_filename, "w")
    makef.write("# Chombo and chombo-streamer directories \n")
    makef.write("CHOMBO_HOME   := " + args.chombo_home + "/lib\n")
    makef.write("STREAMER_HOME := " + args.streamer_home + "\n")
    
    # Make rules
    makef.write("\n")
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
    makef.write("\tAMRTimeDependent BaseTools BoxTools Workshop\n")
    makef.write("\n")
 	
    makef.write("\n")
    makef.write('# Target\n')
    makef.write('all: all-test')
    makef.write('\n')

    makef.write('\n')
    makef.write('base_dir = .\n')
    makef.write('src_dirs = $(STREAMER_HOME)/src \\\n')
    makef.write('\t$(STREAMER_HOME)/src/amr_mesh \\\n')
    makef.write('\t$(STREAMER_HOME)/src/cdr_solver \\\n')
    makef.write('\t$(STREAMER_HOME)/src/elliptic \\\n')
    makef.write('\t$(STREAMER_HOME)/src/geometry \\\n')
    makef.write('\t$(STREAMER_HOME)/src/global \\\n')
    makef.write('\t$(STREAMER_HOME)/src/plasma_solver \\\n')
    makef.write('\t$(STREAMER_HOME)/src/poisson_solver \\\n')
    makef.write('\t$(STREAMER_HOME)/src/rte_solver \\\n')
    makef.write('\t$(STREAMER_HOME)/src/sigma_solver \\\n')
    makef.write('\t$(STREAMER_HOME)/geometries_prebuilt/' + args.geometry + '\\\n')
    makef.write('\t$(STREAMER_HOME)/cell_taggers/'  + args.cell_tagger + '\\\n')
    makef.write('\t$(STREAMER_HOME)/time_steppers/' + args.time_stepper + '\\\n')
    makef.write('\t$(STREAMER_HOME)/plasma_models/' + args.plasma_kinetics + '\n')
    makef.write('\n')

    makef.write('\n')
    makef.write('include $(CHOMBO_HOME)/mk/Make.example\n')
    makef.write('\n')

    makef.close()
