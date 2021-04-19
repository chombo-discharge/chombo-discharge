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
    makef.write('src_dirs = $(DISCHARGE_HOME)/src \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/amr_mesh \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/cell_tagger \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/driver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/elliptic \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/geometry \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/rte_solver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/sigma_solver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/global \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/poisson_solver \\\n')
    makef.write('\t$(DISCHARGE_HOME)/src/particle \\\n')
    makef.write('\t$(DISCHARGE_HOME)/geometries/' + args.geometry + '\\\n')
    makef.write('\t$(DISCHARGE_HOME)/physics/poisson \\\n')
    makef.write('\n')

    makef.write('\n')
    makef.write('include $(CHOMBO_HOME)/mk/Make.example\n')
    makef.write('\n')

    makef.close()
