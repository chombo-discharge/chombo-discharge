

USE_EB := TRUE
# the base name(s) of the application(s) in this directory
ebase = richardsonConvTest

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib
#CHOMBO_HOME = /Users/dharshid/svn/cell_ave_chombo_branch/lib
USE_NEW_EB = TRUE

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRElliptic EBAMRTools Workshop EBTools AMRTimeDependent AMRElliptic AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../src   

# input file for 'run' target
INPUT = regression.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

XTRALDFLAGS = -llapack -lblas
# application-specific variables

# application-specific targets
