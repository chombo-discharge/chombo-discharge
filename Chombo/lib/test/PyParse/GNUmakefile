# -*- Mode: Makefile -*- 

USE_PYTHON=TRUE

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += lib_test_PyParse

# the base name(s) of the application(s) in this directory
ebase = testPyParse

# We need the location of Chombo's lib dir
CHOMBO_HOME := ../..

# names of Chombo libraries needed by this program, in order of search.
ifeq ($(USE_EB), TRUE)
LibNames = PyParse Workshop BoxTools 
else
LibNames = PyParse BoxTools
endif

include $(CHOMBO_HOME)/mk/Make.test

# application-specific variables

# application-specific targets

# Python reference manual -- requires executable to be built.
pydoc: $(ebase)
	@echo "import os, os.path, pydoc\nos.chdir('../doc')\npydoc.writedoc('Workshop')\nos.unlink('../../exec/gendoc.py')\nif os.path.exists('time.table'):\n    os.unlink('time.table')\n" > gendoc.py
	@./$(ebase)$(config).ex gendoc.py

# C++ documentation. Uses Doxygen.
cppdoc: 
	@cd ../doc/cpp ; doxygen doxygen.config

# All documentation.
doc: pydoc cppdoc

