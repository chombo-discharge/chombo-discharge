# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

lbase := pyparse

# Architecture.
UNAME := $(shell uname)

# names of Chombo libraries needed by this program, in order of search.
ifdef USE_EB
LibNames = Workshop EBTools BoxTools 
else
LibNames = BoxTools 
endif

# Mechanism to retrieve stuff from the IntarWebs.
RETRIEVE := sh $(CHOMBO_HOME)/util/retrieve.sh

# Python library.
PYVERSION := 2.7
PYDIR := Python-$(PYVERSION)
PYFILE := $(PYDIR).tgz
PYTHON := .$(PYDIR)
PYURL := http://www.python.org/ftp/python/$(PYVERSION)/$(PYFILE)
#LIBPYTHON := $(CHOMBO_HOME)/src/PyParse/$(PYDIR)/libpython$(PYVERSION).a
#PYTHON_CPPFLAGS := -I$(PYDIR) -I$(PYDIR)/Include
PYINSTALLDIR := $(CHOMBO_HOME)/src/PyParse/python
LIBPYTHON := $(PYINSTALLDIR)/lib/libpython$(PYVERSION).a
PYTHON_CPPFLAGS := -I$(PYINSTALLDIR)/include/python$(PYVERSION)
PYCONFIGOPTIONS := --disable-shared --without-threads --prefix=$(PYINSTALLDIR)
ifeq ($(value UNAME), Darwin) 
EXTRAPYLIBS := 
else
EXTRAPYLIBS := /usr/lib/libutil.a /usr/lib/libdl.a
endif

$(PYFILE): 
	@$(RETRIEVE) $(PYURL) $(PYFILE)

$(PYTHON): $(PYFILE) 
	@echo "Expanding Python-2.7.tgz..."
	@tar xzf $(PYFILE) 
	@mkdir -p $(PYDIR)/Modules
	cp $(PYDIR)_Modules_Setup.dist $(PYDIR)/Modules/Setup.dist
	@touch $(PYTHON)

# Build the Python interpreter and library. We have to do gymnastics here 
# to keep Chombo's OPT from propagating into Python's build, since the 
# two systems ascribe different meanings to this variable.
$(LIBPYTHON): $(PYTHON)
	@echo "Configuring and building Python..."
	@cd $(PYDIR) && cp configure configure.old && sed -e 's/OPT=\"\$$OPT/OPT=\"/' -e 's/\$$OPT\"/"/' < configure.old > configure
	@cd $(PYDIR) && ./configure $(PYCONFIGOPTIONS) && make clean && make OPT='-O2' && make install

# Extra include directories
XTRACPPFLAGS := $(PYTHON_CPPFLAGS) -DCHOMBO_HOME="\"$(CHOMBO_HOME)\""

# Make sure the USE_EB environment variable is passed along appropriately.
ifeq ($(USE_EB),TRUE)
XTRACPPFLAGS += -DCH_USE_EB
endif

# Suppress parallel builds in this directory.
.NOTPARALLEL: $(LIBPYTHON)

# The Python library must be built before the objects in this directory.
# In fact, it must be built even before dependencies are generated!
XTRATARGETS = $(LIBPYTHON) 

# Mash the objects from the Python library into this one.
# We also need libutil.a for forkpty() and openpty(), and 
# libdl.a for dynamic loading stubs (don't ask).
XTRALIBS = $(LIBPYTHON) $(EXTRAPYLIBS)

include $(CHOMBO_HOME)/mk/Make.lib

