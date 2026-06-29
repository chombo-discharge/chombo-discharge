# These are the names for the various chombo-discharge components. Because we build various components into libraries,
# it makes sense to give the libraries a name through the makefile system. Here, SOURCE_LIB is the library name for the
# the library compiled from $DISCHARGE_HOME/Source, GEOMETRY_LIB is the library name for the files compiled from
# $DISCHARGE_HOME/Geometries and so on.
SOURCE_LIB     = Source
GEOMETRIES_LIB = Geometries

# Physics libraries
ADVDIFF_LIB            = AdvectionDiffusionPhysics
BROWNIAN_LIB           = BrownianWalkerPhysics
CDRPLASMA_LIB          = CdrPlasmaPhysics
ELECTROSTATICS_LIB     = ElectrostaticPhysics
GEOMETRYONLY_LIB       = GeometryPhysics
ITOKMC_LIB             = ItoKMCPhysics
MESHODE_LIB            = MeshODEPhysics
RADTRANSFER_LIB        = RadiativeTransferPhysics
DISCHARGEINCEPTION_LIB = DischargeInceptionPhysics
TRACERPARTICLE_LIB     = TracerParticlePhysics

# Headers where the chombo-discharge source code is located. This should all folders
# under $DISCHARGE_HOME/Source
SOURCE_DIRS     := $(shell find $(DISCHARGE_HOME)/Source                     -type d -print)
GEOMETRIES_DIRS := $(shell find $(DISCHARGE_HOME)/Geometries                 -type d -print)

SOURCE_INCLUDE          := $(foreach dir, $(SOURCE_DIRS),          $(addprefix -I, $(dir)))
GEOMETRIES_INCLUDE      := $(foreach dir, $(GEOMETRIES_DIRS),      $(addprefix -I, $(dir)))

# Same as for source and geometries dirs/includes, but for the physics modules. 
ADVDIFF_DIRS            := $(shell find $(DISCHARGE_HOME)/Physics/AdvectionDiffusion -type d -print)
BROWNIAN_DIRS           := $(shell find $(DISCHARGE_HOME)/Physics/BrownianWalker     -type d -print)
CDRPLASMA_DIRS          := $(shell find $(DISCHARGE_HOME)/Physics/CdrPlasma          -type d -print)
ELECTROSTATICS_DIRS     := $(shell find $(DISCHARGE_HOME)/Physics/Electrostatics     -type d -print)
GEOMETRYONLY_DIRS       := $(shell find $(DISCHARGE_HOME)/Physics/Geometry           -type d -print)
ITOKMC_DIRS             := $(shell find $(DISCHARGE_HOME)/Physics/ItoKMC             -type d -print)
MESHODE_DIRS            := $(shell find $(DISCHARGE_HOME)/Physics/MeshODE            -type d -print)
RADTRANSFER_DIRS        := $(shell find $(DISCHARGE_HOME)/Physics/RadiativeTransfer  -type d -print)
DISCHARGEINCEPTION_DIRS := $(shell find $(DISCHARGE_HOME)/Physics/DischargeInception -type d -print)
TRACERPARTICLE_DIRS     := $(shell find $(DISCHARGE_HOME)/Physics/TracerParticle     -type d -print)

ADVDIFF_INCLUDE            := $(foreach dir, $(ADVDIFF_DIRS),           $(addprefix -I, $(dir)))
BROWNIAN_INCLUDE           := $(foreach dir, $(BROWNIAN_DIRS),          $(addprefix -I, $(dir)))
CDRPLASMA_INCLUDE          := $(foreach dir, $(CDRPLASMA_DIRS),         $(addprefix -I, $(dir)))
ELECTROSTATICS_INCLUDE     := $(foreach dir, $(ELECTROSTATICS_DIRS),    $(addprefix -I, $(dir)))
GEOMETRYONLY_INCLUDE       := $(foreach dir, $(GEOMETRYONLY_DIRS),      $(addprefix -I, $(dir)))
ITOKMC_INCLUDE             := $(foreach dir, $(ITOKMC_DIRS),            $(addprefix -I, $(dir)))
MESHODE_INCLUDE            := $(foreach dir, $(MESHODE_DIRS),           $(addprefix -I, $(dir)))
RADTRANSFER_INCLUDE        := $(foreach dir, $(RADTRANSFER_DIRS),       $(addprefix -I, $(dir)))
DISCHARGEINCEPTION_INCLUDE := $(foreach dir, $(DISCHARGEINCEPTION_DIRS), $(addprefix -I, $(dir)))
TRACERPARTICLE_INCLUDE     := $(foreach dir, $(TRACERPARTICLE_DIRS),    $(addprefix -I, $(dir)))

# Source and Geometries headers should always be visible. 
XTRACPPFLAGS += $(SOURCE_INCLUDE) 
XTRACPPFLAGS += $(GEOMETRIES_INCLUDE)

# EBGeometry submodule needs to be visible.
XTRACPPFLAGS += -I$(DISCHARGE_HOME)/Submodules/EBGeometry

# Particle payload precision (sets ParticleReal via -DCD_PARTICLE_REAL). Analogous to Chombo's
# PRECISION flag, but it governs ONLY the SoA particle payload columns -- position and weight are
# always double regardless. Set it on the command line (make PARTICLE_PRECISION=FLOAT) or in
# Lib/Local/Make.defs.local; it defaults to DOUBLE. The value is resolved lazily (deferred) because
# this file is included before Make.defs.local, and it is appended only to our own XTRACPPFLAGS --
# so no modification of Chombo's makefile system is required.
PARTICLE_PRECISION ?= DOUBLE
XTRACPPFLAGS += -DCD_PARTICLE_REAL=$(if $(filter FLOAT,$(PARTICLE_PRECISION)),float,double)

# Fold BOTH floating-point precisions into Chombo's configuration string, in the order
# <Chombo PRECISION>.<payload PARTICLE_PRECISION>, so every build artifact (o/<config>/,
# lib<name><config>.a, main<config>.ex) is uniquely named by the precision of the mesh data AND of
# the particle payload -- e.g. ...MPI.DOUBLE.FLOAT.ex. Distinct configs keep the object trees
# separate, so a precision switch rebuilds cleanly and reuses each tree on switch-back. Without this,
# Chombo's dependency tracking -- which keys only on source/header files, never on -D flag VALUES --
# would never rebuild on a PARTICLE_PRECISION change: the switch would be silently ignored, or a
# partial rebuild would mix objects with different sizeof(ParticleReal) into one ODR-violating binary.
# Both halves are always explicit (DOUBLE included), unlike Chombo's native _precision which omits it.
#
# PRECISION must be given a default here, for the SAME reason PARTICLE_PRECISION is: the top-level
# orchestrator (GNUmakefile) includes THIS file and exports XTRACONFIG to its recursive Chombo
# sub-build, but it never reads Chombo's Make.defs.local where PRECISION is normally set. Without a
# default $(PRECISION) is empty there, so the orchestrator builds/names the Chombo libraries as
# ...OPENMPCC..DOUBLE while every app build (which DOES read Make.defs.local) computes
# ...OPENMPCC.DOUBLE.DOUBLE -- a permanent mismatch that makes the prebuilt Chombo libraries
# invisible and forces a second, redundant Chombo rebuild. `?=' yields to a real Make.defs.local or
# command-line setting (e.g. PRECISION=FLOAT), so this only fills in the orchestrator's blank.
PRECISION ?= DOUBLE
#
# The append is resolved lazily, so $(PRECISION)/$(PARTICLE_PRECISION) pick up values from
# Make.defs.local, which is read after this file. It must also be IDEMPOTENT: Chombo exports
# XTRACONFIG (mk/Make.rules), so recursive $(MAKE) sub-builds -- e.g. an app's `dependencies' rule
# building discharge-lib, or the top-level orchestrator building Chombo -- inherit it through the
# environment. A bare `+=' would re-append the token at every recursion level (...DOUBLE.FLOAT.DOUBLE
# .FLOAT) and break library/object lookup. The one-shot sentinel below, exported so nested makes see
# it, makes the append happen exactly once across the whole recursive build.
ifndef CD_CONFIG_PRECISION_DONE
XTRACONFIG += .$(PRECISION).$(PARTICLE_PRECISION)
CD_CONFIG_PRECISION_DONE := 1
endif
export CD_CONFIG_PRECISION_DONE

# Source and Geometries libraries should always be visible.
XTRALIBFLAGS += $(addprefix -l, $(SOURCE_LIB))$(config)
XTRALIBFLAGS += $(addprefix -l, $(GEOMETRIES_LIB))$(config)
XTRALIBFLAGS += -L/$(DISCHARGE_HOME)/Lib

# As a rule we always use EB (embedded boundaries) and MF (multi-fluid) from
# Chombo
USE_EB = TRUE
USE_MF = TRUE

# Chombo libraries needed for building chombo-discharge
LibNames = MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop ParticleTools
