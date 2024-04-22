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
