if SMALLBUILD_AM
      UNIDIM_DIRS = BoxTools
      DIST_UNIDIM_DIRS = BoxTools
else
      BASIC_DIRS = BoxTools AMRTools AMRElliptic \
                      AMRTimeDependent EllipticDesign ParticleTools
      EB_DIRS    = EBTools EBAMRTools EBAMRElliptic EBAMRTimeDependent \
                      Workshop
      MF_DIRS = MFTools MFElliptic 
      DIST_UNIDIM_DIRS = $(BASIC_DIRS) $(EB_DIRS) $(MF_DIRS)

# If USE_MF_AM==TRUE then we proceed as though USE_EB_AM is also TRUE.
if    USE_MF_AM
          UNIDIM_DIRS = $(BASIC_DIRS) $(EB_DIRS) $(MF_DIRS)
else
if        USE_EB_AM
              UNIDIM_DIRS = $(BASIC_DIRS) $(EB_DIRS)
else
              UNIDIM_DIRS = $(BASIC_DIRS)
endif     !USE_EB_AM
endif !USE_MF_AM

endif !SMALLBUILD_AM

MULTIDIM_DIR = MultiDim

#
# If MULTIDIM_ONLY_AM, we still go into every directory under lib/src, but
# we don't do anything in it; we go straight into lib/src/*/multidim (which
# is where we _didn't_ go under !MULTIDIM_ONLY_AM).
#
if MULTIDIM_ONLY_AM
  SUBDIRS = $(UNIDIM_DIRS) $(MULTIDIM_DIR)
else
  SUBDIRS = $(UNIDIM_DIRS)
endif !MULTIDIM_ONLY_AM

DIST_SUBDIRS = $(DIST_UNIDIM_DIRS) $(MULTIDIM_DIR)

# InfiniteDomain doesn't build -- 3dsine.F wants to include fftw3.f, which
# is nowhere to be found.
