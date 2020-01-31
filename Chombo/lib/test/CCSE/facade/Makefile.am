CCSE = ccse
include $(top_srcdir)/lib/mk/autoconf/Automake.rules

if !SMALLBUILD_AM
bin_PROGRAMS = testFacade
endif !SMALLBUILD_AM

bindir = $(libexecdir)/Chombo/tests/$(CCSE)/facade

nodist_testFacade_SOURCES = testFacade.cpp
testFacade_CPPFLAGS = $(CH_KPPFLAGS) $(DASH_I_CHOMBOFACADE)
testFacade_CXXFLAGS =
testFacade_CFLAGS = $(CH_KPPFLAGS) $(DASH_I_CHOMBOFACADE)

AM_LDFLAGS += -lboxlib$(DIM)D -L$(pkglibdir) -lboxtools$(DIM)D -lfacade$(DIM)D

testit_SCRIPTS = testit.sh
testitdir = $(libexecdir)/Chombo/tests/ccse/facade
canonical_DATA = canonical.out
canonicaldir = $(libexecdir)/Chombo/tests/ccse/facade

EXTRA_DIST += $(srcdir)/*.cpp testit.sh canonical.out
include $(srcdir)/MakefileAm.extra
