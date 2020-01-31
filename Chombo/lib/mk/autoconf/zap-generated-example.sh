#!/bin/sh -x
#
# Deletes all autoconf-generated stuff that even
# make distclean misses.
#
rm -f `find . -name Makefile.am`
rm -f `find . -name Makefile.in`
rm -f `find . -name Makefile`
rm -f  aclocal.m4 \
     configure \
     install-sh \
     missing \
     mkinstalldirs \
     config.guess \
     config.sub \
     ltmain.sh \
     depcomp \
     conftest* \
     confdefs* \
     config.log config.status libtool .config/*

rm -rf autom4te.cache
rm -rf `find . -name .deps`
rm -rf `find . -name .libs`
for f in `find . -name GNUmakefile.bak`; do mv $f `dirname $f`/GNUmakefile; done
