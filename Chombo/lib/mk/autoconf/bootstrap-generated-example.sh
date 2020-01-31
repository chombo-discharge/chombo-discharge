#!/bin/sh
#
# Generates metamakefiles for autoconf build.
#
# To delete use ./zap.
#

bootstrapdir=`dirname $0`
if ! test $bootstrapdir = "."; then
    echo "*****************************************************"
    echo "Error: you must run bootstrap from its own directory."
    echo "*****************************************************"
    exit 1
fi

set -x

for f in `find . -name GNUmakefile`; do mv $f ${f}.bak; done

#
# Process configure.in and Makefile.am's.
#
mkdir -p ./config
libtoolize --force
aclocal                   # Produces aclocal.m4, so autoconf can understand
                          # AM_* statements found in Makefile.am
automake --add-missing --foreign # Produces Makefile.in, and various symlinks.
autoconf                  # Produces configure.
