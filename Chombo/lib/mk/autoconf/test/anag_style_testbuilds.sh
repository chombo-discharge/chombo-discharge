#!/bin/sh

cd `dirname $0`/../../../..
./bootstrap
./configure
make dist

VERSION=`grep AM_INIT_AUTOMAKE configure.in | cut -d, -f2 | awk '{print $1}'`
tar xvfz Chombo-${VERSION}.tar.gz

cd Chombo-${VERSION}

make -f makefile.anag install

## echo "running ./builds/defaults/install/bin/main..."
## ./builds/defaults/install/bin/main
## if ! test $? -eq 0; then
##     echo "Error: main returned nonzero."
##     exit 1
## fi
## echo "running ./builds/defaults/install/bin/subdirmain..."
## ./builds/defaults/install/bin/subdirmain
## if ! test $? -eq 0; then
##     echo "Error: subdirmain returned nonzero."
##     exit 2
## fi

MAKEDEFS=./lib/mk/autoconf/Make.defs
grep SMALLBUILD $MAKEDEFS | grep TRUE > /dev/null
if test $? -eq 0; then
  subdir=lib/src/smallBoxTools
else
  subdir=lib/src/BoxTools
fi
cd $subdir
make -f makefile.anag DEBUG=FALSE PROFILE=TRUE install

## echo "running ../builds/DEBUG_FALSE.PROFILE_TRUE/install/bin/subdirmain..."
## ../builds/DEBUG_FALSE.PROFILE_TRUE/install/bin/subdirmain
## if ! test $? -eq 0; then
##     echo "Error: subdirmain (DEBUG_FALSE.PROFILE_TRUE) returned nonzero."
##     exit 3
## fi
## 
## if test -f ../builds/DEBUG_FALSE.PROFILE_TRUE/install/bin/main; then
##     echo "Error: main from topsrcdir was built, too."
##     exit 4
## fi
