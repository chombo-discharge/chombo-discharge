#!/bin/sh

#
# Run this on a fresh checkout of Chombo.  It'll build in-source, from a "make dist",
# and out-of-source against the dist.
#

#
# Usage: $0 [1|0]
# Set the optional arg to 1 if you only want the outofsource build.
#
if test $# -eq 1 ; then
  outofsource_only=1
else
  outofsource_only=0
fi

cd `dirname $0`/../../../..

VERSION=`grep AM_INIT_AUTOMAKE configure.pre | cut -d, -f2 | awk '{print $1}'`

#
# In-source: don't bother, as we test it all the time.
#
./bootstrap
ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi

install_prefix=myinstall

if test $multidim = TRUE ; then
  withccse="--without-ccse"
else
  withccse="--with-ccse"
fi
./configure --prefix=`pwd`/$install_prefix --enable-explicit $withccse
ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi

if test $outofsource_only -eq 0; then
  make -j4 install
  ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi
  
  if test -d ./${install_prefix}/libexec/Chombo/tests; then
    echo "@@@ Running in-source tests..."
    ( cd ./${install_prefix}/libexec/Chombo/tests
      outfile=/tmp/runtests.insource.out.$$
      time ./runtests.sh 2>&1 | tee $outfile
      ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi
    
      echo "@@@ number of test failures:"
      grep "finished with status" $outfile | grep -v '0$' | wc -l
    )
  fi
fi

#
# Dist
#
make dist
ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi

tar xvfz Chombo-${VERSION}.tar.gz

if test $outofsource_only -eq 0; then
  ( cd Chombo-${VERSION}
    ./configure DEBUG=FALSE DIM=3 --prefix=`pwd`/${install_prefix} $withccse
    ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi
    make -j4 install
    ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi
  
    if test -d ./${install_prefix}/libexec/Chombo/tests; then
      echo "@@@ Running tests from dist build..."
      ( cd ./${install_prefix}/libexec/Chombo/tests
        outfile=/tmp/runtests.dist-insource.out.$$
        time ./runtests.sh 2>&1 | tee $outfile
        ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi
    
        echo "@@@ number of test failures:"
        grep "finished with status" $outfile | grep -v '0$' | wc -l
      )
    fi
  )
fi
  
#
# Out of source
#
( cd Chombo-${VERSION}
  if test -f lib/Makefile; then
    make distclean
  fi
  ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi
  mkdir -p outofsource
)
( cd Chombo-${VERSION}/outofsource
  ../configure --enable-explicit --prefix=`pwd`/$install_prefix $withccse
  ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi
  make -j4 install
  ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi

  if test -d ./${install_prefix}/libexec/Chombo/tests; then
    echo "@@@ Running tests from out-of-source build from dist..."
    ( cd ./${install_prefix}/libexec/Chombo/tests
      outfile=/tmp/runtests.dist-outofsource.out.$$
      time ./runtests.sh 2>&1 | tee $outfile
      ret=$?; if [ $ret -ne 0 ]; then exit $ret; fi

      echo "@@@ number of test failures:"
      grep "finished with status" $outfile | grep -v '0$' | wc -l
    )
  fi
)
