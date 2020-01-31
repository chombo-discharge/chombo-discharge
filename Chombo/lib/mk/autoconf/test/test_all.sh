#!/bin/sh

cd `dirname $0`/../../../..
TESTHOME=lib/mk/autoconf/test

for multidim in FALSE TRUE; do
  export multidim

  echo "@@ MULTIDIM=$multidim"
  cat $TESTHOME/../Make.defs | sed "s/\(MULTIDIM=\).*/\1$multidim/" > /tmp/t.t
  mv /tmp/t.t $TESTHOME/../Make.defs

  echo "@@ anag_style_testbuilds.sh"  
  ./$TESTHOME/anag_style_testbuilds.sh
  if test $? -ne 0; then
      echo "@@ Failure: anag_style_testbuilds.sh"
      exit 1
  fi
  
  echo "@@ testbuilds.sh 1"
  ./$TESTHOME/testbuilds.sh 1
  if test $? -ne 0; then
      echo "@@ Failure: testbuilds.sh"
      exit 2
  fi
done
