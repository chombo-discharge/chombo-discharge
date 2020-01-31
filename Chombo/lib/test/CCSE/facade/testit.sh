#!/bin/sh

path_to_exec=.

txt_outfile=/tmp/facade.out
$path_to_exec/testFacade > $txt_outfile

txt_canonical=canonical.out
cmp $txt_canonical $txt_outfile > /dev/null
if test $? -eq 0 ; then
    echo "passed text test"
else
    echo "failed text test.  diff $txt_canonical $txt_outfile"
fi
