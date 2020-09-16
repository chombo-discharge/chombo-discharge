#!/bin/csh -f
if ($#argv != 5) then
  echo "Usage: threadrun minthread maxthread timedir exec inputs "
  exit 1
endif

set minProcs = $1
set maxProcs = $2

set timedir = $3
set exec    = $4
set inputs  = $5

set procs = $minProcs
echo "executable name = $exec"
echo "input file name = $inputs"
echo "time directory  = $timedir"

if (! -e $timedir) then
    mkdir $timedir
endif

while ($procs <= $maxProcs)
 echo "number of threads = $procs"
 echo "setenv OMP_NUM_THREADS $procs"
 setenv OMP_NUM_THREADS $procs
 echo "$exec $inputs"
 $exec $inputs
 set pprocs = `echo $procs | awk '{printf("%06d",$1);}'`    
 echo "mv time.table $timedir/amrg.time.$pprocs.thread"
 mv time.table $timedir/amrg.time.$pprocs.thread

#  echo $dir

  set procs = `expr 1 +  $procs`
end
