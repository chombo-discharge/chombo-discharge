# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_muscl.H\"/\#include <CD_CdrMuscl.H>/g' $i
    sed -i 's/\#include <cdr_muscl.H>/\#include <CD_CdrMuscl.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_muscl/CdrMuscl/g' $i
done

# Move files
mv src/CdrSolver/cdr_muscl.H       src/CdrSolver/CD_CdrMuscl.H
mv src/CdrSolver/cdr_muscl.cpp     src/CdrSolver/CD_CdrMuscl.cpp
mv src/CdrSolver/cdr_muscl.options src/CdrSolver/CD_CdrMuscl.options

# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_muscl/CdrMuscl/g' $i
    sed -i 's/compute_slopes/computeSlopes/g' $i
    sed -i 's/compute_bndry_outflow/computeBoundaryOutflow/g' $i
done
