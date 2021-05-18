# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_tga.H\"/\#include <CD_CdrTGA.H>/g' $i
    sed -i 's/\#include <cdr_tga.H>/\#include <CD_CdrTGA.H>/g' $i
done


# Move files
mv src/CdrSolver/cdr_tga.H    src/CdrSolver/CD_CdrTGA.H
mv src/CdrSolver/cdr_tga.cpp  src/CdrSolver/CD_CdrTGA.cpp

# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_tga/CdrTGA/g' $i
done
