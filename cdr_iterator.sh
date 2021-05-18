# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_iterator.H\"/\#include <CD_CdrIterator.H>/g' $i
    sed -i 's/\#include <cdr_iterator.H>/\#include <CD_CdrIterator.H>/g' $i
    sed -i 's/\#include <cdr_iteratorI.H>/\#include <CD_CdrIterator.H>/g' $i
    sed -i 's/\#include \"cdr_iteratorI.H\"/\#include <CD_CdrIterator.H>/g' $i
done


# Move files
mv src/CdrSolver/cdr_iterator.H    src/CdrSolver/CD_CdrIterator.H
mv src/CdrSolver/cdr_iteratorI.H   src/CdrSolver/CD_CdrIteratorImplem.H

# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_iterator/CdrIterator/g' $i
    sed -i 's/species_iteration/SpeciesIteration/g' $i
done
