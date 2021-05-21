# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ito_iterator.H\"/\#include <CD_ItoIterator.H>/g' $i
    sed -i 's/\#include <ito_iterator.H>/\#include <CD_ItoIterator.H>/g' $i
    sed -i 's/\#include \"ito_iteratorI.H\"/\#include <CD_ItoIteratorImplem.H>/g' $i
    sed -i 's/\#include <ito_iteratorI.H>/\#include <CD_ItoIteratorImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/ito_iterator/ItoIterator/g' $i
done

# Move files
mv Source/ItoDiffusion/ito_iterator.H        Source/ItoDiffusion/CD_ItoIterator.H
mv Source/ItoDiffusion/ito_iteratorI.H       Source/ItoDiffusion/CD_ItoIteratorImplem.H

