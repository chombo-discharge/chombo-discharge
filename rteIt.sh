#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rte_iterator.H\"/\#include <CD_RtIterator.H>/g' $i
    sed -i 's/\#include <rte_iterator.H>/\#include <CD_RtIterator.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rte_iterator/RtIterator/g' $i
done

# Move files
mv Source/RadiativeTransfer/rte_iterator.H       Source/RadiativeTransfer/CD_RtIterator.H

