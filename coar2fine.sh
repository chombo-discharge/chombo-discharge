# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"EBFastCoarToFineRedist.H\"/\#include \"CD_EbFastCoarToFineRedist.H\"/g' $i
    sed -i 's/\#include \"EBFastCoarToFineRedistImplem.H\"/\#include \"CD_EbFastCoarToFineRedistImplem.H\"/g' $i
    
    sed -i 's/\#include <EBFastCoarToFineRedist.H>/\#include <CD_EbFastCoarToFineRedist.H>/g' $i
    sed -i 's/\#include <EBFastCoarToFineRedistImplem.H>/\#include <CD_EbFastCoarToFineRedistImplem.H>/g' $i
done

# Find and replace all instances of EBFastCoarToFineRedist. This is now EbCentroidInterpolation
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/EBFastCoarToFineRedist/EbFastCoarToFineRedist/g' $i
done

# Move files
mv src/AmrMesh/EBFastCoarToFineRedist.H    src/AmrMesh/CD_EbFastCoarToFineRedist.H
mv src/AmrMesh/EBFastCoarToFineRedist.cpp  src/AmrMesh/CD_EbFastCoarToFineRedist.cpp
    
