# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"EBFastCoarToCoarRedist.H\"/\#include \"CD_EbFastCoarToCoarRedist.H\"/g' $i
    sed -i 's/\#include \"EBFastCoarToCoarRedistImplem.H\"/\#include \"CD_EbFastCoarToCoarRedistImplem.H\"/g' $i
    
    sed -i 's/\#include <EBFastCoarToCoarRedist.H>/\#include <CD_EbFastCoarToCoarRedist.H>/g' $i
    sed -i 's/\#include <EBFastCoarToCoarRedistImplem.H>/\#include <CD_EbFastCoarToCoarRedistImplem.H>/g' $i
done

# Find and replace all instances of EBFastCoarToCoarRedist. This is now EbCentroidInterpolation
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/EBFastCoarToCoarRedist/EbFastCoarToCoarRedist/g' $i
done

# Move files
mv src/AmrMesh/EBFastCoarToCoarRedist.H    src/AmrMesh/CD_EbFastCoarToCoarRedist.H
mv src/AmrMesh/EBFastCoarToCoarRedist.cpp  src/AmrMesh/CD_EbFastCoarToCoarRedist.cpp
    
