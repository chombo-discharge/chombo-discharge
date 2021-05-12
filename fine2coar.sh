# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"EBFastFineToCoarRedist.H\"/\#include \"CD_EbFastFineToCoarRedist.H\"/g' $i
    sed -i 's/\#include \"EBFastFineToCoarRedistImplem.H\"/\#include \"CD_EbFastFineToCoarRedistImplem.H\"/g' $i
    
    sed -i 's/\#include <EBFastFineToCoarRedist.H>/\#include <CD_EbFastFineToCoarRedist.H>/g' $i
    sed -i 's/\#include <EBFastFineToCoarRedistImplem.H>/\#include <CD_EbFastFineToCoarRedistImplem.H>/g' $i
done

# Find and replace all instances of EBFastFineToCoarRedist. This is now EbCentroidInterpolation
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/EBFastFineToCoarRedist/EbFastFineToCoarRedist/g' $i
done

# Move files
mv src/AmrMesh/EBFastFineToCoarRedist.H    src/AmrMesh/CD_EbFastFineToCoarRedist.H
mv src/AmrMesh/EBFastFineToCoarRedist.cpp  src/AmrMesh/CD_EbFastFineToCoarRedist.cpp
    
