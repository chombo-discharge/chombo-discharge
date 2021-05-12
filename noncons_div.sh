# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"noncons_div.H\"/\#include \"CD_NonConservativeDivergenceStencil.H\"/g' $i
    sed -i 's/\#include <noncons_div.H>/\#include <CD_NonConservativeDivergenceStencil.H>/g' $i
done

# Find and replace all instances of noncons_div. This is now EbCentroidInterpolation
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/noncons_div/NonConservativeDivergenceStencil/g' $i
done

# Move files
mv src/AmrMesh/noncons_div.H        src/AmrMesh/CD_NonConservativeDivergenceStencil.H
mv src/AmrMesh/noncons_div.cpp      src/AmrMesh/CD_NonConservativeDivergenceStencil.cpp
    
