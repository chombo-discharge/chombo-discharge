# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"nwoebquadcfinterp.H\"/\#include \"CD_NwoEbQuadCfInterp.H\"/g' $i
    sed -i 's/\#include <nwoebquadcfinterp.H>/\#include <CD_NwoEbQuadCfInterp.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/nwoebquadcfinterp/NwoEbQuadCfInterp/g' $i
done

# Move files
mv src/AmrMesh/nwoebquadcfinterp.cpp src/AmrMesh/CD_NwoEbQuadCfInterp.cpp
mv src/AmrMesh/nwoebquadcfinterp.H   src/AmrMesh/CD_NwoEbQuadCfInterp.H

