# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"centroid_interp.H\"/\#include \"CD_CentroidInterpolationStencil.H\"/g' $i
    sed -i 's/\#include <centroid_interp.H>/\#include <CD_CentroidInterpolationStencil.H>/g' $i
done

# Find and replace all instances of centroid_interp. This is now EbCentroidInterpolation
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/centroid_interp/CentroidInterpolationStencil/g' $i
done

# Move files
mv src/AmrMesh/centroid_interp.H        src/AmrMesh/CD_CentroidInterpolationStencil.H
mv src/AmrMesh/centroid_interp.cpp      src/AmrMesh/CD_CentroidInterpolationStencil.cpp
    
