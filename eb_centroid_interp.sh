# Find and replace all instances
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/EbCentroidInterpolation/EbCentroidInterpolationStencil/g' $i
done

# Move files
mv src/AmrMesh/CD_EbCentroidInterpolation.H      src/AmrMesh/CD_EbCentroidInterpolationStencil.H
mv src/AmrMesh/CD_EbCentroidInterpolation.cpp    src/AmrMesh/CD_EbCentroidInterpolationStencil.cpp

    
