# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"TiledMeshRefine.H\"/\#include \"CD_TiledMeshRefine.H\"/g' $i
    
    sed -i 's/\#include <TiledMeshRefine.H>/\#include <CD_TiledMeshRefine.H>/g' $i
done

# # Move files
mv src/AmrMesh/TiledMeshRefine.H    src/AmrMesh/CD_TiledMeshRefine.H
mv src/AmrMesh/TiledMeshRefine.cpp  src/AmrMesh/CD_TiledMeshRefine.cpp
