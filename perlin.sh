#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"perlin_slab_if.H\"/\#include <CD_PerlinSlabSdf.H>/g' $i
    sed -i 's/\#include <perlin_slab_if.H>/\#include <CD_PerlinSlabSdf.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/perlin_slab_if/PerlinSlabSdf/g' $i
done

# Move files

mv Source/Geometry/perlin_slab_if.H   Source/Geometry/CD_PerlinSlabSdf.H
mv Source/Geometry/perlin_slab_if.cpp Source/Geometry/CD_PerlinSlabSdf.cpp
