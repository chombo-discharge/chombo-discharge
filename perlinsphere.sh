#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"perlin_sphere_if.H\"/\#include <CD_PerlinSphereSdf.H>/g' $i
    sed -i 's/\#include <perlin_sphere_if.H>/\#include <CD_PerlinSphereSdf.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/perlin_sphere_if/PerlinSphereSdf/g' $i
done

# Move files

mv Source/Geometry/perlin_sphere_if.H   Source/Geometry/CD_PerlinSphereSdf.H
mv Source/Geometry/perlin_sphere_if.cpp Source/Geometry/CD_PerlinSphereSdf.cpp
