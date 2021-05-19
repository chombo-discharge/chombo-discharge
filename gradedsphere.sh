#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"graded_perlin_sphere_if.H\"/\#include <CD_GradedPerlinSphereSdf.H>/g' $i
    sed -i 's/\#include <graded_perlin_sphere_if.H>/\#include <CD_GradedPerlinSphereSdf.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/graded_perlin_sphere_if/GradedPerlinSphereSdf/g' $i
done

# Move files

mv Source/Geometry/graded_perlin_sphere_if.H   Source/Geometry/CD_GradedPerlinSphereSdf.H
mv Source/Geometry/graded_perlin_sphere_if.cpp Source/Geometry/CD_GradedPerlinSphereSdf.cpp
