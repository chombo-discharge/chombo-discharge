#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"new_sphere_if.H\"/\#include <CD_SphereSdf.H>/g' $i
    sed -i 's/\#include <new_sphere_if.H>/\#include <CD_SphereSdf.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/new_sphere_if/SphereSdf/g' $i
done

# Move files

mv Source/Geometry/new_sphere_if.H   Source/Geometry/CD_SphereSdf.H
mv Source/Geometry/new_sphere_if.cpp Source/Geometry/CD_SphereSdf.cpp
