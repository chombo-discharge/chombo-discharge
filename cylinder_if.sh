#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cylinder_if.H\"/\#include <CD_CylinderSdf.H>/g' $i
    sed -i 's/\#include <cylinder_if.H>/\#include <CD_CylinderSdf.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/cylinder_if/CylinderSdf/g' $i
done

# Move files

mv Source/Geometry/cylinder_if.H   Source/Geometry/CD_CylinderSdf.H
mv Source/Geometry/cylinder_if.cpp Source/Geometry/CD_CylinderSdf.cpp
