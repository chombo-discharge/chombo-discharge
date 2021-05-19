#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"hollow_cylinder_if.H\"/\#include <CD_HollowCylinderIF.H>/g' $i
    sed -i 's/\#include <hollow_cylinder_if.H>/\#include <CD_HollowCylinderIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/hollow_cylinder_if/HollowCylinderIF/g' $i
done

# Move files

mv Source/Geometry/hollow_cylinder_if.H   Source/Geometry/CD_HollowCylinderIF.H
mv Source/Geometry/hollow_cylinder_if.cpp Source/Geometry/CD_HollowCylinderIF.cpp
