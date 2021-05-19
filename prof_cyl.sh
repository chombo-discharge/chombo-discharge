#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"profile_cylinder_if.H\"/\#include <CD_ProfileCylinderIF.H>/g' $i
    sed -i 's/\#include <profile_cylinder_if.H>/\#include <CD_ProfileCylinderIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/profile_cylinder_if/ProfileCylinderIF/g' $i
done

# Move files

mv Source/Geometry/profile_cylinder_if.H   Source/Geometry/CD_ProfileCylinderIF.H
mv Source/Geometry/profile_cylinder_if.cpp Source/Geometry/CD_ProfileCylinderIF.cpp
