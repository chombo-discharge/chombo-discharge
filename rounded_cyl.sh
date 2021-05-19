#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rounded_cylinder_if.H\"/\#include <CD_RoundedCylinderIF.H>/g' $i
    sed -i 's/\#include <rounded_cylinder_if.H>/\#include <CD_RoundedCylinderIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rounded_cylinder_if/RoundedCylinderIF/g' $i
done

# Move files

mv Source/Geometry/rounded_cylinder_if.H   Source/Geometry/CD_RoundedCylinderIF.H
mv Source/Geometry/rounded_cylinder_if.cpp Source/Geometry/CD_RoundedCylinderIF.cpp
