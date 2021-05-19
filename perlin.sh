#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"torus_if.H\"/\#include <CD_TorusSdf.H>/g' $i
    sed -i 's/\#include <torus_if.H>/\#include <CD_TorusSdf.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/torus_if/TorusSdf/g' $i
done

# Move files

mv Source/Geometry/torus_if.H   Source/Geometry/CD_TorusSdf.H
mv Source/Geometry/torus_if.cpp Source/Geometry/CD_TorusSdf.cpp
