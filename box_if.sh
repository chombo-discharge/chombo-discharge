#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"box_if.H\"/\#include <CD_BoxSdf.H>/g' $i
    sed -i 's/\#include <box_if.H>/\#include <CD_BoxSdf.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/box_if/BoxSdf/g' $i
done

# Move files

mv Source/Geometry/box_if.H   Source/Geometry/CD_BoxSdf.H
mv Source/Geometry/box_if.cpp Source/Geometry/CD_BoxSdf.cpp
