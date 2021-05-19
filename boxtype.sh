#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"BoxType.H\"/\#include <CD_BoxType.H>/g' $i
    sed -i 's/\#include <BoxType.H>/\#include <CD_BoxType.H>/g' $i
done


# Move files
mv Source/Geometry/BoxType.H   Source/Geometry/CD_BoxType.H
mv Source/Geometry/BoxType.cpp Source/Geometry/CD_BoxType.cpp
