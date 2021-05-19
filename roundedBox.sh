#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rounded_box_if.H\"/\#include <CD_RoundedBoxIF.H>/g' $i
    sed -i 's/\#include <rounded_box_if.H>/\#include <CD_RoundedBoxIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rounded_box_if/RoundedBoxIF/g' $i
done

# Move files

mv Source/Geometry/rounded_box_if.H   Source/Geometry/CD_RoundedBoxIF.H
mv Source/Geometry/rounded_box_if.cpp Source/Geometry/CD_RoundedBoxIF.cpp
