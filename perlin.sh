#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rod_if.H\"/\#include <CD_RodIF.H>/g' $i
    sed -i 's/\#include <rod_if.H>/\#include <CD_RodIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rod_if/RodIF/g' $i
done

# Move files

mv Source/Geometry/rod_if.H   Source/Geometry/CD_RodIF.H
mv Source/Geometry/rod_if.cpp Source/Geometry/CD_RodIF.cpp
