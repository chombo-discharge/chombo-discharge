#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"hyperboloid_one_if.H\"/\#include <CD_HyperboloidOneIF.H>/g' $i
    sed -i 's/\#include <hyperboloid_one_if.H>/\#include <CD_HyperboloidOneIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/hyperboloid_one_if/HyperboloidOneIF/g' $i
done

# Move files

mv Source/Geometry/hyperboloid_one_if.H   Source/Geometry/CD_HyperboloidOneIF.H
mv Source/Geometry/hyperboloid_one_if.cpp Source/Geometry/CD_HyperboloidOneIF.cpp
