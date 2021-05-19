#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"hyperboloid_if.H\"/\#include <CD_HyperboloidIF.H>/g' $i
    sed -i 's/\#include <hyperboloid_if.H>/\#include <CD_HyperboloidIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/hyperboloid_if/HyperboloidIF/g' $i
done

# Move files

mv Source/Geometry/hyperboloid_if.H   Source/Geometry/CD_HyperboloidIF.H
mv Source/Geometry/hyperboloid_if.cpp Source/Geometry/CD_HyperboloidIF.cpp
