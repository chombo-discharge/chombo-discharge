#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"hyperboloid_two_if.H\"/\#include <CD_HyperboloidTwoIF.H>/g' $i
    sed -i 's/\#include <hyperboloid_two_if.H>/\#include <CD_HyperboloidTwoIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/hyperboloid_two_if/HyperboloidTwoIF/g' $i
done

# Move files
mv Source/Geometry/hyperboloid_two_if.H   Source/Geometry/CD_HyperboloidTwoIF.H
mv Source/Geometry/hyperboloid_two_if.cpp Source/Geometry/CD_HyperboloidTwoIF.cpp
