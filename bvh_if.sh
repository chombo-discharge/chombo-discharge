#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"bvh_if.H\"/\#include <CD_BvhSdf.H>/g' $i
    sed -i 's/\#include <bvh_if.H>/\#include <CD_BvhSdf.H>/g' $i
    sed -i 's/\#include \"bvh_ifI.H\"/\#include <CD_BvhSdfImplem.H>/g' $i
    sed -i 's/\#include <bvh_ifI.H>/\#include <CD_BvhSdfImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/bvh_if/BvhSdf/g' $i
done

# Move files

mv Source/Geometry/bvh_if.H   Source/Geometry/CD_BvhSdf.H
mv Source/Geometry/bvh_ifI.H  Source/Geometry/CD_BvhSdfImplem.H
