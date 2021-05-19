#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_BVH.H\"/\#include <CD_DcelBVH.H>/g' $i
    sed -i 's/\#include <dcel_BVH.H>/\#include <CD_DcelBVH.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/dcel_BVH/DcelBVH/g' $i
    sed -i 's/namespace dcel/namespace Dcel/g' $i
    sed -i 's/dcel::/Dcel::/g' $i
done

# Move files

mv Source/Geometry/dcel_BVH.H   Source/Geometry/CD_DcelBVH.H
