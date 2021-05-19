#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"BVH.H\"/\#include <CD_BVH.H>/g' $i
    sed -i 's/\#include <BVH.H>/\#include <CD_BVH.H>/g' $i
    sed -i 's/\#include \"BVHI.H\"/\#include <CD_BVHImplem.H>/g' $i
    sed -i 's/\#include <BVHI.H>/\#include <CD_BVHImplem.H>/g' $i
done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
#     sed -i 's/BVH/Bvh/g' $i
# done

# Move files
mv Source/Geometry/BVH.H   Source/Geometry/CD_BVH.H
mv Source/Geometry/BVHI.H  Source/Geometry/CD_BVHImplem.H
