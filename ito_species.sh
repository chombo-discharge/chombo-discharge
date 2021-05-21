# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"bvh.H\"/\#include <CD_ItoMerge.H>/g' $i
    sed -i 's/\#include <bvh.H>/\#include <CD_ItoMerge.H>/g' $i
    sed -i 's/\#include \"bvhI.H\"/\#include <CD_ItoMergeImplem.H>/g' $i
    sed -i 's/\#include <bvhI.H>/\#include <CD_ItoMergeImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/bvh_node/Node/g' $i
    sed -i 's/bvh_tree/Tree/g' $i
    sed -i 's/build_tree/buildTree/g' $i
    sed -i 's/can_split/canSplit/g' $i
    sed -i 's/is_leaf/isLeaf/g' $i
    sed -i 's/set_leaf/setLeaf/g' $i
    sed -i 's/set_parent/setParent/g' $i
    sed -i 's/set_left/setLeft/g' $i
    sed -i 's/set_right/setRight/g' $i
    sed -i 's/set_data/setData/g' $i
    sed -i 's/set_mass/setMass/g' $i
    sed -i 's/get_parent/getParent/g' $i
    sed -i 's/get_left/getLeft/g' $i
    sed -i 's/get_right/getRight/g' $i
    sed -i 's/get_leaves/getLeaves/g' $i

done

# Move files
mv Source/ItoDiffusion/bvh.H        Source/ItoDiffusion/CD_ItoMerge.H
mv Source/ItoDiffusion/bvhI.H       Source/ItoDiffusion/CD_ItoMergeImplem.H

