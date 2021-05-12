# # Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"load_balance.H\"/\#include \"CD_LoadBalancing.H\"/g' $i
#     sed -i 's/\#include \"load_balanceI.H\"/\#include \"CD_LoadBalancingImplem.H\"/g' $i
    
#     sed -i 's/\#include <load_balance.H>/\#include <CD_LoadBalancing.H>/g' $i
#     sed -i 's/\#include <load_balanceI.H>/\#include <CD_LoadBalancingImplem.H>/g' $i
# done

# # Find and replace all instances of load_balance. This is now EbCentroidInterpolation
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/load_balance/LoadBalancing/g' $i
# done

# # Replace function names
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/make_balance/makeBalance/g' $i
#     sed -i 's/level_by_level/balanceLevelByLevel/g' $i
#     sed -i 's/round_robin/roundRobin/g' $i
#     sed -i 's/gather_boxes/gatherBoxes/g' $i
#     sed -i 's/gather_loads/gatherLoads/g' $i
#     sed -i 's/gather_boxes_loads/gatherBoxesAndLoads/g' $i
#     sed -i 's/std_sort/standardSort/g' $i
#     sed -i 's/morton_sort/mortonSort/g' $i
# done

# # Move files
# mv src/AmrMesh/load_balance.H    src/AmrMesh/CD_LoadBalancing.H
# mv src/AmrMesh/load_balance.cpp  src/AmrMesh/CD_LoadBalancing.cpp
# mv src/AmrMesh/load_balanceI.H  src/AmrMesh/CD_LoadBalancingImplem.H

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/gatherBoxes_loads/gatherBoxesAndLoads/g' $i
    sed -i 's/BoxSorting::none/BoxSorting::None/g' $i
    sed -i 's/BoxSorting::std/BoxSorting::Std/g' $i
    sed -i 's/BoxSorting::shuffle/BoxSorting::Shuffle/g' $i
    sed -i 's/BoxSorting::morton/BoxSorting::Morton/g' $i
done
    
