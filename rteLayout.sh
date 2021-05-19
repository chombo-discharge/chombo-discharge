#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rte_layout.H\"/\#include <CD_RtLayout.H>/g' $i
    sed -i 's/\#include <rte_layout.H>/\#include <CD_RtLayout.H>/g' $i
    sed -i 's/\#include \"rte_layoutI.H\"/\#include <CD_RtLayoutImplem.H>/g' $i
    sed -i 's/\#include <rte_layoutI.H>/\#include <CD_RtLayoutImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rte_layout/RtLayout/g' $i
    sed -i 's/rte_factory/RtFactory/g' $i
    sed -i 's/add_solver/addSolver/g' $i
    sed -i 's/get_solvers/getSolvers/g' $i
    sed -i 's/get_phase/getPhase/g' $i
    sed -i 's/new_layout/newLayout/g' $i
done

# Move files
mv Source/RadiativeTransfer/rte_layout.H       Source/RadiativeTransfer/CD_RtLayout.H
mv Source/RadiativeTransfer/rte_layoutI.H      Source/RadiativeTransfer/CD_RtLayoutImplem.H
