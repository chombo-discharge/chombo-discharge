#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rte_solver.H\"/\#include <CD_RtSolver.H>/g' $i
    sed -i 's/\#include <rte_solver.H>/\#include <CD_RtSolver.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rte_solver/RtSolver/g' $i
    sed -i 's/set_rte_species/setRtSpecies/g' $i
    sed -i 's/set_stationary/setStationary/g' $i
    sed -i 's/is_stationary/isStationary/g' $i
    sed -i 's/get_time/getTime/g' $i
    sed -i 's/get_kappa/getKappa/g' $i
    sed -i 's/get_kappa_eb/getKappaEb/g' $i
done

# Move files
mv Source/RadiativeTransfer/rte_solver.H       Source/RadiativeTransfer/CD_RtSolver.H
mv Source/RadiativeTransfer/rte_solver.cpp     Source/RadiativeTransfer/CD_RtSolver.cpp
