#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rte_species.H\"/\#include <CD_RtSpecies.H>/g' $i
    sed -i 's/\#include <rte_species.H>/\#include <CD_RtSpecies.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rte_species/RtSpecies/g' $i
    sed -i 's/constant_kappa/isKappaConstant/g' $i
    sed -i 's/get_scatter/getScatteringCoefficient/g' $i
done

# Move files
mv Source/RadiativeTransfer/rte_species.H       Source/RadiativeTransfer/CD_RtSpecies.H
mv Source/RadiativeTransfer/rte_species.cpp     Source/RadiativeTransfer/CD_RtSpecies.cpp
