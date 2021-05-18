# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_species.H\"/\#include <CD_CdrSpecies.H>/g' $i
    sed -i 's/\#include <cdr_species.H>/\#include <CD_CdrSpecies.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_species/CdrSpecies/g' $i
done

# Move files
mv src/CdrSolver/cdr_species.H       src/CdrSolver/CD_CdrSpecies.H
mv src/CdrSolver/cdr_species.cpp     src/CdrSolver/CD_CdrSpecies.cpp

# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_species/CdrSpecies/g' $i
    sed -i 's/get_charge/getChargeNumber/g' $i
    sed -i 's/init_with_particles/initializeWithParticles/g' $i
    sed -i 's/init_with_function/initializeWithFunction/g' $i
    sed -i 's/get_deposition/getDeposition/g' $i
    sed -i 's/get_initial_particles/getInitialParticles/g' $i
    sed -i 's/m_charge/m_chargeNumber/g' $i
    sed -i 's/get_unit/getUnit/g' $iw
    sed -i 's/force_output/forceOutput/g' $i
    sed -i 's/m_initial_particles/m_initialParticles/g' $i
done
