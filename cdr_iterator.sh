# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_iterator.H\"/\#include <CD_CdrIterator.H>/g' $i
    sed -i 's/\#include <cdr_iterator.H>/\#include <CD_CdrIterator.H>/g' $i
done


# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/species_iteration/SpeciesIteration/g' $i
    sed -i 's/num_solvers/getNumberOfSolvers/g' $i
    sed -i 's/get_species/getSpecies/g' $i
    sed -i 's/get_diffco_eb/getEbCenteredDiffusionCoefficient/g' $i
    sed -i 's/getEbDiffusionCoefficients/getEbCenteredDiffusionCoefficient/g' $i
done
