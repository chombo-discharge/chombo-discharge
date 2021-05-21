# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ito_species.H\"/\#include <CD_ItoSpecies.H>/g' $i
    sed -i 's/\#include <ito_species.H>/\#include <CD_ItoSpecies.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/ito_species/ItoSpecies/g' $i
done

# Move files
mv Source/ito_solver                      Source/ItoDiffusion
mv Source/ItoDiffusion/ito_species.H      Source/ItoDiffusion/CD_ItoSpecies.H
mv Source/ItoDiffusion/ito_species.cpp    Source/ItoDiffusion/CD_ItoSpecies.cpp


# Update makefiles
for i in `find . -type f \( -iname \*GNUmakefile  -o -iname \*.py \)`; do
    sed -i 's/Source\/ito_solver/Source\/ItoDiffusion/g' $i
done
