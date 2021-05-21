# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ito_layout.H\"/\#include <CD_ItoLayout.H>/g' $i
    sed -i 's/\#include <ito_layout.H>/\#include <CD_ItoLayout.H>/g' $i
    sed -i 's/\#include \"ito_layoutI.H\"/\#include <CD_ItoLayoutImplem.H>/g' $i
    sed -i 's/\#include <ito_layoutI.H>/\#include <CD_ItoLayoutImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/ito_layout/ItoLayout/g' $i
    sed -i 's/ito_factory/itoFactory/g' $i
    sed -i 's/get_densities/getDensities/g' $i
done

# Move files
mv Source/ItoDiffusion/ito_layout.H        Source/ItoDiffusion/CD_ItoLayout.H
mv Source/ItoDiffusion/ito_layoutI.H       Source/ItoDiffusion/CD_ItoLayoutImplem.H

