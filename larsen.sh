#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"larsen_coefs.H\"/\#include <CD_LarsenCoefficients.H>/g' $i
    sed -i 's/\#include <larsen_coefs.H>/\#include <CD_LarsenCoefficients.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/larsen_coefs/LarsenCoefficients/g' $i
done

# Move files
mv Source/RadiativeTransfer/larsen_coefs.H       Source/RadiativeTransfer/CD_LarsenCoefficients.H
mv Source/RadiativeTransfer/larsen_coefs.cpp     Source/RadiativeTransfer/CD_LarsenCoefficients.cpp
