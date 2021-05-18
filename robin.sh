# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include <robin_coef.H>/\#include <CD_RobinCoefficients.H>/g' $i
    sed -i 's/\#include \"robin_coef.H\"/\#include <CD_RobinCoefficients.H>/g' $i
    sed -i 's/\#include <robin_coefI.H>/\#include <CD_RobinCoefficientsImplem.H>/g' $i
    sed -i 's/\#include \"robin_coefI.H\"/\#include <CD_RobinCoefficientsImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/robin_coef/RobinCoefficients/g' $i
done

# Move files
mv Source/Elliptic/robin_coef.H     Source/Elliptic/CD_RobinCoefficients.H
mv Source/Elliptic/robin_coef.cpp   Source/Elliptic/CD_RobinCoefficients.cpp

