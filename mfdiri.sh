# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include <mfdirichletconductivityebbc.H>/\#include <CD_MfElectrostaticDirichletEbBc.H>/g' $i
    sed -i 's/\#include \"mfdirichletconductivityebbc.H\"/\#include <CD_MfElectrostaticDirichletEbBc.H>/g' $i
    sed -i 's/\#include <mfdirichletconductivityebbcI.H>/\#include <CD_MfElectrostaticDirichletEbBcImplem.H>/g' $i
    sed -i 's/\#include \"mfdirichletconductivityebbcI.H\"/\#include <CD_MfElectrostaticDirichletEbBcImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/mfdirichletconductivityebbc/MfElectrostaticDirichletEbBc/g' $i
    sed -i 's/setJump_object/setJumpObject/g' $i
done

# Move files
mv Source/Elliptic/mfdirichletconductivityebbc.H     Source/Elliptic/CD_MfElectrostaticDirichletEbBc.H
mv Source/Elliptic/mfdirichletconductivityebbc.cpp   Source/Elliptic/CD_MfElectrostaticDirichletEbBc.cpp

