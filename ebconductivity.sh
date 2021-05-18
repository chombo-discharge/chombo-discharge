# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include <CD_EbConductivityOp.H>/\#include <CD_EbHelmholtzOp.H>/g' $i
    sed -i 's/\#include <CD_EbConductivityOpFactory.H>/\#include <CD_EbHelmholtzOpFactory.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/EbConductivityOp/EbHelmholtzOp/g' $i
done

# Move files
mv Source/Elliptic/CD_EbConductivityOp.H           Source/Elliptic/CD_EbHelmholtzOp.H
mv Source/Elliptic/CD_EbConductivityOp.cpp         Source/Elliptic/CD_EbHelmholtzOp.cpp
mv Source/Elliptic/CD_EbConductivityOpFactory.H    Source/Elliptic/CD_EbHelmholtzOpFactory.H
mv Source/Elliptic/CD_EbConductivityOpFactory.cpp  Source/Elliptic/CD_EbHelmholtzOpFactory.cpp

