#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"photon.H\"/\#include <CD_Photon.H>/g' $i
    sed -i 's/\#include <photon.H>/\#include <CD_Photon.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/photon/Photon/g' $i
done

# Move files
mv Source/RadiativeTransfer/photon.H       Source/RadiativeTransfer/CD_Photon.H
mv Source/RadiativeTransfer/photon.cpp     Source/RadiativeTransfer/CD_Photon.cpp
