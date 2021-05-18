#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"EBGhostCloud.H\"/\#include <CD_EbGhostCloud.H>/g' $i
    sed -i 's/\#include <EBGhostCloud.H>/\#include <CD_EbGhostCloud.H>/g' $i
    sed -i 's/\#include <EBGhostCloudF_F.H>/\#include <CD_EbGhostCloudF_F.H>/g' $i
    sed -i 's/\#include \"EBGhostCloudF_F.H\"/\#include <CD_EbGhostCloudF_F.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/EBGhostCloud/EbGhostCloud/g' $i
done

# Move files
mv Source/Particle/EBGhostCloud.H   Source/Particle/CD_EbGhostCloud.H
mv Source/Particle/EBGhostCloud.cpp Source/Particle/CD_EbGhostCloud.cpp
mv Source/Particle/EBGhostCloudF_F.H  Source/Particle/CD_EbGhostCloudF_F.H
mv Source/Particle/EBGhostCloudF.ChF  Source/Particle/CD_EbGhostCloudF.ChF


