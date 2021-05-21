#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"MFLevelGrid.H\"/\#include <CD_MFLevelGrid.H>/g' $i
    sed -i 's/\#include <MFLevelGrid.H>/\#include <CD_MFLevelGrid.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/get_domain/getDomain/g' $i
    sed -i 's/interface_pair/interfacePair/g' $i
done

# Move files
mv Source/Utilities/MFLevelGrid.H      Source/Multifluid/CD_MFLevelGrid.H
mv Source/Utilities/MFLevelGrid.cpp    Source/Multifluid/CD_MFLevelGrid.cpp
