#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"MFInterfaceFAB.H\"/\#include <CD_MFInterfaceFAB.H>/g' $i
    sed -i 's/\#include <MFInterfaceFAB.H>/\#include <CD_MFInterfaceFAB.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/get_ivs/getIVS/g' $i
done

# Move files
mv Source/Utilities/MFInterfaceFAB.H    Source/Utilities/CD_MFInterfaceFAB.H
