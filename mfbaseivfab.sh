#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"MFBaseIVFAB.H\"/\#include <CD_MFBaseIVFAB.H>/g' $i
    sed -i 's/\#include <MFBaseIVFAB.H>/\#include <CD_MFBaseIVFAB.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/get_ivfab/getIVFAB/g' $i
    sed -i 's/getPhase_ptr/getPhasePtr/g' $i
done

# Move files
mv Source/Utilities/MFBaseIVFAB.H    Source/Utilities/CD_MFBaseIVFAB.H
