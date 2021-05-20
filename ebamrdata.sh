#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"EBAMRData.H\"/\#include <CD_EBAMRData.H>/g' $i
    sed -i 's/\#include <EBAMRData.H>/\#include <CD_EBAMRData.H>/g' $i
    sed -i 's/\#include \"EBAMRDataI.H\"/\#include <CD_EBAMRDataImplem.H>/g' $i
    sed -i 's/\#include <EBAMRDataI.H>/\#include <CD_EBAMRDataImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/get_data/getData/g' $i
done


# Move files
mv Source/Utilities/EBAMRData.H     Source/Utilities/CD_EBAMRData.H
mv Source/Utilities/EBAMRData.cpp   Source/Utilities/CD_EBAMRData.cpp
mv Source/Utilities/EBAMRDataI.H    Source/Utilities/CD_EBAMRDataImplem.H
