# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ebcoarseaverage.H\"/\#include \"CD_EbCoarAve.H\"/g' $i
    sed -i 's/\#include <ebcoarseaverage.H>/\#include <CD_EbCoarAve.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/ebcoarseaverage/EbCoarAve/g' $i
done

# Move files
mv Source/Utilities/ebcoarseaverage.cpp Source/Utilities/CD_EbCoarAve.cpp
mv Source/Utilities/ebcoarseaverage.H   Source/Utilities/CD_EbCoarAve.H

