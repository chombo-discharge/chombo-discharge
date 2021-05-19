#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_if.H\"/\#include <CD_DcelSdf.H>/g' $i
    sed -i 's/\#include <dcel_if.H>/\#include <CD_DcelSdf.H>/g' $i
    sed -i 's/\#include \"dcel_ifI.H\"/\#include <CD_DcelSdfImplem.H>/g' $i
    sed -i 's/\#include <dcel_ifI.H>/\#include <CD_DcelSdfImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/dcel_if/DcelSdf/g' $i
done

# Move files

mv Source/Geometry/dcel_if.H    Source/Geometry/CD_DcelSdf.H
mv Source/Geometry/dcel_ifI.H   Source/Geometry/CD_DcelSdfImplem.H
