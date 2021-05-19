#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"Vec.H\"/\#include <CD_Vec.H>/g' $i
    sed -i 's/\#include <Vec.H>/\#include <CD_Vec.H>/g' $i
    sed -i 's/\#include \"VecI.H\"/\#include <CD_VecImplem.H>/g' $i
    sed -i 's/\#include <VecI.H>/\#include <CD_VecImplem.H>/g' $i
done


# Move files
mv Source/Geometry/Vec.H    Source/Geometry/CD_Vec.H
mv Source/Geometry/VecI.H   Source/Geometry/CD_VecImplem.H

