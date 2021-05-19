#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"wedge_if.H\"/\#include <CD_WedgeIF.H>/g' $i
    sed -i 's/\#include <wedge_if.H>/\#include <CD_WedgeIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/wedge_if/WedgeIF/g' $i
done

# Move files

mv Source/Geometry/wedge_if.H   Source/Geometry/CD_WedgeIF.H
mv Source/Geometry/wedge_if.cpp Source/Geometry/CD_WedgeIF.cpp
