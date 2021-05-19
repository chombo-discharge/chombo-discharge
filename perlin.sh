#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"polygon_rod_if.H\"/\#include <CD_PolygonRodIF.H>/g' $i
    sed -i 's/\#include <polygon_rod_if.H>/\#include <CD_PolygonRodIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/polygon_rod_if/PolygonRodIF/g' $i
done

# Move files

mv Source/Geometry/polygon_rod_if.H   Source/Geometry/CD_PolygonRodIF.H
mv Source/Geometry/polygon_rod_if.cpp Source/Geometry/CD_PolygonRodIF.cpp
