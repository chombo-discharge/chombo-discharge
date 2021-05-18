# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include <wall_bc.H>/\#include <CD_WallBc.H>/g' $i
    sed -i 's/\#include \"wall_bc.H\"/\#include <CD_WallBc.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/wall_bc/WallBc/g' $i
done

# Move files
mv Source/Elliptic/wall_bc.H     Source/Elliptic/CD_WallBc.H
mv Source/Elliptic/wall_bc.cpp   Source/Elliptic/CD_WallBc.cpp

