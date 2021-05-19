#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"profile_plane_if.H\"/\#include <CD_ProfilePlaneIF.H>/g' $i
    sed -i 's/\#include <profile_plane_if.H>/\#include <CD_ProfilePlaneIF.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/profile_plane_if/ProfilePlaneIF/g' $i
done

# Move files

mv Source/Geometry/profile_plane_if.H   Source/Geometry/CD_ProfilePlaneIF.H
mv Source/Geometry/profile_plane_if.cpp Source/Geometry/CD_ProfilePlaneIF.cpp
