# #Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"mushroom_if.H\"/\#include <CD_MushroomIF.H>/g' $i
#     sed -i 's/\#include <mushroom_if.H>/\#include <CD_MushroomIF.H>/g' $i
# done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/Source\/global/Source\/Utilities/g' $i
done

# Move files

# mv Source/Geometry/mushroom_if.H   Source/Geometry/CD_MushroomIF.H
# mv Source/Geometry/mushroom_if.cpp Source/Geometry/CD_MushroomIF.cpp
mv Source/global Source/Utilities
