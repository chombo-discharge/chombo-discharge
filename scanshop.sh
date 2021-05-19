#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ScanShop.H\"/\#include <CD_ScanShop.H>/g' $i
    sed -i 's/\#include <ScanShop.H>/\#include <CD_ScanShop.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rte_physics_species/RtPhysicsSpecies/g' $i
    sed -i 's/rte_stepper/RtPhysicsStepper/g' $i
done

# Move files

mv Source/geometry Source/Geometry
mv Source/Geometry/ScanShop.H Source/Geometry/CD_ScanShop.H
mv Source/Geometry/ScanShop.cpp Source/Geometry/CD_ScanShop.cpp

for i in `find . -type f \( -iname \*.py -o -iname *GNUmakefile \)`; do
    sed -i 's/Source\/geometry/Source\/Geometry/g' $i
done
