#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"point_mass.H\"/\#include <CD_PointMass.H>/g' $i
    sed -i 's/\#include <point_mass.H>/\#include <CD_PointMass.H>/g' $i
done



# Move files
mv Source/Particle/point_mass.H     Source/Particle/CD_PointMass.H
mv Source/Particle/point_mass.cpp   Source/Particle/CD_PointMass.cpp



for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/point_mass/PointMass/g' $i
done
