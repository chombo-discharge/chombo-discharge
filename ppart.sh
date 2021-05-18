#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"point_particle.H\"/\#include <CD_PointParticle.H>/g' $i
    sed -i 's/\#include <point_particle.H>/\#include <CD_PointParticle.H>/g' $i
done



# Move files
mv Source/Particle/point_particle.H     Source/Particle/CD_PointParticle.H
mv Source/Particle/point_particle.cpp   Source/Particle/CD_PointParticle.cpp



for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/point_particle/PointParticle/g' $i
done
