#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"particle_ops.H\"/\#include <CD_ParticleOps.H>/g' $i
    sed -i 's/\#include <particle_ops.H>/\#include <CD_ParticleOps.H>/g' $i
    sed -i 's/\#include \"particle_opsI.H\"/\#include <CD_ParticleOpsImplem.H>/g' $i
    sed -i 's/\#include <particle_opsI.H>/\#include <CD_ParticleOpsImplem.H>/g' $i
done



# Move files
mv Source/Particle/particle_ops.H     Source/Particle/CD_ParticleOps.H
mv Source/Particle/particle_opsI.H    Source/Particle/CD_ParticleOpsImplem.H



for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/particle_ops/ParticleOps/g' $i
    sed -i 's/domain_intersection/domainIntersection/g' $i
    sed -i 's/eb_intersection_bisect/ebIntersectionBisect/g' $i
    sed -i 's/eb_intersection_raycast/ebIntersectionRaycast/g' $i
done
