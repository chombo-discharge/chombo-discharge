#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"BoundingVolumes.H\"/\#include <CD_BoundingVolumes.H>/g' $i
    sed -i 's/\#include <BoundingVolumes.H>/\#include <CD_BoundingVolumes.H>/g' $i
    sed -i 's/\#include \"BoundingVolumesI.H\"/\#include <CD_BoundingVolumesImplem.H>/g' $i
    sed -i 's/\#include <BoundingVolumesI.H>/\#include <CD_BoundingVolumesImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rte_physics_species/RtPhysicsSpecies/g' $i
    sed -i 's/rte_stepper/RtPhysicsStepper/g' $i
done

# Move files

mv Source/Geometry/BoundingVolumes.H  Source/Geometry/CD_BoundingVolumes.H
mv Source/Geometry/BoundingVolumesI.H Source/Geometry/CD_BoundingVolumesImplem.H
