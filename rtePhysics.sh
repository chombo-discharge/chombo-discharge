#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"rte_physics_species.H\"/\#include <CD_RtPhysicsSpecies.H>/g' $i
    sed -i 's/\#include <rte_physics_species.H>/\#include <CD_RtPhysicsSpecies.H>/g' $i
    sed -i 's/\#include \"rte_stepper.H\"/\#include <CD_RtPhysicsStepper.H>/g' $i
    sed -i 's/\#include <rte_stepper.H>/\#include <CD_RtPhysicsStepper.H>/g' $i
    sed -i 's/\#include \"rte_stepperI.H\"/\#include <CD_RtPhysicsStepperImplem.H>/g' $i
    sed -i 's/\#include <rte_stepperI.H>/\#include <CD_RtPhysicsStepperImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/rte_physics_species/RtPhysicsSpecies/g' $i
    sed -i 's/rte_stepper/RtPhysicsStepper/g' $i
done

# Move files

mv Physics/rte/rte_physics_species.H   Physics/rte/CD_RtPhysicsSpecies.H
mv Physics/rte/rte_physics_species.cpp Physics/rte/CD_RtPhysicsSpecies.cpp
mv Physics/rte/rte_stepper.H           Physics/rte/CD_RtPhysicsStepper.H
mv Physics/rte/rte_stepperI.H          Physics/rte/CD_RtPhysicsStepperImplem.H

mv Physics/rte Physics/RadiativeTransfer


for i in `find . -type f \( -iname \*.py -o -iname *GNUmakefile \)`; do
    sed -i 's/Physics\/rte/Physics\/RadiativeTransfer/g' $i
done
