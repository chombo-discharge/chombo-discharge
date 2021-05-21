# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"simple_ito_particle.H\"/\#include <CD_SimpleItoParticle.H>/g' $i
    sed -i 's/\#include <simple_ito_particle.H>/\#include <CD_SimpleItoParticle.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/simple_ito_particle/SimpleItoParticle/g' $i
done

# Move files
mv Source/ItoDiffusion/simple_ito_particle.H    Source/ItoDiffusion/CD_SimpleItoParticle.H
mv Source/ItoDiffusion/simple_ito_particle.cpp  Source/ItoDiffusion/CD_SimpleItoParticle.cpp

