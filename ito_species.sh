# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ito_particle.H\"/\#include <CD_ItoParticle.H>/g' $i
    sed -i 's/\#include <ito_particle.H>/\#include <CD_ItoParticle.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/ito_particle/ItoParticle/g' $i
    sed -i 's/set_num_runtime_scalars/setNumRunTimeScalars/g' $i
    sed -i 's/set_num_runtime_vectors/setNumRunTimeVectors/g' $i
    sed -i 's/runtime_scalar/runtimeScalar/g' $i
    sed -i 's/runtime_vector/runtimeVector/g' $i
    sed -i 's/simple_ItoParticle/simple_ito_particle/g' $i
done

# Move files
mv Source/ItoDiffusion/ito_particle.H      Source/ItoDiffusion/CD_ItoParticle.H
mv Source/ItoDiffusion/ito_particle.cpp    Source/ItoDiffusion/CD_ItoParticle.cpp

