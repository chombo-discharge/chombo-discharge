# Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"sigma_solver.H\"/\#include \"CD_SigmaSolver.H\"/g' $i
#     sed -i 's/\#include <sigma_solver.H>/\#include <CD_SigmaSolver.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
#     sed -i 's/particle_container/ParticleContainer/g' $i
# done

# Move files
mv Source/particle Source/Particle

for i in `find . -type f \( -iname \*.py -o -iname *GNUmakefile \)`; do
    sed -i 's/\/particle/\/Particle/g' $i
done
