#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"particle_container.H\"/\#include <CD_ParticleContainer.H>/g' $i
    sed -i 's/\#include <particle_container.H>/\#include <CD_ParticleContainer.H>/g' $i
    sed -i 's/\#include \"particle_containerI.H\"/\#include <CD_ParticleContainerImplem.H>/g' $i
    sed -i 's/\#include <particle_containerI.H>/\#include <CD_ParticleContainerImplem.H>/g' $i
done



# Move files
mv Source/Particle/particle_container.H     Source/Particle/CD_ParticleContainer.H
mv Source/Particle/particle_containerI.H    Source/Particle/CD_ParticleContainerImplem.H



for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/particle_container/ParticleContainer/g' $i
    sed -i 's/copy_halo_particles/copyMaskParticles/g' $i
    sed -i 's/copy_non_halo_particles/copyNonMaskParticles/g' $i
    sed -i 's/clear_halo_particles/clearMaskParticles/g' $i
    sed -i 's/clear_non_halo_particles/clearNonMaskParticles/g' $i
    sed -i 's/clear_particles/clearParticles/g' $i
    sed -i 's/get_particles/getParticles/g' $i
    sed -i 's/get_halo_particles/getMaskParticles/g' $i
    sed -i 's/get_non_halo_particles/getNonMaskParticles/g' $i
    sed -i 's/get_cache_particles/getCacheParticles/g' $i
    sed -i 's/get_pvr/getPVR/g' $i
    sed -i 's/get_cell_particles_destructive/getCellParticlesDestructive/g' $i
    sed -i 's/get_cell_particles/getCellParticles/g' $i
    sed -i 's/sort_particles_by_cell/sortParticlesByCell/g' $i
    sed -i 's/sort_particles_by_patch/sortParticlesByPatch/g' $i
    sed -i 's/add_particles_destructive/addParticlesDestructive/g' $i
    sed -i 's/add_particles/addParticles/g' $i
    sed -i 's/cache_particles/cacheParticles/g' $i
    sed -i 's/discard_invalid_particles/discardInvalidParticles/g' $i
    sed -i 's/get_num_valid_local/getNumberOfValidParticesLocal/g' $i
    sed -i 's/get_num_valid_global/getNumberOfValidParticesGlobal/g' $i
    sed -i 's/get_num_outcast_local/getNumberOfOutcastParticesLocal/g' $i
    sed -i 's/get_num_outcast_global/getNumberOfOutcastParticesGlobal/g' $i
    sed -i 's/get_num_halo_local/getNumberOfMaskParticesLocal/g' $i
    sed -i 's/get_num_halo_global/getNumberOfMaskParticesGlobal/g' $i
    sed -i 's/get_weight_valid_local/getWeightValidLocal/g' $i
    sed -i 's/get_weight_valid_global/getWeightValidGlobal/g' $i
    sed -i 's/setup_pvr/setupPVR/g' $i
    sed -i 's/setup_particle_data/setupParticleData/g' $i
    sed -i 's/remap_lost_particles/remapLostParticles/g' $i
done
