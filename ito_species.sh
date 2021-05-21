# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ito_solver.H\"/\#include <CD_ItoSolver.H>/g' $i
    sed -i 's/\#include <ito_solver.H>/\#include <CD_ItoSolver.H>/g' $i
    sed -i 's/\#include \"ito_solverI.H\"/\#include <CD_ItoSolverImplem.H>/g' $i
    sed -i 's/\#include <ito_solverI.H>/\#include <CD_ItoSolverImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/ito_solver/ItoSolver/g' $i
    sed -i 's/which_container/WhichContainer/g' $i
    sed -i 's/set_mass_to_conductivity/setMassToConductivity/g' $i
    sed -i 's/unset_mass_to_conductivity/unsetMassToConductivity/g' $i
    sed -i 's/set_mass_to_diffusivity/setMassToDiffusivity/g' $i
    sed -i 's/unset_mass_to_diffusivity/unsetMassToDiffusivity/g' $i
    sed -i 's/set_mass_to_energy/setMassToEnergy/g' $i
    sed -i 's/unset_mass_to_energy/unsetMassToEnergy/g' $i
    sed -i 's/deposit_conductivity/depositConductivity/g' $i
    sed -i 's/deposit_diffusivity/depositDiffusivity/g' $i
    sed -i 's/deposit_energy_density/depositEnergyDensity/g' $i
    sed -i 's/deposit_energy_density/depositEnergyDensity/g' $i
    sed -i 's/compute_average_mobility/computeAverageMobility/g' $i
    sed -i 's/compute_average_diffusion/computeAverageDiffusion/g' $i
    sed -i 's/compute_average_energy/computeAverageEnergy/g' $i
    sed -i 's/deposit_particles/depositParticles/g' $i
    sed -i 's/deposit_weights/depositWeights/g' $i
    sed -i 's/remove_covered_particles/removeCoveredParticles/g' $i
    sed -i 's/transfer_covered_particles/transferCoveredParticles/g' $i
    sed -i 's/intersect_particles_if/intersectParticlesIF/g' $i
    sed -i 's/intersect_particles/intersectParticles/g' $i
    sed -i 's/compute_loads/computeLoads/g' $i
    sed -i 's/get_num_particles/getNumParticles/g' $i
    sed -i 's/writeCheckpointLevel_particles/writeCheckPointLevelParticles/g' $i
    sed -i 's/writeCheckpointLevel_fluid/writeCheckPointLevelFluid/g' $i
    sed -i 's/get_velo_func/getVelocityFunction/g' $i
    sed -i 's/get_diffco_func/getDiffusionFunction/g' $i
    sed -i 's/get_scratch/getScratch/g' $i
    sed -i 's/get_mobility_func/getMobilityFunction/g' $i
    sed -i 's/setDiffusionCoefficient_func/setDiffusionFunction/g' $i
    sed -i 's/setVelocity_func/setVelocityFunction/g' $i
    sed -i 's/interpolate_mobilities_mu/interpolateMobilitiesMu/g' $i
    sed -i 's/interpolate_mobilities_vel/interpolateMobilitiesVel/g' $i
    sed -i 's/interpolate_velocities/interpolateVelocities/g' $i
    sed -i 's/interpolate_mobilities/interpolateMobilities/g' $i
    sed -i 's/update_mobilities/updateMobilities/g' $i
    sed -i 's/interpolate_diffusion/interpolateDiffusion/g' $i
    sed -i 's/update_diffusion/updateDiffusion/g' $i
    sed -i 's/make_superparticles/makeSuperparticles/g' $i
    sed -i 's/random_gaussian/randomGaussian/g' $i
    sed -i 's/compute_min_dt/computeMinDt/g' $i
    sed -i 's/compute_min_drift_dt/computeMinDriftDt/g' $i
    sed -i 's/compute_advective_dt/computeAdvectiveDt/g' $i
    sed -i 's/compute_drift_dt/computeDriftDt/g' $i
    sed -i 's/compute_min_diffusion_dt/computeMinDiffusionDt/g' $i
    sed -i 's/compute_diffusive_dt/computeDiffusiveDt/g' $i
    sed -i 's/which_checkpoint/WhichCheckpoint/g' $i
    sed -i 's/mobility_interp/WhichMobilityInterpolation/g' $i
    sed -i 's/parse_superparticles/parseSuperParticles/g' $i
    sed -i 's/parse_diffusion_hop/parseDiffusionHop/g' $i
    sed -i 's/parse_redistribution/parseRedistribution/g' $i
    sed -i 's/parse_checkpointing/parseCheckpointing/g' $i
    sed -i 's/bvh_merge/mergeBVH/g' $i
    sed -i 's/restart_particles/restartParticles/g' $i
    sed -i 's/remove_covered_particles_if/removeCoveredParticlesIF/g' $i
    sed -i 's/remove_covered_particles_discrete/removeCoveredParticlesDiscrete/g' $i
    sed -i 's/remove_covered_particles_voxels/removeCoveredParticlesVoxels/g' $i
    sed -i 's/transfer_covered_particles_if/transferCoveredParticlesIF/g' $i

    sed -i 's/random_position/randomPosition/g' $i
done

# Move files
mv Source/ItoDiffusion/ito_solver.H        Source/ItoDiffusion/CD_ItoSolver.H
mv Source/ItoDiffusion/ito_solverI.H       Source/ItoDiffusion/CD_ItoSolverImplem.H
mv Source/ItoDiffusion/ito_solver.cpp      Source/ItoDiffusion/CD_ItoSolver.cpp
mv Source/ItoDiffusion/ito_solver.options  Source/ItoDiffusion/CD_ItoSolver.options

