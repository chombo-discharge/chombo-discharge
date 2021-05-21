# #Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"cdr_plasma_physics.H\"/\#include <CD_CdrPlasmaPhysics.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
#     sed -i 's/cdr_plasma_physics/CdrPlasmaPhysics/g' $i
#     sed -i 's/compute_alpha/computeAlpha/g' $i
#     sed -i 's/advance_reaction_network/advanceReactionNetwork/g' $i
#     sed -i 's/compute_cdr_velocities/computeCdrDriftVelocities/g' $i
#     sed -i 's/compute_cdr_diffusion_coefficients/computeCdrDiffusionCoefficients/g' $i
#     sed -i 's/compute_cdr_electrode_fluxes/computeCdrElectrodeFluxes/g' $i
#     sed -i 's/compute_cdr_dielectric_fluxes/computeCdrDielectricFluxes/g' $i
#     sed -i 's/compute_cdr_domain_fluxes/computeCdrDomainFluxes/g' $i
#     sed -i 's/initial_sigma/initialSigma/g' $i
#     sed -i 's/init_eed/initEed/g' $i
#     sed -i 's/get_CdrSpecies/getCdrSpecies/g' $i
#     sed -i 's/get_RtSpecies/getRtSpecies/g' $i
#     sed -i 's/get_num_CdrSpecies/getNumCdrSpecies/g' $i
#     sed -i 's/get_num_RtSpecies/getNumRtSpecies/g' $i
#     sed -i 's/solve_eed/solveEed/g' $i
#     sed -i 's/get_eed_index/getEedIndex/g' $i
#     sed -i 's/set_dt/setDt/g' $i
#     sed -i 's/set_finest_dx/setFinestDx/g' $i
#     sed -i 's/m_num_CdrSpecies/m_numCdrSpecies/g' $i
#     sed -i 's/m_num_RtSpecies/m_numRtSpecies/g' $i
# done

# # Move files
# mv Physics/cdr_plasma Physics/CdrPlasma
# mv Physics/CdrPlasma/cdr_plasma_physics.H   Physics/CdrPlasma/CD_CdrPlasmaPhysics.H

# for i in `find . -type f \( -iname \*.py -o -iname *GNUmakefile \)`; do
#     sed -i 's/Physics\/cdr_plasma/Physics\/CdrPlasma/g' $i
# done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/using namespace physics::cdr_plasma/using namespace Physics::CdrPlasma/g' $i
    sed -i 's/namespace physics/namespace Physics/g' $i
    sed -i 's/namespace cdr_plasma/namespace CdrPlasma/g' $i
done
