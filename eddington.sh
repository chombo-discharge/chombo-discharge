#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"eddington_sp1.H\"/\#include <CD_EddingtonSP1.H>/g' $i
    sed -i 's/\#include <eddington_sp1.H>/\#include <CD_EddingtonSP1.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/eddington_sp1/EddingtonSP1/g' $i
    sed -i 's/allocate_WallBc/allocateWallBc/g' $i
    sed -i 's/set_reflection_coefficients/setReflectionCoefficients/g' $i
    sed -i 's/compute_boundary_flux/computeBoundaryFlux/g' $i
    sed -i 's/compute_domain_flux/computeDomainFlux/g' $i
    sed -i 's/compute_density/computeDensity/g' $i
    sed -i 's/m_gmg_relax_type/m_multigridRelaxMethod/g' $i
    sed -i 's/m_gmg_type/m_multigridType/g' $i
    sed -i 's/m_needs_setup/m_needsMultigridSetup/g' $i
    sed -i 's/m_has_mg_stuff/m_hasDeeperMultigridLevels/g' $i
    sed -i 's/m_use_tga/m_useTGA/g' $i
    sed -i 's/m_gmg_verbosity/m_multigridVerbosity/g' $i
    sed -i 's/m_gmg_coarsen/m_numCoarseningsBeforeAggregation/g' $i
    sed -i 's/m_gmg_pre_smooth/m_multigridPreSmooth/g' $i
    sed -i 's/m_gmg_post_smooth/m_multigridPostSmooth/g' $i
    sed -i 's/m_gmg_bot_smooth/m_multigridBottomSmooth/g' $i
    sed -i 's/m_gmg_max_iter/m_multigridMaxIterations/g' $i
    sed -i 's/m_gmg_min_iter/m_multigridMinIterations/g' $i
    sed -i 's/m_bottomsolver/m_bottomSolver/g' $i
    sed -i 's/m_numsmooth/m_numSmoothingsForSimpleSolver/g' $i
    sed -i 's/m_bottom_drop/m_numCellsBottomDrop/g' $i
    sed -i 's/m_gmg_eps/m_multigridTolerance/g' $i
    sed -i 's/m_gmg_hang/m_multigridHang/g' $i
    sed -i 's/m_r1/m_reflectionCoefficientOne/g' $i
    sed -i 's/m_r2/m_reflectionCoefficientTwo/g' $i
    sed -i 's/m_gmg_solver/m_multigridSolver/g' $i
    sed -i 's/m_tgasolver/m_tgaSolver/g' $i
    sed -i 's/m_eulersolver/m_eulerSolver/g' $i
    sed -i 's/m_opfact/m_operatorFactory/g' $i
    sed -i 's/m_robinco/m_robinCoefficients/g' $i
    sed -i 's/m_ebfact/m_robinEbBcFactory/g' $i
    sed -i 's/m_domfact/m_robinDomainBcFactory/g' $i
    sed -i 's/m_wallbc/m_wallBc/g' $i
    sed -i 's/m_simple_solver/m_simpleSolver/g' $i
    sed -i 's/m_aCoefficient/m_aCoef/g' $i
    sed -i 's/setup_gmg/setupMultigrid/g' $i
    sed -i 's/setCoefficientsficients/setMultigridCoefficients/g' $i
    sed -i 's/set_aco_and_bco/setACoefAndBCoef/g' $i
    sed -i 's/set_aco_and_bco_box/setACoefAndBCoefBox/g' $i
    sed -i 's/define_mg_levels/defineDeeperMultigridLevels/g' $i
    sed -i 's/setup_operator_factory/setupOperatorFactory/g' $i
    sed -i 's/setup_multigrid/setupMultigridSolver/g' $i
    sed -i 's/set_neumann_WallBc/setNeumannWallBc/g' $i
    sed -i 's/set_robin_WallBc/setRobinWallBc/g' $i
    sed -i 's/setup_tga/setupTGA/g' $i
    sed -i 's/setup_euler/setupEuler/g' $i
    sed -i 's/parse_stationary/parseStationary/g' $i
    sed -i 's/parse_gmg_settings/parseMultigridSettings/g' $i
    sed -i 's/parse_reflection/parseReflection/g' $i
done

# Move files
mv Source/rte_solver                              Source/RadiativeTransfer
mv Source/RadiativeTransfer/eddington_sp1.H       Source/RadiativeTransfer/CD_EddingtonSP1.H
mv Source/RadiativeTransfer/eddington_sp1.cpp     Source/RadiativeTransfer/CD_EddingtonSP1.cpp
mv Source/RadiativeTransfer/eddington_sp1.options Source/RadiativeTransfer/CD_EddingtonSP1.options

# Update makefiles to use 
for i in `find . -type f \( -iname \*.py -o -iname *GNUmakefile \)`; do
    sed -i 's/Source\/rte_solver/Source\/RadiativeTransfer/g' $i
done
