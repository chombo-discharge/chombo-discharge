# Change inclusion guards. 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"amr_mesh.H\"/\#include \"CD_AmrMesh.H\"/g' $i
#     sed -i 's/\#include <amr_mesh.H>/\#include <CD_AmrMesh.H>/g' $i
    
#     sed -i 's/\#include \"amr_meshImplem.H\"/\#include \<CD_AmrMeshImplem.H\>/g' $i
#     sed -i 's/\#include <amr_meshImplem.H>/\#include <CD_AmrMeshImplem.H>/g' $i
# done

# # Find and replace all instances of amr_mesh everywhere. This is now AmrMesh. 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/amr_mesh/AmrMesh/g' $i
# done

# # Move files
# mv src/amr_mesh/amr_mesh.H        src/amr_mesh/CD_AmrMesh.H
# mv src/amr_mesh/amr_meshImplem.H  src/amr_mesh/CD_AmrMeshImplem.H
# mv src/amr_mesh/amr_mesh.cpp      src/amr_mesh/CD_AmrMesh.cpp
# mv src/amr_mesh/amr_mesh.options  src/amr_mesh/CD_AmrMesh.options

# # Move folder
# mv src/amr_mesh src/AmrMesh

# # Update makefiles
# for i in `find . -name "*GNUmakefile" -type f`; do
#     sed -i 's/amr_mesh/AmrMesh/g' $i
# done


# # Update function signatures
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/register_realm/registerRealm/g' $i
#     sed -i 's/register_operator/registerOperator/g' $i
#     sed -i 's/average_down/averageDown/g' $i
#     sed -i 's/interp_ghost/interpGhost/g' $i
#     sed -i 's/allocate_ptr/allocatePointer/g' $i
#     sed -i 's/register_mask/registerMask/g' $i
#     sed -i 's/query_realm/queryRealm/g' $i
#     sed -i 's/compute_gradient/computeGradient/g' $i
#     sed -i 's/conservative_average/conservativeAverage/g' $i
#     sed -i 's/build_grids/buildGrids/g' $i
#     sed -i 's/interp_ghost_quad/InterpGhostQuad/g' $i
#     sed -i 's/interp_ghost_pwl/InterpGhostPwl/g' $i
#     sed -i 's/interpolate_to_centroids/InterpToCentroids/g' $i
#     sed -i 's/set_mfis/setMultifluidIndexSpace/g' $i
#     sed -i 's/set_baseif/setBaseImplicitFunction/g' $i
#     sed -i 's/parse_domain/parseDomain/g' $i
#     sed -i 's/parse_ghost_interpolation/parseGhostInterpolation/g' $i
#     sed -i 's/build_domains/buildDomains/g' $i
#     sed -i 's/parse_grid_generation/parsegridGeneration/g' $i
#     sed -i 's/parse_verbosity/parseVerbosity/g' $i
#     sed -i 's/parse_coarsest_num_cells/parseCoarsestLevelNumCells/g' $i
#     sed -i 's/parse_max_amr_depth/parseMaxAmrDepth/g' $i
#     sed -i 's/parse_max_sim_depth/parseMaxSimulationDepth/g' $i
#     sed -i 's/parse_ebcf/parseEbCf/g' $i
#     sed -i 's/parse_refinement_ratio/parseRefinementRatios/g' $i
#     sed -i 's/set_refinement_ratios/setRefinementRatios/g' $i
#     sed -i 's/parse_buffer_size/parseBrBufferSize/g' $i
#     sed -i 's/parse_irreg_growth/parseIrregTagGrowth/g' $i
#     sed -i 's/parse_fill_ratio/parseBrFillRatio/g' $i
#     sed -i 's/set_finest_level/setFinestLevel/g' $i
#     sed -i 's/parse_max_box_size/parseMaxBoxSize/g' $i
#     sed -i 's/parse_max_ebis_box_size/parseMaxEbisBoxSize/g' $i
#     sed -i 's/parse_blocking_factor/parseBlockingFactor/g' $i
#     sed -i 's/parse_eb_ghost/parseEbGhostCells/g' $i
#     sed -i 's/parse_num_ghost/parseNumGhostCells/g' $i
#     sed -i 's/parse_redist_rad/parseRedistributionRadius/g' $i
#     sed -i 's/parse_centroid_stencils/parseCentroidStencils/g' $i
#     sed -i 's/parse_eb_stencils/parseEbCentroidStencils/g' $i
#     sed -i 's/set_irreg_sten_type/setIrregularInterpolationStencilType/g' $i
#     sed -i 's/set_irreg_sten_order/setIrregularInterpolationStencilOrder/g' $i
#     sed -i 's/set_irreg_sten_radius/setIrregularInterpolationStencilRadius/g' $i
#     sed -i 's/regrid_amr/regridAmr/g' $i
#     sed -i 's/regrid_realm/regridRealm/g' $i
#     sed -i 's/set_grids/setGrids/g' $i
#     sed -i 's/regrid_operators/regridOperators/g' $i
#     sed -i 's/sanity_check/sanityCheck/g' $i
#     sed -i 's/get_ebcf/getEbCf/g' $i
#     sed -i 's/get_finest_level/getFinestLevel/g' $i
#     sed -i 's/get_irreg_growth/getIrregTagGrowth/g' $i
#     sed -i 's/get_max_amr_depth/getMaxAmrDepth/g' $i
#     sed -i 's/get_max_sim_depth/getMaxSimulationDepth/g' $i
#     sed -i 's/get_refine_all_depth/getRefineAllLevelsDepth/g' $i
#     sed -i 's/get_blocking_factor/getBlockingFactor/g' $i
#     sed -i 's/get_max_box_size/getMaxBoxSize/g' $i
#     sed -i 's/get_num_ghost/getNumberOfGhostCells/g' $i
#     sed -i 's/get_eb_ghost/getNumberOfEbGhostCells/g' $i
#     sed -i 's/get_redist_rad/getRedistributionRadius/g' $i
#     sed -i 's/get_finest_dx/getFinestDx/g' $i
#     sed -i 's/get_prob_lo/getProbLo/g' $i
#     sed -i 's/get_prob_hi/getProbHi/g' $i
#     sed -i 's/get_finest_domain/getFinestDomain/g' $i
#     sed -i 's/get_dx/getDx/g' $i
#     sed -i 's/get_ref_rat/getRefinementRatios/g' $i
#     sed -i 's/get_refinement_ratio/getRefinementRatio/g' $i
#     sed -i 's/get_baseif/getBaseImplicitFunction/g' $i
#     sed -i 's/get_irreg_tags/getIrregularTags/g' $i
#     sed -i 's/get_proxy_grids/getProxyGrids/g' $i
#     sed -i 's/get_grids/getGrids/g' $i
#     sed -i 's/get_domains/getDomains/g' $i
#     sed -i 's/get_mask/getMask/g' $i
#     sed -i 's/get_eblg/getEBLevelGrid/g' $i
#     sed -i 's/get_ebisl/getEBISLayout/g' $i
#     sed -i 's/get_mflg/getMFLevelGrid/g' $i
#     sed -i 's/get_vofit/getVofIterator/g' $i
#     sed -i 's/get_neighbors/getNeighbors/g' $i
#     sed -i 's/get_levelset/getLevelset/g' $i
#     sed -i 's/get_coarave/getCoarseAverage/g' $i
#     sed -i 's/get_ghostcloud/getGhostCloud/g' $i
#     sed -i 's/get_quadcfi/getNWOEBQuadCFInterp/g' $i
#     sed -i 's/get_old_quadcfi/getEBQuadCFInterp/g' $i
#     sed -i 's/get_fillpatch/getFillPatch/g' $i
#     sed -i 's/get_eb_pwl_interp/getPwlInterpolator/g' $i
#     sed -i 's/get_eb_mg_interp/getEBMGInterp/g' $i
#     sed -i 's/get_flux_reg/getFluxRegister/g' $i
#     sed -i 's/get_level_redist/getLevelRedist/g' $i
#     sed -i 's/get_coar_to_fine_redist/getCoarToFineRedist/g' $i
#     sed -i 's/get_coar_to_coar_redist/getCoarToCoarRedist/g' $i
#     sed -i 's/get_fine_to_coar_redist/getFineToCoarRedist/g' $i
#     sed -i 's/get_centroid_interp_stencils/getCentroidInterpolationStencils/g' $i
#     sed -i 's/get_eb_centroid_interp_stencils/getEbCentroidInterpolationStencils/g' $i
#     sed -i 's/get_noncons_div_stencils/getNonConservativeDivergenceStencils/g' $i
#     sed -i 's/get_copier/getCopier/g' $i
#     sed -i 's/get_reverse_copier/getReverseCopier/g' $i
#     sed -i 's/make_tiles/makeTiles/g' $i
#     sed -i 's/get_realms/getRealms/g' $i
#     sed -i 's/get_box_sorting/getBoxSorting/g' $i
#     sed -i 's/m_boxsort/m_boxSort/g' $i
#     sed -i 's/m_gridgen/m_gridGenerationMethod/g' $i
#     sed -i 's/m_mfis/m_multifluidIndexSpace/g' $i
#     sed -i 's/m_stencil_type/m_stencilType/g' $i
#     sed -i 's/m_centroid_stencil/m_centroidStencilType/g' $i
#     sed -i 's/m_eb_stencil/m_ebCentroidStencilType/g' $i
#     sed -i 's/m_interp_type/m_ghostCellInterpolationMethod/g' $i
#     sed -i 's/m_num_cells/m_numCells/g' $i
#     sed -i 's/m_fill_ratio/m_fillRatioBR/g' $i
#     sed -i 's/m_prob_lo/m_probLo/g' $i
#     sed -i 's/a_prob_lo/a_probLo/g' $i
#     sed -i 's/m_prob_Hi/m_probHi/g' $i
#     sed -i 's/a_prob_Hi/a_probHi/g' $i
#     sed -i 's/m_ref_ratios/m_refinementRatios/g' $i
#     sed -i 's/m_ref_ratio/m_refRatio/g' $i
#     sed -i 's/m_finest_level/m_finestLevel/g' $i
#     sed -i 's/m_max_amr_depth/m_maxAmrDepth/g' $i
#     sed -i 's/m_max_sim_depth/m_maxSimulationDepth/g' $i
#     sed -i 's/m_refine_all_depth/m_refineAllLevelsDepth/g' $i
#     sed -i 's/m_max_box_size/m_maxBoxSize/g' $i
#     sed -i 's/m_max_ebis_box_size/m_maxEbisBoxSize/g' $i
#     sed -i 's/m_buffer_size/m_bufferSizeBR/g' $i
#     sed -i 's/m_irreg_growth/m_irregTagGrowth/g' $i
#     sed -i 's/m_blocking_factor/m_blockingFactor/g' $i
#     sed -i 's/m_ebghost/m_numEbGhostsCells/g' $i
#     sed -i 's/m_num_ghost/m_numGhostCells/g' $i
#     sed -i 's/m_lsf_ghost/m_numLsfGhostCells/g' $i
#     sed -i 's/m_redist_rad/m_redistributionRadius/g' $i
#     sed -i 's/m_centroid_sten_order/m_centroidStencilOrder/g' $i
#     sed -i 's/m_centroid_sten_radius/m_centroidStencilRadius/g' $i
#     sed -i 's/m_eb_sten_order/m_ebCentroidStencilOrder/g' $i
#     sed -i 's/m_eb_sten_radius/m_ebCentroidStencilRadius/g' $i
#     sed -i 's/m_irreg_sten_order/m_irregStenOrder/g' $i
#     sed -i 's/m_irreg_sten_radius/m_irregStenRadius/g' $i
#     sed -i 's/m_ebcf/m_hasEbCf/g' $i
#     sed -i 's/m_ebcf/m_hasEbCf/g' $i
#     sed -i 's/m_has_grids/m_hasGrids/g' $i
#     sed -i 's/define_realms/defineRealms/g' $i
# done


for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/define_realms/defineRealms/g' $i
    sed -i 's/parse_options/parseOptions/g' $i
    sed -i 's/parse_runtime_options/parseRuntimeOptions/g' $i
    sed -i 's/interpGhost_quad/interpGhostQuad/g' $i
    sed -i 's/interpGhost_pwl/interpGhostPwl/g' $i
    sed -i 's/InterpToCentroids/interpToCentroids/g' $i
    sed -i 's/box_sorting/BoxSorting/g' $i
    sed -i 's/grid_generation/GridGenerationMethod/g' $i
done
    
