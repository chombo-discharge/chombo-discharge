# # Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include <cdr_solverF_F.H>/\#include <CD_CdrSolverF_F.H>/g' $i
#     sed -i 's/\#include \"cdr_solver.H\"/\#include \"CD_CdrSolver.H\"/g' $i
#     sed -i 's/\#include \"cdr_solverF_F.H\"/\#include <CD_CdrSolverF_F.H>/g' $i
#     sed -i 's/\#include <cdr_solver.H>/\#include <CD_CdrSolver.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/cdr_solver/CdrSolver/g' $i

# done

# # Move files
# mv src/cdr_solver/cdr_solver.cpp     src/cdr_solver/CD_CdrSolver.cpp
# mv src/cdr_solver/cdr_solver.H       src/cdr_solver/CD_CdrSolver.H
# mv src/cdr_solver/cdr_solverF_F.H    src/cdr_solver/CD_CdrSolverF_F.H
# mv src/cdr_solver/cdr_solverF.ChF    src/cdr_solver/CD_CdrSolverF.ChF
# mv src/cdr_solver/cdr_solver.options src/cdr_solver/CD_CdrSolver.options


# # Move folder
# mv src/cdr_solver src/CdrSolver

# # 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/cdr_solver/CdrSolver/g' $i
#     sed -i 's/advance_tga/advanceTGA/g' $i
#     sed -i 's/advance_euler/advanceEuler/g' $i
#     sed -i 's/compute_div/computeDiv/g' $i
#     sed -i 's/GWN_diffusion_source/gwnDiffusionSource/g' $i
#     sed -i 's/smooth_heaviside_faces/smoothHeavisideFaces/g' $i
#     sed -i 's/fill_GWN/fillGwn/g' $i
#     sed -i 's/average_velo_to_faces/averageVelocityToFaces/g' $i
#     sed -i 's/set_Realm/setRealm/g' $i
#     sed -i 's/set_species/setSpecies/g' $i
#     sed -i 's/a_Realm/a_realm/g' $i
#     sed -i 's/set_domain_bc/setDomainBc/g' $i
#     sed -i 's/set_phase/setPhase/g' $i
#     sed -i 's/set_verbosity/setVerbosity/g' $i
#     sed -i 's/set_time/setTime/g' $i
#     sed -i 's/set_velocity/setVelocity/g' $i
#     sed -i 's/set_diffco/setDiffusionCoefficient/g' $i
#     sed -i 's/set_source/setSource/g' $i
#     sed -i 's/set_ebflux/setEbFlux/g' $i
#     sed -i 's/set_domain_flux/setDomainFlux/g' $i
#     sed -i 's/initialData_distribution/initialDataDistribution/g' $i
#     sed -i 's/initialData_particles/initialDataParticles/g' $i
#     sed -i 's/inject_ebflux/injectEbFlux/g' $i
#     sed -i 's/write_data/writeData/g' $i
#     sed -i 's/reset_redist_weights/setRedistWeights/g' $i
#     sed -i 's/get_Realm/getRealm/g' $i
#     sed -i 's/get_name/getName/g' $i
#     sed -i 's/get_plotVariablesNames/getPlotVariableNames/g' $i
#     sed -i 's/query_ghost/queryGhost/g' $i
#     sed -i 's/get_num_plotvars/getNumberOfPlotVariables/g' $i
#     sed -i 's/compute_advection_dt/computeAdvectionDt/g' $i
#     sed -i 's/compute_diffusion_dt/computeDiffusionDt/g' $i
#     sed -i 's/compute_advection_diffusion_dt/computeAdvectionDiffusionDt/g' $i
#     sed -i 's/compute_source_dt/computeSourceDt/g' $i
#     sed -i 's/compute_mass/computeMass/g' $i
#     sed -i 's/compute_charge/computeCharge/g' $i
#     sed -i 's/is_diffusive/isDiffusive/g' $i
#     sed -i 's/is_mobile/isMobile/g' $i
#     sed -i 's/extrap_source/extrapolateSourceTerm/g' $i
#     sed -i 's/get_state/getPhi/g' $i
#     sed -i 's/get_source/getSource/g' $i
#     sed -i 's/get_velo_cell/getCellCenteredVelocity/g' $i
#     sed -i 's/get_velo_face/getFaceCenteredVelocity/g' $i
#     sed -i 's/get_velo_eb/getEbCenteredVelocity/g' $i
#     sed -i 's/get_diffco_face/getFaceCenteredDiffusionCoefficient/g' $i
#     sed -i 's/get_eb_face/getEbCenteredDiffusionCoefficient/g' $i
#     sed -i 's/get_ebflux/getEbFlux/g' $i
#     sed -i 's/get_domainflux/getDomainFlux/g' $i
#     sed -i 's/m_interp_stencils/m_interpStencils/g' $i
#     sed -i 's/m_interp_sets/m_interpSets/g' $i
#     sed -i 's/m_dombc/m_domainBc/g' $i
#     sed -i 's/m_class_name/m_className/g' $i
#     sed -i 's/m_state/m_phi/g' $i
#     sed -i 's/m_velo_cell/m_cellVelocity/g' $i
#     sed -i 's/m_face_states/m_faceStates/g' $i
#     sed -i 's/m_divG_nc/m_nonConservativeDivG/g' $i
#     sed -i 's/m_mass_diff/m_massDifference/g' $i
#     sed -i 's/m_eb_zero/m_EbZero/g' $i
#     sed -i 's/m_cache_state/m_cachePhi/g' $i
#     sed -i 's/m_cache_source/m_cacheSource/g' $i
#     sed -i 's/m_velo_face/m_faceVelocity/g' $i
#     sed -i 's/m_velo_eb/m_ebVelocity/g' $i
#     sed -i 's/m_ebflux/m_ebFlux/g' $i
#     sed -i 's/m_domainflux/m_domainFlux/g' $i
#     sed -i 's/m_aco/m_aCoefficient/g' $i
#     sed -i 's/m_diffco/m_faceCenteredDiffusionCoefficient/g' $i
#     sed -i 's/m_diffco_eb/m_ebCenteredDiffusionCoefficient/g' $i
#     sed -i 's/m_step/m_timeStep/g' $i
#     sed -i 's/m_redist_mass_weighted/m_useMassWeightedRedistribution/g' $i
#     sed -i 's/m_blend_conservation/m_blendConservation/g' $i
#     sed -i 's/m_diffusive/m_isDiffusive/g' $i
#     sed -i 's/m_mobile/m_isMobile/g' $i
#     sed -i 's/m_extrap_source/m_useSourceExtrapolation/g' $i
#     sed -i 's/m_plot_phi/m_plotPhi/g' $i
#     sed -i 's/m_plot_vel/m_plotVelocity/g' $i
#     sed -i 's/m_plot_dco/m_plotDiffusionCoefficient/g' $i
#     sed -i 's/m_plot_ebf/m_plotEbFlux/g' $i
#     sed -i 's/m_plot_src/m_plotSource/g' $i
#     sed -i 's/m_plot_numbers/m_plotNumbers/g' $i
#     sed -i 's/compute_flux/computeFlux/g' $i
#     sed -i 's/compute_diffusion_flux/computeDiffusionFlux/g' $i
#     sed -i 's/conservative_divergence/conservativeDivergenceNoKappaDivision/g' $i
#     sed -i 's/conservative_divergence_eb/conservativeDivergenceNoKappaDivisionAndNoEbFlux/g' $i
#     sed -i 's/nonconservative_divergence/nonConservativeDivergence/g' $i
#     sed -i 's/hybrid_divergence/hybridDivergence/g' $i
#     sed -i 's/setup_flux_interpolant/setupFluxInterpolant/g' $i
#     sed -i 's/interpolate_flux_to_centroids/interpolateFluxToFaceCentroids/g' $i
#     sed -i 's/computeDivG_irreg/computeDivergenceIrregular/g' $i
#     sed -i 's/increment_flux_register/incrementFluxRegister/g' $i
#     sed -i 's/rest_flux_register/resetFluxRegister/g' $i
#     sed -i 's/coarse_fine_increment/coarseFineIncrement/g' $i
#     sed -i 's/hyperbolic_redistribution/hyperbolicRedistribution/g' $i
#     sed -i 's/increment_redist_flux/incrementRedistFlux/g' $i
#     sed -i 's/increment_redist/incrementRedist/g' $i
#     sed -i 's/coarse_fine_redistribution/coarseFineRedistribution/g' $i
#     sed -i 's/set_ebis/setEbIndexSpace/g' $i
#     sed -i 's/define_interp_stencils/defineInterpolationStencils/g' $i
#     sed -i 's/define_interpolant/defineInterpolant/g' $i
#     sed -i 's/parseDomain_bc/parseDomainBc/g' $i
#     sed -i 's/parse_rng_seed/parseRngSeed/g' $i
#     sed -i 's/parse_plotmode/parsePlotMode/g' $i
#     sed -i 's/parse_conservation/parseDivergenceComputation/g' $i
#     sed -i 's/parse_extrap_source/parseSourceExtrapolation/g' $i
#     sed -i 's/a_new_state/a_newPhi/g' $i
#     sed -i 's/a_old_state/a_oldPhi/g' $i
#     sed -i 's/a_extrap_dt/a_extrapDt/g' $i
#     sed -i 's/a_ebflux/a_ebFlux/g' $i
#     sed -i 's/a_diffusive_term/a_divD/g' $i
#     sed -i 's/a_ransource/a_noiseSource/g' $i
#     sed -i 's/a_cell_states/a_cellPhi/g' $i
#     sed -i 's/a_face_states/a_facePhi/g' $i
#     sed -i 's/a_diffco/a_diffusionCoefficient/g' $i
#     sed -i 's/a_ebG/a_ebFlux/g' $i
#     sed -i 's/a_domain_flux/a_domainFlux/g' $i
#     sed -i 's/a_velo_face/a_faceVelocity/g' $i
#     sed -i 's/a_velo_cell/a_cellVelocity/g' $i
#     sed -i 's/a_face_state/a_facePhi/g' $i
#     sed -i 's/a_face_vel/a_faceVelocity/g' $i
#     sed -i 's/a_cons_div/a_conservativeDivergence/g' $i
#     sed -i 's/a_consdiv/a_conservativeDivergence/g' $i
#     sed -i 's/a_div_nc/a_nonConservativeDivergence/g' $i
#     sed -i 's/a_mass_diff/a_massDifference/g' $i
#     sed -i 's/a_divF_nc/a_nonConservativeDivergence/g' $i
#     sed -i 's/a_divF_H/a_hybridDivergence/g' $i
#     sed -i 's/a_del_vel_rho/a_stableDivergence/g' $i
#     sed -i 's/a_state/a_phi/g' $i
#     sed -i 's/a_diffusionCoefficient_eb/a_ebDiffusionCoefficient/g' $i
#     sed -i 's/parse_extrapolateSourceTerm/parseExtrapolateSourceTerm/g' $i
# done

# # Update makefiles
# for i in `find . -name "*GNUmakefile" -type f`; do
#     sed -i 's/cdr_solver/CdrSolver/g' $i
# done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    # sed -i 's/m_Realm/m_realm/g' $i
    # sed -i 's/m_EbZero/m_ebZero/g' $i
    # sed -i 's/m_faceCenteredDiffusionCoefficient_eb/m_ebCenteredDiffusionCoefficient/g' $i
    # sed -i 's/m_faceCenteredDiffusionCoefficient_eb/m_ebCenteredDiffusionCoefficient/g' $i
    # sed -i 's/conservativeDivergenceNoKappaDivision_eb/conservativeDivergenceNoKappaDivisionOnlyEbFlux/g' $i
    # sed -i 's/nonconservativeDivergenceNoKappaDivision/nonConservativeDivergence/g' $i
    # sed -i 's/a_NonConservativeDivergenceStencil/a_nonConservativeDivergence/g' $i
    sed -i 's/consdiv_regular/conservativeDivergenceRegular/g' $i
done
