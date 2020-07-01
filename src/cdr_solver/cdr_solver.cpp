/*!
  @file   cdr_solver.cpp
  @brief  Implementation of cdr_solver.H
  @author Robert Marskar
  @date   Nov. 2017
  @todo   The diffusive dt computations use a faceiterator box. This should be replaced by fortran routines (or internal ebcell fa functions)
*/

#include "cdr_solver.H"
#include "cdr_solverF_F.H"
#include "cdr_fhdF_F.H"
#include "data_ops.H"

#include <ParmParse.H>
#include <EBAMRIO.H>
#include <EBArith.H>
#include <EBAlias.H>

cdr_solver::cdr_solver(){
  m_name       = "cdr_solver";
  m_class_name = "cdr_solver";
}

cdr_solver::~cdr_solver(){

}

std::string cdr_solver::get_name(){
  return m_name;
}

Vector<std::string> cdr_solver::get_plotvar_names() const {
  CH_TIME("cdr_solver::get_plotvar_names");
  if(m_verbosity > 5){
    pout() << m_name + "::get_plotvar_names" << endl;
  }
  
  Vector<std::string> names(0);
  
  if(m_plot_phi) names.push_back(m_name + " phi");
  if(m_plot_dco && m_diffusive) names.push_back(m_name + " diffusion_coefficient");
  if(m_plot_src) names.push_back(m_name + " source");
  if(m_plot_vel && m_mobile){
    names.push_back("x-Velocity " + m_name);
    names.push_back("y-Velocity " + m_name);
    if(SpaceDim == 3){
      names.push_back("z-Velocity " + m_name);
    }
  }
  if(m_plot_ebf && m_mobile){
    names.push_back(m_name + " eb_flux");
  }
  
  return names;
}

int cdr_solver::query_ghost() const {
  CH_TIME("cdr_solver::query_ghost");
  if(m_verbosity > 5){
    pout() << m_name + "::query_ghost" << endl;
  }

  return 3;
}

int cdr_solver::get_num_plotvars() const {
  CH_TIME("cdr_solver::get_num_plotvars");
  if(m_verbosity > 5){
    pout() << m_name + "::get_num_plotvars" << endl;
  }

  int num_output = 0;

  if(m_plot_phi)                num_output = num_output + 1;
  if(m_plot_dco && m_diffusive) num_output = num_output + 1;
  if(m_plot_src)                num_output = num_output + 1;
  if(m_plot_vel && m_mobile)    num_output = num_output + SpaceDim;
  if(m_plot_ebf && m_mobile)    num_output = num_output + 1;

  return num_output;
}

void cdr_solver::allocate_internals(){
  CH_TIME("cdr_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  CH_assert(m_phase == phase::gas);

  const int sca = 1;
  const int vec = SpaceDim;

  // This is allocated no matter what. 
  m_amr->allocate(m_state,      m_phase, sca);
  m_amr->allocate(m_source,     m_phase, sca);
  m_amr->allocate(m_scratch,    m_phase, sca);
  data_ops::set_value(m_state,      0.0);
  data_ops::set_value(m_source,     0.0);
  data_ops::set_value(m_scratch,    0.0);

  // Only allocate memory for cell-centered and face-centered velocities if the solver is mobile. Otherwise, allocate
  // a NULL pointer that we can pass around in time_stepper in order to handle special cases
  if(m_mobile){
    m_amr->allocate(m_velo_face,     m_phase, sca);
    m_amr->allocate(m_velo_cell,     m_phase, vec);
    m_amr->allocate(m_face_states, m_phase, sca);
    
    data_ops::set_value(m_velo_face,  0.0);
    data_ops::set_value(m_velo_cell,  0.0);
  }
  else{
    m_amr->allocate_ptr(m_velo_face);
    m_amr->allocate_ptr(m_velo_cell);
  }

  // Only allocate memory for diffusion coefficients if we need it. Otherwise, allocate a NULL pointer that we can
  // pass around in time_stepper in order to handle special cases
  if(m_diffusive){
    m_amr->allocate(m_aco,        m_phase, sca);
    m_amr->allocate(m_diffco,     m_phase, sca);
    m_amr->allocate(m_diffco_eb,  m_phase, sca);
    
    data_ops::set_value(m_aco,        0.0);
    data_ops::set_value(m_diffco,     0.0);
    data_ops::set_value(m_diffco_eb,  0.0);
  }
  else{
    m_amr->allocate_ptr(m_aco);
    m_amr->allocate_ptr(m_diffco);
    m_amr->allocate_ptr(m_diffco_eb);
  }

  // Allocate stuff for holding fluxes
  if(m_diffusive || m_mobile){
    m_amr->allocate(m_scratchFluxOne, m_phase, sca);
    m_amr->allocate(m_scratchFluxTwo, m_phase, sca);
  }

  // These don't consume (much) memory so just allocate them 
  m_amr->allocate(m_ebflux,     m_phase, sca);
  m_amr->allocate(m_eb_zero,    m_phase, sca);
  m_amr->allocate(m_domainflux, m_phase, sca);
  m_amr->allocate(m_mass_diff,  m_phase, sca);
  m_amr->allocate(m_divG_nc,    m_phase, sca);
  
  data_ops::set_value(m_ebflux,     0.0);
  data_ops::set_value(m_eb_zero,     0.0);
  data_ops::set_value(m_domainflux, 0.0);

  // This defines interpolation stencils and space for interpolants
  this->define_interp_stencils();
  this->define_interpolant();
}

void cdr_solver::deallocate_internals(){
  CH_TIME("cdr_solver::deallocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::deallocate_internals" << endl;
  }

  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
  m_amr->deallocate(m_velo_face);
  m_amr->deallocate(m_velo_cell);
  m_amr->deallocate(m_ebflux);
  m_amr->deallocate(m_diffco);
  m_amr->deallocate(m_diffco_eb);
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_scratchFluxOne);
  m_amr->deallocate(m_scratchFluxTwo);
  m_amr->deallocate(m_face_states);
}

void cdr_solver::average_velo_to_faces(){
  CH_TIME("cdr_solver::average_velo_to_faces(public, full)");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces(public, full)" << endl;
  }

  this->average_velo_to_faces(m_velo_face, m_velo_cell); // Average velocities to face centers for all levels
}

void cdr_solver::average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell){
  CH_TIME("cdr_solver::average_velo_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::average_cell_to_face(*a_velo_face[lvl], *a_velo_cell[lvl], m_amr->get_domains()[lvl]);
    a_velo_face[lvl]->exchange();
  }
}

void cdr_solver::pre_regrid(const int a_lmin, const int a_old_finest_level){
  CH_TIME("cdr_solver::pre_regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::pre_regrid" << endl;
  }

  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  
  m_amr->allocate(m_cache_state,  m_phase, ncomp);
  m_amr->allocate(m_cache_source, m_phase, ncomp);
  
  for (int lvl = 0; lvl <= a_old_finest_level; lvl++){
    m_state[lvl]->localCopyTo(*m_cache_state[lvl]);
    m_source[lvl]->localCopyTo(*m_cache_source[lvl]);
  }
}

void cdr_solver::coarse_fine_increment(const EBAMRIVData& a_mass_diff){
  CH_TIME("cdr_solver::coarse_fine_increment");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_increment" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_phase)[lvl];

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < 0;

    if(has_coar){
      fine2coar_redist->setToZero();

    }
    if(has_fine){
      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      if(has_coar){
	fine2coar_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
      }

      if(has_fine){
	coar2fine_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
	coar2coar_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
      }
    }
  }
}

void cdr_solver::coarse_fine_redistribution(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::coarse_fine_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_redistribution" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(m_phase)[lvl];
    
    if(has_coar){
      fine2coar_redist->redistribute(*a_state[lvl-1], interv);
      fine2coar_redist->setToZero();
    }

    if(has_fine){
      coar2fine_redist->redistribute(*a_state[lvl+1], interv);
      coar2coar_redist->redistribute(*a_state[lvl],   interv);

      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }
  }
}

void cdr_solver::compute_divG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebG){
  CH_TIME("cdr_solver::compute_divG");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divG" << endl;
  }

  data_ops::set_value(a_divG, 0.0);
  
  this->conservative_divergence(a_divG, a_G, a_ebG);       // Make the conservative divergence.
  this->nonconservative_divergence(m_divG_nc, a_divG);     // Non-conservative divergence
  this->hybrid_divergence(a_divG, m_mass_diff, m_divG_nc); // a_divG becomes hybrid divergence. Mass diff computed. 
  this->increment_flux_register(a_G);                      // Increment flux register
  this->increment_redist(m_mass_diff);                     // Increment level redistribution register

  const bool ebcf = m_amr->get_ebcf();
  if(ebcf){ // If we have EBCF, much more work with the CF interface
    this->coarse_fine_increment(m_mass_diff);             // Compute C2F, F2C, and C2C mass transfers
    this->increment_redist_flux();                        // Tell flux register about whats going on
    this->hyperbolic_redistribution(a_divG, m_mass_diff); // Level redistribution. Weights is a dummy parameter
    this->coarse_fine_redistribution(a_divG);             // Do the coarse-fine redistribution
    this->reflux(a_divG);                                 // Reflux
  }
  else{ // Much simpler if we don't have EBCF
    this->hyperbolic_redistribution(a_divG, m_mass_diff); // Level redistribution. Weights is a dummy parameter
    this->reflux(a_divG);                                 // Reflux
  }
}

void cdr_solver::inject_ebflux(EBAMRCellData& a_phi, const EBAMRIVData& a_ebG, const Real a_dt){
  CH_TIME("cdr_solver::inject_ebflux");
  if(m_verbosity > 5){
    pout() << m_name + "::inject_ebflux" << endl;
  }

  MayDay::Warning("cdr_solver::inject_ebflux - routine has not been wetted!");

  if(m_redist_mass_weighted){
    this->reset_redist_weights(a_phi);
  }

  this->conservative_divergence_eb(m_scratch, a_ebG);         // Compute conservative divergence, but only EB
  this->nonconservative_divergence(m_divG_nc, m_scratch);     // Blend with volume fraction
  this->hybrid_divergence(m_scratch, m_mass_diff, m_divG_nc); // Hybrid divergence
  this->increment_redist(m_mass_diff);                        // Increment redistribution register

  const bool ebcf = m_amr->get_ebcf();
  if(ebcf){ // Much more work with the EBCF interface. *Sigh*
    this->coarse_fine_increment(m_mass_diff);                // Compute C2F, F2C, and C2C mass transfers
    this->increment_redist_flux();                           // Tell flux register about whats going on
    this->hyperbolic_redistribution(m_scratch, m_mass_diff); // Level redistribution. 
    this->coarse_fine_redistribution(m_scratch);             // Do the coarse-fine redistribution
    this->reflux(m_scratch);                                 // Reflux
  }
  else{
    this->hyperbolic_redistribution(m_scratch, m_mass_diff);
  }

  // Now do the increment. 
  data_ops::incr(a_phi, m_scratch, -a_dt);
}

void cdr_solver::conservative_divergence_eb(EBAMRCellData& a_consdiv, const EBAMRIVData& a_ebflux){
  CH_TIME("cdr_solver::inject_ebflux");
  if(m_verbosity > 5){
    pout() << m_name + "::inject_ebflux" << endl;
  }

  // TLDR: This sets a_consdiv = a_ebflux*area/dx

  const int comp = 0;

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& divG            = (*a_consdiv[lvl])[dit()];
      const EBISBox& ebisbox     = ebisl[dit()];
      const BaseIVFAB<Real>& flx = (*a_ebflux[lvl])[dit()];

      divG.setVal(0.0);
      VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real area     = ebisbox.bndryArea(vof);

	divG(vof, comp) = flx(vof, comp)*area/dx;
      }
    }
  }
}

void cdr_solver::compute_divG_irreg(LevelData<EBCellFAB>&              a_divG,
				    const LevelData<BaseIVFAB<Real> >& a_ebflux,
				    const int                          a_lvl){
  CH_TIME("cdr_solver::compute_divG_irreg");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divG_irreg" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];
  const Real dx                = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];

    EBCellFAB& divG               = a_divG[dit()];
    const BaseIVFAB<Real>& ebflux = a_ebflux[dit()];

    VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real area     = ebisbox.bndryArea(vof);

      // EB flux
      divG(vof, comp) = ebflux(vof,comp)*area;

      // Face fluxes
      for (int dir = 0; dir < SpaceDim; dir++){

	const BaseIFFAB<Real>& flux = (*(m_interpolant[dir])[a_lvl])[dit()];

	for (SideIterator sit; sit.ok(); ++sit){
	  const int isign = sign(sit());
	  const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	  for (int iface = 0; iface < faces.size(); iface++){
	    const FaceIndex face = faces[iface];
	    const Real face_area = ebisbox.areaFrac(face);
	    divG(vof, comp) += isign*face_area*flux(face, comp);
	  }
	}
      }

      // Scale divG by dx but not by kappa.
      divG(vof, comp) *= 1./dx;
    }
  }
}

void cdr_solver::compute_flux(EBAMRFluxData&       a_flux,
			      const EBAMRFluxData& a_face_state,
			      const EBAMRFluxData& a_face_vel,
			      const EBAMRIFData&   a_domain_flux){
  CH_TIME("cdr_solver::compute_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_flux" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){

#if 1 // New code
    compute_flux(*a_flux[lvl], *a_face_state[lvl], *a_face_vel[lvl], *a_domain_flux[lvl], lvl);
#else // Old code (that we know works)
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& flx       = (*a_flux[lvl])[dit()][dir];
	const EBFaceFAB& phi = (*a_face_state[lvl])[dit()][dir];
	const EBFaceFAB& vel = (*a_face_vel[lvl])[dit()][dir];

	flx.setVal(0.0, comp);
	flx += phi;
	flx *= vel;

	// Irregular faces
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;
	for (FaceIterator faceit(ebisbox.getIrregIVS(box), ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  flx(face, comp) = vel(face, comp)*phi(face, comp);
	}

	// Domain faces
	for (SideIterator sit; sit.ok(); ++sit){
	  BaseIFFAB<Real>& domflux = (*a_domain_flux[lvl])[dit()](dir, sit());

	  const IntVectSet& ivs  = domflux.getIVS();
	  const EBGraph& ebgraph = domflux.getEBGraph();

	  const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
	  for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	    const FaceIndex& face = faceit();
	    
	    if(m_dombc == cdr_bc::external){
	      flx(face, comp) = domflux(face, comp);
	    }
	    else if(m_dombc == cdr_bc::wall){
	      flx(face, comp) = 0.0;
	    }
	    else if(m_dombc == cdr_bc::outflow){

	      flx(face, comp) = Max(0.0, sign(sit())*flx(face, comp));
	    }
	    else if(m_dombc == cdr_bc::extrap){
	      // Don't do anything, the solver should have extrapolated the face-centered state
	    }
	    else {
	      MayDay::Abort("cdr_solver::compute_flux - stop this madness!");
	    }
	  }
	}
      }
    }
#endif
  }
}

void cdr_solver::compute_flux(LevelData<EBFluxFAB>&              a_flux,
			      const LevelData<EBFluxFAB>&        a_face_state,
			      const LevelData<EBFluxFAB>&        a_face_vel,
			      const LevelData<DomainFluxIFFAB>&  a_domain_flux,
			      const int                          a_lvl){

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flx       = a_flux[dit()][dir];
      const EBFaceFAB& phi = a_face_state[dit()][dir];
      const EBFaceFAB& vel = a_face_vel[dit()][dir];

      flx.setVal(0.0, comp);
      flx += phi;
      flx *= vel;

      // Irregular faces
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ebisbox.getIrregIVS(box), ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	flx(face, comp) = vel(face, comp)*phi(face, comp);
      }

      // Domain faces
      for (SideIterator sit; sit.ok(); ++sit){
	const BaseIFFAB<Real>& domflux = a_domain_flux[dit()](dir, sit());

	const IntVectSet& ivs  = domflux.getIVS();
	const EBGraph& ebgraph = domflux.getEBGraph();

	const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
	for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	    
	  if(m_dombc == cdr_bc::external){
	    flx(face, comp) = domflux(face, comp);
	  }
	  else if(m_dombc == cdr_bc::wall){
	    flx(face, comp) = 0.0;
	  }
	  else if(m_dombc == cdr_bc::outflow){
	    flx(face, comp) = Max(0.0, sign(sit())*flx(face, comp));
	  }
	  else if(m_dombc == cdr_bc::extrap){
	    // Don't do anything, the solver should have extrapolated the face-centered state
	  }
	  else {
	    MayDay::Abort("cdr_solver::compute_flux - stop this madness!");
	  }
	}
      }
    }
  }
}

void cdr_solver::compute_diffusion_flux(EBAMRFluxData& a_flux, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::compute_diffusion_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_flux" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    this->compute_diffusion_flux(*a_flux[lvl], *a_state[lvl], lvl);
  }
}



void cdr_solver::compute_diffusion_flux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_state, const int a_lvl){
  CH_TIME("cdr_solver::compute_diffusion_flux(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_flux(level)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const Real dx                = m_amr->get_dx()[a_lvl];
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& cellbox = dbl.get(dit());

    const EBCellFAB& state = a_state[dit()];
    const EBISBox& ebisbox = state.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux = a_flux[dit()][dir];
      const EBFaceFAB& dco = (*m_diffco[a_lvl])[dit()][dir];

      // Do regular cells
      Box facebox = cellbox;
      facebox.surroundingNodes(dir);

      BaseFab<Real>& regFlux        = flux.getSingleValuedFAB();
      const BaseFab<Real>& regState = state.getSingleValuedFAB();
      const BaseFab<Real>& regDco   = dco.getSingleValuedFAB();
      FORT_DFLUX_REG(CHF_FRA1(regFlux, comp),
		     CHF_CONST_FRA1(regState, comp),
		     CHF_CONST_FRA1(regDco, comp),
		     CHF_CONST_INT(dir),
		     CHF_CONST_REAL(dx),
		     CHF_BOX(facebox));


      // Now redo the irregular faces
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ebisbox.getIrregIVS(cellbox), ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();

	if(face.isBoundary()){ // No boundary flux for diffusion stuff. Could be changed but that's the current status. 
	  flux(face, comp) = 0.0;
	}
	else{

	  const VolIndex hiVoF = face.getVoF(Side::Hi);
	  const VolIndex loVoF = face.getVoF(Side::Lo);

	  flux(face, comp) = dco(face,comp)*(state(hiVoF,comp) - state(loVoF, comp))/dx;
	}
      }
    }
  }
  
}

void cdr_solver::conservative_divergence(EBAMRCellData& a_cons_div, EBAMRFluxData& a_flux, const EBAMRIVData& a_ebflux){
  CH_TIME("cdr_solver::conservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::conservative_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_flux[lvl]->exchange();
    
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    this->consdiv_regular(*a_cons_div[lvl], *a_flux[lvl], lvl);
    this->setup_flux_interpolant(*a_flux[lvl], lvl);                 // Copy face-centered fluxes in a_flux to m_interpolant
    this->interpolate_flux_to_centroids(*a_flux[lvl], lvl);          // Interpolate fluxes w m_interpolant. Copy 2 a_flux
    this->compute_divG_irreg(*a_cons_div[lvl], *a_ebflux[lvl], lvl); // Recompute divergence on irregular cells

    a_cons_div[lvl]->exchange();
  }
}

void cdr_solver::consdiv_regular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("cdr_solver::consdiv_regular");
  if(m_verbosity > 5){
    pout() << m_name + "::consdiv_regular" << endl;
  }

  CH_assert(a_divJ.nComp() == 1);
  CH_assert(a_flux.nComp() == 1);

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain domain  = m_amr->get_domains()[a_lvl];
  const Real dx               = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& divJ         = a_divJ[dit()];
    BaseFab<Real>& divJ_fab = divJ.getSingleValuedFAB();
    const Box box = dbl.get(dit());

    divJ.setVal(0.0);
    for (int dir = 0; dir < SpaceDim; dir++){
      const EBFaceFAB& flx         = a_flux[dit()][dir];
      const BaseFab<Real>& flx_fab = flx.getSingleValuedFAB();

      ParmParse pp ("cache_blocking");

      int block_loops = 0;
      pp.query("tile_loops", block_loops);
      if(block_loops==0){
	FORT_CONSDIV_REG(CHF_FRA1(divJ_fab, comp),
			 CHF_CONST_FRA1(flx_fab, comp),
			 CHF_CONST_INT(dir),
			 CHF_CONST_REAL(dx),
			 CHF_BOX(box));
      }
      else{
	Vector<int> tilesize(SpaceDim);
	pp.getarr("tile_size", tilesize, 0, SpaceDim);
	pout() << tilesize << endl;
	Vector<Box> boxes = m_amr->make_tiles(box, IntVect(D_DECL(tilesize[0], tilesize[1], tilesize[2])));
	for (int ibox = 0; ibox < boxes.size(); ibox++){
	  FORT_CONSDIV_REG(CHF_FRA1(divJ_fab, comp),
			   CHF_CONST_FRA1(flx_fab, comp),
			   CHF_CONST_INT(dir),
			   CHF_CONST_REAL(dx),
			   CHF_BOX(boxes[ibox]));
	}
      }
    }

    VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      divJ(vof, comp) = 0.0;
    }
  }

  a_divJ.exchange();
}

void cdr_solver::define_interp_stencils(){
  CH_TIME("cdr_solver::define_interp_stencils");
  if(m_verbosity > 5){
    pout() << m_name + "::define_interp_stencils" << endl;
  }

  const int comp            = 0;
  const int ncomp           = 1;
  const int finest_level    = m_amr->get_finest_level();
  FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    (m_interp_stencils[dir]).resize(1 + finest_level);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const ProblemDomain& domain  = m_amr->get_domains()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

      m_interp_stencils[dir][lvl] = RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > >
	(new LayoutData<BaseIFFAB<FaceStencil> >(dbl));

      LayoutData<IntVectSet> cfivs(dbl);
      EBArith::defineCFIVS(cfivs, dbl, domain);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BaseIFFAB<FaceStencil>& sten = (*m_interp_stencils[dir][lvl])[dit()];
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet ivs   = ebisbox.getIrregIVS(box);

	sten.define(ivs, ebisbox.getEBGraph(), dir, ncomp);
	
	for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  FaceStencil& facesten = sten(face,comp);

	  facesten = EBArith::getInterpStencil(face, cfivs[dit()], ebisbox, domain.domainBox());
	}
      }
    }
  }
}

void cdr_solver::increment_redist_flux(){
  CH_TIME("cdr_solver::increment_redist_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist_flux" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;
    
    if(has_fine){
      RefCountedPtr<EBFluxRegister>& fluxreg = m_amr->get_flux_reg(m_phase)[lvl];
      RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_phase)[lvl];
      RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_phase)[lvl];
      
      const Real scale = -dx;

      fluxreg->incrementRedistRegister(*coar2fine_redist, interv, scale);
      fluxreg->incrementRedistRegister(*coar2coar_redist, interv, scale);

    }
  }
}

void cdr_solver::initial_data(){
  CH_TIME("cdr_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }

  const bool deposit_function  = m_species->init_with_function();
  const bool deposit_particles = m_species->init_with_particles();

  data_ops::set_value(m_state, 0.0);
  
  if(deposit_particles){
    initial_data_particles();
  }

  // Increment with function values if this is also called for
  if(deposit_function){
    initial_data_distribution();
  }

  m_amr->average_down(m_state, m_phase);
  m_amr->interp_ghost(m_state, m_phase);
}

void cdr_solver::initial_data_distribution(){
  CH_TIME("cdr_solver::initial_data_distribution");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data_distribution" << endl;
  }

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  // Copy this
  data_ops::copy(m_scratch, m_state);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx = m_amr->get_dx()[lvl];

    for (DataIterator dit = m_state[lvl]->dataIterator(); dit.ok(); ++dit){
      EBCellFAB& state         = (*m_state[lvl])[dit()];
      const EBCellFAB& scratch = (*m_scratch[lvl])[dit()];
      const Box box            = m_state[lvl]->disjointBoxLayout().get(dit());
      const EBISBox& ebisbox   = state.getEBISBox();

      BaseFab<Real>& reg_state   = state.getSingleValuedFAB();
      const BaseFab<Real>& reg_scratch = scratch.getSingleValuedFAB();

      // Regular cells
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	const RealVect pos = origin + RealVect(iv)*dx + 0.5*dx*RealVect::Unit;

	for (int comp = 0; comp < state.nComp(); comp++){
	  reg_state(iv, comp) = reg_scratch(iv, comp) + m_species->initial_data(pos, m_time);
	}
      }

      // Irreg and multicells
      VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa    = ebisbox.volFrac(vof);
	const RealVect pos  = EBArith::getVofLocation(vof, m_amr->get_dx()[lvl]*RealVect::Unit, origin);
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = scratch(vof, comp) + m_species->initial_data(pos, m_time);
	}
      }
    }
  }

  m_amr->average_down(m_state, m_phase);
  m_amr->interp_ghost(m_state, m_phase);

  data_ops::set_covered_value(m_state, 0, 0.0);
}

void cdr_solver::initial_data_particles(){
  CH_TIME("cdr_solver::initial_data_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data_particles" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  const RealVect origin  = m_amr->get_prob_lo();
  DepositionType::Which deposition  = m_species->get_deposition();

  const int comp = 0;
  const Interval interv(comp, comp);

  // Need buffer for everything that has non-zero cloud width
  int pvr_buffer = 1;
  if(deposition == DepositionType::CIC){
    pvr_buffer = 0;
  }

  const List<Particle>& initial_particles = m_species->get_initial_particles();
  const int num_particles = initial_particles.length();


  // Allocate some stuff
  if(num_particles > 0){
    Vector<RefCountedPtr<ParticleData<Particle> > > amrparticles;
    Vector<RefCountedPtr<ParticleValidRegion> > pvr;
    m_amr->allocate(amrparticles);
    m_amr->allocate(pvr, pvr_buffer);

    // Have not added the particles to the 
    const DisjointBoxLayout& dbl_coar = m_amr->get_grids()[0];
    for (DataIterator dit = dbl_coar.dataIterator(); dit.ok(); ++dit){
      const Box box = dbl_coar.get(dit());
      (*amrparticles[0])[dit()].addItems(m_species->get_initial_particles());
    }

    // Remap amr particles
    for (int lvl = 1; lvl <= m_amr->get_finest_level(); lvl++){

      // 1. Collect coarser level particles into this levels PVR
      collectValidParticles(amrparticles[lvl]->outcast(),
			    *amrparticles[lvl-1],
			    pvr[lvl]->mask(),
			    m_amr->get_dx()[lvl]*RealVect::Unit,
			    m_amr->get_ref_rat()[lvl-1],
			    false, 
			    origin);
      amrparticles[lvl]->remapOutcast();
    }

    // We will deposit onto m_state, using m_scratch as a scratch holder for interpolation stuff
    data_ops::set_value(m_state, 0.0);
    data_ops::set_value(m_scratch, 0.0);

    // Deposit onto mseh
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const Real dx                = m_amr->get_dx()[lvl];
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const ProblemDomain& dom     = m_amr->get_domains()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

      const bool has_coar = (lvl > 0);
      const bool has_fine = (lvl < finest_level);

      // 1. If we have a coarser level whose cloud hangs into this level, interpolate the coarser level here first
      if(has_coar && deposition != DepositionType::NGP){
	RefCountedPtr<EBPWLFineInterp>& interp = m_amr->get_eb_pwl_interp(m_phase)[lvl];
	interp->interpolate(*m_state[lvl], *m_scratch[lvl-1], interv);
      }
    
      // 2. Deposit this levels particles and exchange ghost cells
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	EBParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin);
	interp.deposit((*amrparticles[lvl])[dit()].listItems(), (*m_state[lvl])[dit()].getFArrayBox(), deposition);
      }


      // Exchange ghost cells
      const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_phase)[lvl];
      LDaddOp<FArrayBox> addOp;
      LevelData<FArrayBox> aliasFAB;
      aliasEB(aliasFAB, *m_state[lvl]);
      aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

      // 3. If we have a finer level, copy contributions from this level to the temporary holder that is used for
      //    interpolation of "hanging clouds"
      if(has_fine){
	m_state[lvl]->localCopyTo(*m_scratch[lvl]);
      }
    }

#if CH_SPACEDIM==2 // Only do this scaling for planar cartesian
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const Real dx = m_amr->get_dx()[lvl];
      data_ops::scale(*m_state[lvl], 1./dx);
    }
#endif
  }
}

void cdr_solver::hybrid_divergence(EBAMRCellData&     a_hybrid_div,
				   EBAMRIVData&       a_mass_diff,
				   const EBAMRIVData& a_noncons_div){
  CH_TIME("cdr_solver::hybrid_divergence(AMR)");
  if(m_verbosity > 5){
    pout() << m_name + "::hybrid_divergence(AMR)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    hybrid_divergence(*a_hybrid_div[lvl], *a_mass_diff[lvl], *a_noncons_div[lvl], lvl);
  }
}

void cdr_solver::hybrid_divergence(LevelData<EBCellFAB>&              a_divF_H,
				   LevelData<BaseIVFAB<Real> >&       a_mass_diff,
				   const LevelData<BaseIVFAB<Real> >& a_divF_nc,
				   const int                          a_lvl){
  CH_TIME("cdr_solver::hybrid_divergence(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::hybrid_divergence(level)" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];
    
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& divH               = a_divF_H[dit()];  // On input, this contains kappa*div(F)
    BaseIVFAB<Real>& deltaM       = a_mass_diff[dit()];
    const BaseIVFAB<Real>& divNC  = a_divF_nc[dit()]; 

    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];

    VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real kappa    = ebisbox.volFrac(vof);
      const Real dc       = divH(vof, comp);
      const Real dnc      = divNC(vof, comp);

      // Note to self: deltaM = (1-kappa)*(dc - kappa*dnc) because dc was not divided by kappa,
      // which it would be otherwise. 
      divH(vof, comp)   = dc + (1-kappa)*dnc;          // On output, contains hybrid divergence
      deltaM(vof, comp) = (1-kappa)*(dc - kappa*dnc);
    }
  }
}

void cdr_solver::reset_redist_weights(const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::reset_redist_weights");
  if(m_verbosity > 5){
    pout() << m_name + "::reset_redist_weights" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& redist = *m_amr->get_level_redist(m_phase)[lvl];
    redist.resetWeights(*a_state[lvl], comp);

    if(m_amr->get_ebcf()){
      const bool has_coar = lvl > 0;
      const bool has_fine = lvl < finest_level;

      if(has_coar){
	EBFineToCoarRedist& fine2coar = *m_amr->get_fine_to_coar_redist(m_phase)[lvl];
	fine2coar.resetWeights(*a_state[lvl], comp);
      }
      if(has_fine){
	EBCoarToCoarRedist& coar2coar = *m_amr->get_coar_to_coar_redist(m_phase)[lvl];
	EBCoarToFineRedist& coar2fine = *m_amr->get_coar_to_fine_redist(m_phase)[lvl];

	coar2coar.resetWeights(*a_state[lvl], comp);
	coar2fine.resetWeights(*a_state[lvl], comp);
      }
    }

  }
}

void cdr_solver::hyperbolic_redistribution(EBAMRCellData& a_divF, const EBAMRIVData&   a_mass_diff) {
  CH_TIME("cdr_solver::hyberbolic_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::hyperbolic_redistribution" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[lvl]);
    level_redist.redistribute(*a_divF[lvl], interv);
    level_redist.setToZero();
  }
}

void cdr_solver::concentration_redistribution(EBAMRCellData& a_phi, const EBAMRIVData& a_mass_diff){
  CH_TIME("cdr_solver::concentration_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::concentration_redistribution" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelConcentrationRedist& redist = *(m_amr->get_concentration_redist(m_phase)[lvl]);
    redist.redistribute(*a_phi[lvl], interv);
    redist.setToZero();
  }
}

void cdr_solver::interpolate_flux_to_centroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("cdr_solver::interpolate_flux_to_centroids");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_flux_to_centroids" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];

  const FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);
    
    for (int dir = 0; dir < SpaceDim; dir++){
      BaseIFFAB<Real>& interpolant = (*(m_interpolant[dir])[a_lvl])[dit()]; // Holds the face-centered flux
      BaseIFFAB<Real> centroid_flux;
      EBFaceFAB& faceFlux = a_flux[dit()][dir];

      // Compute face centroid flux
      centroid_flux.define(ivs, ebgraph, dir, ncomp);
      for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	const FaceIndex& face   = faceit();
	const FaceStencil& sten = (*m_interp_stencils[dir][a_lvl])[dit()](face, comp);

	centroid_flux(face, comp) = 0.;
	Real sum = 0.0;
	for (int i = 0; i < sten.size(); i++){
	  const FaceIndex& iface = sten.face(i);
	  const Real iweight     = sten.weight(i);
	  sum += iweight;
	  
	  centroid_flux(face, comp) += iweight*interpolant(iface, comp);
	}

	faceFlux(face, comp) = centroid_flux(face, comp);
      }

      // Copy centroid flux into a_flux
      interpolant.setVal(0.0);
      interpolant.copy(box, interv, box, centroid_flux, interv);
    }
  }
}

void cdr_solver::reset_flux_register(){
  CH_TIME("cdr_solver::reset_flux_register");
  if(m_verbosity > 5){
    pout() << m_name + "::reset_flux_register" << endl;
  }
  
  const int finest_level = m_amr->get_finest_level();
  
  Vector<RefCountedPtr<EBFluxRegister> >& fluxreg = m_amr->get_flux_reg(m_phase);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const bool has_fine = lvl < finest_level;
    if(has_fine){
      fluxreg[lvl]->setToZero();
    }
  }
  
}

void cdr_solver::increment_flux_register(const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_velo_face){
  CH_TIME("cdr_solver::increment_flux_register");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_flux_register" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister> >& fluxreg = m_amr->get_flux_reg(m_phase);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      fluxreg[lvl]->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const Box box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
	const Real scale     = 1.;
	const EBFaceFAB& phi = (*a_face_state[lvl])[dit()][dir];
	const EBFaceFAB& vel = (*a_velo_face[lvl])[dit()][dir];

	// Compute flux
	EBFaceFAB flux(ebisbox, box, dir, ncomp);
	flux.setVal(0.0);
	flux += phi;
	flux *= vel;

	// Increment flux register for irregular/regular. Add both from coarse to fine and from fine to coarse
	if(has_fine){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl]->incrementCoarseBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
	if(has_coar){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl-1]->incrementFineBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
      }
    }
  }
}

void cdr_solver::increment_flux_register(const EBAMRFluxData& a_flux){
  CH_TIME("cdr_solver::increment_flux_register(flux)");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_flux_register(flux)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister> >& fluxreg = m_amr->get_flux_reg(m_phase);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      fluxreg[lvl]->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const Box box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
	const Real scale      = 1.0;
	const EBFaceFAB& flux = (*a_flux[lvl])[dit()][dir];

	// Increment flux register for irregular/regular. Add both from coarse to fine and from fine to coarse
	if(has_fine){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl]->incrementCoarseBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
	if(has_coar){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl-1]->incrementFineBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
      }
    }
  }
}

void cdr_solver::increment_redist(const EBAMRIVData& a_mass_diff){
  CH_TIME("cdr_solver::increment_redist");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      level_redist.increment((*a_mass_diff[lvl])[dit()], dit(), interv);
    }
  }
}

void cdr_solver::increment_concentration_redist(const EBAMRIVData& a_mass_diff){
  CH_TIME("cdr_solver::increment_redist");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    
    EBLevelConcentrationRedist& redist = *(m_amr->get_concentration_redist(m_phase)[lvl]);
    redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      redist.increment((*a_mass_diff[lvl])[dit()], dit(), interv);
    }
  }
}

void cdr_solver::nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divG){
  CH_TIME("cdr_solver::nonconservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::nonconservative_divergence" << endl;
  }

  if(m_blend_conservation){
    irreg_amr_stencil<noncons_div>& stencils = m_amr->get_noncons_div_stencils(m_phase);
    stencils.apply(a_div_nc, a_divG);
  }
  else{
    data_ops::set_value(a_div_nc, 0.0);
  }
}

void cdr_solver::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("cdr_solver::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocate_internals();

  Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->get_eb_pwl_interp(m_phase);

  // These levels have not changed
  for (int lvl = 0; lvl <= Max(0,a_lmin-1); lvl++){
    m_cache_state[lvl]->copyTo(*m_state[lvl]); // Base level should never change, but ownership can.
    m_cache_source[lvl]->copyTo(*m_source[lvl]); // Base level should never change, but ownership can.
  }

  // These levels have changed
  for (int lvl = a_lmin; lvl <= a_new_finest_level; lvl++){
    interpolator[lvl]->interpolate(*m_state[lvl], *m_state[lvl-1], interv);
    interpolator[lvl]->interpolate(*m_source[lvl], *m_source[lvl-1], interv);
    if(lvl <= Min(a_old_finest_level, a_new_finest_level)){
      m_cache_state[lvl]->copyTo(*m_state[lvl]);
      m_cache_source[lvl]->copyTo(*m_source[lvl]);
    }
  }

  //  data_ops::floor(m_state, 0.0);
  m_amr->average_down(m_state, m_phase);
  m_amr->interp_ghost(m_state, m_phase);

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);

}

void cdr_solver::reflux(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::reflux");
  if(m_verbosity > 5){
    pout() << m_name + "::reflux" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister > >& fluxreg = m_amr->get_flux_reg(m_phase);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      const Real scale = 1.0/dx;
      
      fluxreg[lvl]->reflux(*a_state[lvl], interv, scale);
      fluxreg[lvl]->setToZero();
    }
  }
}

void cdr_solver::sanity_check(){
  CH_TIME("cdr_solver::sanity_check");
  if(m_verbosity > 5){
    pout() << m_name + "::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_species.isNull());
  CH_assert(!m_ebis.isNull());
}

void cdr_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("cdr_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;

}

void cdr_solver::register_operators(){
  CH_TIME("cdr_solver::register_operators");
  if(m_verbosity > 5){
    pout() << m_name + "::register_operators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("cdr_solver::register_operators - need to set amr_mesh!");
  }
  else{
    m_amr->register_operator(s_eb_coar_ave,     m_phase);
    m_amr->register_operator(s_eb_quad_cfi,     m_phase);
    m_amr->register_operator(s_eb_fill_patch,   m_phase);
    m_amr->register_operator(s_eb_pwl_interp,   m_phase);
    m_amr->register_operator(s_eb_flux_reg,     m_phase);
    m_amr->register_operator(s_eb_redist,       m_phase);
    m_amr->register_operator(s_eb_irreg_interp, m_phase);
    m_amr->register_operator(s_eb_noncons_div,  m_phase);
  }
}

void cdr_solver::set_domain_bc(const cdr_bc::which_bc a_bctype){
  CH_TIME("cdr_solver::set_domain_bc");
  if(m_verbosity > 5){
    pout() << m_name + "::set_domain_bc" << endl;
  }

  m_dombc = a_bctype;
}

void cdr_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("cdr_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;

  const RefCountedPtr<mfis> mfis = m_compgeom->get_mfis();
  
  this->set_ebis(mfis->get_ebis(m_phase));
}

void cdr_solver::set_diffco(const EBAMRFluxData& a_diffco, const EBAMRIVData& a_diffco_eb){
  CH_TIME("cdr_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_diffco[lvl]->localCopyTo(*m_diffco[lvl]);
    a_diffco_eb[lvl]->localCopyTo(*m_diffco_eb[lvl]);
  }

  m_amr->average_down(m_diffco, m_phase);
  m_amr->average_down(m_diffco_eb, m_phase);
}

void cdr_solver::set_diffco(const Real a_diffco){
  CH_TIME("cdr_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_diffco[lvl],    a_diffco);
    data_ops::set_value(*m_diffco_eb[lvl], a_diffco);

    m_diffco[lvl]->exchange();
  }

  m_amr->average_down(m_diffco, m_phase);
  m_amr->average_down(m_diffco_eb, m_phase);
}

void cdr_solver::set_ebflux(const EBAMRIVData& a_ebflux){
  CH_TIME("cdr_solver::set_ebflux(variable)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebflux(variable)" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_ebflux[lvl]->localCopyTo(interv, *m_ebflux[lvl], interv);
  }
}

void cdr_solver::set_ebflux(const Real a_ebflux){
  CH_TIME("cdr_solver::set_ebflux(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebflux(constant)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_ebflux[lvl], a_ebflux);
  }
}

void cdr_solver::set_domain_flux(const Real a_domain_flux){
  CH_TIME("cdr_solver::set_domain_flux(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_domain_flux(constant)" << endl;
  }

  data_ops::set_value(m_domainflux, 0.0);
}

void cdr_solver::set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("cdr_solver::set_ebis");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebis" << endl;
  }

  m_ebis = a_ebis;
}

void cdr_solver::set_species(const RefCountedPtr<cdr_species> a_species){
  CH_TIME("cdr_solver::set_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_species" << endl;
  }

  m_species   = a_species;
  m_name      = m_species->get_name();
  m_diffusive = m_species->is_diffusive();
  m_mobile    = m_species->is_mobile();
}

void cdr_solver::set_source(const EBAMRCellData& a_source){
  CH_TIME("cdr_solver::set_source");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->localCopyTo(*m_source[lvl]);
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void cdr_solver::set_source(const Real a_source){
  CH_TIME("cdr_solver::set_source");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source" << endl;
  }

  const int comp = 0;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_source[lvl], a_source, comp);
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void cdr_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("cdr_solver::set_time");
  if(m_verbosity > 5){
    pout() << m_name + "::set_time" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void cdr_solver::set_velocity(const EBAMRCellData& a_velo){
  CH_TIME("cdr_solver::set_velocity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_velo[lvl]->localCopyTo(*m_velo_cell[lvl]);
  }

  m_amr->average_down(m_velo_cell, m_phase);
  m_amr->interp_ghost(m_velo_cell, m_phase);
}

void cdr_solver::set_velocity(const RealVect a_velo){
  CH_TIME("cdr_solver::set_velocity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (int dir = 0; dir < SpaceDim; dir++){
      data_ops::set_value(*m_velo_cell[lvl], a_velo[dir], dir);
    }

    m_velo_cell[lvl]->exchange();
  }

  m_amr->average_down(m_velo_cell, m_phase);
  m_amr->interp_ghost(m_velo_cell, m_phase);
}

void cdr_solver::set_phase(const phase::which_phase a_phase){
  CH_TIME("cdr_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void cdr_solver::set_verbosity(const int a_verbosity){
  CH_TIME("cdr_solver::set_verbosity");
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
}

void cdr_solver::setup_flux_interpolant(const LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("cdr_solver::setup_flux_interpolant");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_flux_interpolant" << endl;
  }

  const int ncomp = 1;
  const int comp  = 0;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];

  for (int dir = 0; dir < SpaceDim; dir++){
    
    const LayoutData<IntVectSet>& grown_set = (*(m_interp_sets[dir])[a_lvl]);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIFFAB<Real>& interpol   = (*(m_interpolant[dir])[a_lvl])[dit()];
      const Box& box              = dbl.get(dit());
      const EBISBox& ebisbox      = ebisl[dit()];
      const EBGraph& ebgraph      = ebisbox.getEBGraph();
      const EBFaceFAB& flux       = a_flux[dit()][dir];
      const IntVectSet ivs        = grown_set[dit()] & box;
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;

      for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	interpol(face, comp) = flux(face, comp);
      }
    }

    (*(m_interpolant[dir])[a_lvl]).exchange();
  }
}

void cdr_solver::set_plot_variables(){
  CH_TIME("cdr_solver::set_plot_variables");
  if(m_verbosity > 5){
    pout() << m_name + "::set_plot_variables" << endl;
  }

  m_plot_phi = false;
  m_plot_vel = false;
  m_plot_dco = false;
  m_plot_src = false;

  ParmParse pp("cdr_solver");
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi") m_plot_phi = true;
    else if(str[i] == "vel") m_plot_vel = true;
    else if(str[i] == "dco") m_plot_dco = true; 
    else if(str[i] == "src") m_plot_src = true;
  }
}

void cdr_solver::write_plot_file(){
  CH_TIME("cdr_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }


  // Number of output components and their names
  const int ncomps = get_num_plotvars();
  const Vector<std::string> names = get_plotvar_names();

  // Allocate storage
  EBAMRCellData output;
  m_amr->allocate(output, m_phase, ncomps);
  data_ops::set_value(output, 0.0);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0;
  write_plot_data(output, icomp);

  // Filename
  char file_char[100];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);

  // Alias
  Vector<LevelData<EBCellFAB>* > output_ptr;
  m_amr->alias(output_ptr, output);
  
  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->get_grids(),
	      output_ptr,
	      names,
	      m_amr->get_domains()[0].domainBox(),
	      m_amr->get_dx()[0],
	      m_dt,
	      m_time,
	      m_amr->get_ref_rat(),
	      m_amr->get_finest_level() + 1,
	      true,
	      covered_values,
	      IntVect::Unit);
}

void cdr_solver::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("cdr_solver::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  // Plot state
  if(m_plot_phi) {
    write_data(a_output, a_comp, m_state, true);
  }

  // Plot diffusion coefficients
  if(m_plot_dco && m_diffusive) { // Need to compute the cell-centerd stuff first
    data_ops::set_value(m_scratch, 0.0);
    data_ops::average_face_to_cell(m_scratch, m_diffco, m_amr->get_domains());
    write_data(a_output, a_comp, m_scratch,   false);
  }

  // Plot source terms
  if(m_plot_src) {
    write_data(a_output, a_comp, m_source,    false);
  }

  // Plot velocities
  if(m_plot_vel && m_mobile) {
    write_data(a_output, a_comp, m_velo_cell, false);
  }

  // Plot EB fluxes
  if(m_plot_ebf && m_mobile){
    data_ops::set_value(m_scratch, 0.0);
    data_ops::incr(m_scratch, m_ebflux, 1.0);
    write_data(a_output, a_comp, m_scratch, false);
  }
}

void cdr_solver::write_data(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp){
  CH_TIME("cdr_solver::write_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_data" << endl;
  }

  const int comp = 0;
  const int ncomp = a_data[0]->nComp();

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_phase, ncomp);
  data_ops::copy(scratch, a_data);

  if(a_interp){
    m_amr->interpolate_to_centroids(scratch, phase::gas);
  }

  m_amr->average_down(scratch, m_phase);
  m_amr->interp_ghost(scratch, m_phase);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }

  data_ops::set_covered_value(a_output, a_comp, 0.0);

  a_comp += ncomp;
}

void cdr_solver::write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("cdr_solver::write_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::write_checkpoint_level" << endl;
  }

  // Write state vector
  write(a_handle, *m_state[a_level], m_name);
  write(a_handle, *m_source[a_level],   m_name+"_src");
}

void cdr_solver::read_checkpoint_level(HDF5Handle& a_handle, const int a_level){
  CH_TIME("cdr_solver::read_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::read_checkpoint_level" << endl;
  }

  read<EBCellFAB>(a_handle, *m_state[a_level], m_name, m_amr->get_grids()[a_level], Interval(0,0), false);
  read<EBCellFAB>(a_handle, *m_source[a_level], m_name+"_src", m_amr->get_grids()[a_level], Interval(0,0), false);
}

Real cdr_solver::compute_cfl_dt(){
  CH_TIME("cdr_solver::compute_cfl_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_cfl_dt" << endl;
  }

#if 0
  Real min_dt = 1.E99;
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real dx = m_amr->get_dx()[lvl];
    Real maxVal = 0.0;

    for (int dir = 0; dir < SpaceDim; dir++){
      Real maxloc = 0.0;
      Real minloc = 0.0;
      EBLevelDataOps::getMaxMin(maxloc, minloc, *m_velo_cell[lvl], dir, true);
      maxVal = std::max(maxVal, maxloc);
    }
    if(maxVal > 0.0){
      min_dt = std::min(min_dt, dx/maxVal);
    }
  }

  return min_dt;
#else

  Real min_dt = 1.E99;

  if(m_mobile){
    const Real SAFETY = 1.E-10; 

    const int comp  = 0;
    const int ncomp = 1;
    const int finest_level = m_amr->get_finest_level();
    const Interval interv(comp, comp);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

      Real max_vel = 0.0;
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& velo  = (*m_velo_cell[lvl])[dit()];
	const Box box          = dbl.get(dit());

#if 0 // I think this is wrong, the CFL should be dx/(|vx| + |vy| + |vz|)
	int fab_type;
	const BaseFab<char>& mask     = ebgraph.getMask(fab_type);
	const BaseFab<Real>& velo_fab = velo.getSingleValuedFAB();

	if(fab_type != -1){// not all covered
	  if(fab_type == 0){ // Has irregular cells
	    for (int dir = 0; dir < SpaceDim; dir++){
	      FORT_GET_MAX_VEL(CHF_REAL(max_vel),
			       CHF_CONST_FRA1(velo_fab, dir),
			       CHF_BOX(box),
			       CHF_CONST_FBA1(mask, 0));
	    }
	    for (VoFIterator vofit(ebisbox.getMultiCells(box), ebgraph); vofit.ok(); ++vofit){
	      const VolIndex& vof = vofit();
	      for (int dir = 0; dir < SpaceDim; dir++){
		max_vel = Max(max_vel, Abs(velo(vof, dir)));
	      }
	    }
	  }
	  else { // all regular
	    for (int dir = 0; dir < SpaceDim; dir++){
	      FORT_GET_MAXNORM(CHF_REAL(max_vel),
			       CHF_CONST_FRA1(velo_fab, dir),
			       CHF_BOX(box));
	    }
	  }
	}
#else // This is the correct way of computing this (I think)
	const BaseFab<Real>& velo_fab = velo.getSingleValuedFAB();
	FORT_ADVECTIVE_CFL_DT(CHF_CONST_FRA(velo_fab),
			      CHF_CONST_REAL(dx),
			      CHF_BOX(box),
			      CHF_REAL(min_dt));

	VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	
	  Real vel = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    vel += Abs(velo(vof, dir));
	  }
	  const Real thisdt = dx/vel;
	  min_dt = Min(thisdt, min_dt);
	}
#endif
      }

    
      max_vel = Max(max_vel, SAFETY);
      min_dt  = Min(min_dt, dx/max_vel);
    }

  
#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){
      MayDay::Error("cdr_solver::compute_cfl_dt() - communication error on norm");
    }
    min_dt = tmp;
#endif
  }


  return min_dt;
#endif
}

Real cdr_solver::compute_diffusive_dt(){
  CH_TIME("cdr_solver::compute_diffusive_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusive_dt" << endl;
  }

  Real min_dt = 1.E99;

  if(m_diffusive){
    const int comp  = 0;
    const int ncomp = 1;
    const int finest_level = m_amr->get_finest_level();
    const Interval interv(comp, comp);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& velo     = (*m_velo_cell[lvl])[dit()];
	const Box box             = dbl.get(dit());

	// Regular faces
	for (int dir = 0; dir < SpaceDim; dir++){
	  const EBFaceFAB& diffco = (*m_diffco[lvl])[dit()][dir];

	  const BaseFab<Real>& diffco_fab = diffco.getSingleValuedFAB();
	  Box facebox = box;
	  facebox.surroundingNodes(dir);

	  FORT_DIFFUSIVE_DT(CHF_CONST_FRA1(diffco_fab, comp),
			    CHF_CONST_REAL(dx),
			    CHF_BOX(facebox),
			    CHF_REAL(min_dt));
	}

	// Irregular faces
	const BaseIVFAB<Real>& diffco = (*m_diffco_eb[lvl])[dit()];

	VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();
	  const Real thisdt = dx*dx/(2*SpaceDim*diffco(vof, comp));

	  min_dt = Min(min_dt, thisdt);
	}
      }
    }
  
#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){
      MayDay::Error("cdr_solver::compute_diffusive_dt() - communication error on norm");
    }
    min_dt = tmp;
#endif
  }

  return min_dt;
}

Real cdr_solver::compute_source_dt(const Real a_max, const Real a_tolerance){
  CH_TIME("cdr_solver::compute_source_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_source_dt" << endl;
  }

  const int comp = 0;
  
  Real min_dt = 1.E99;

  if(a_max > 0.0){
    const int finest_level = m_amr->get_finest_level();
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& state  = (*m_state[lvl])[dit()];
	const EBCellFAB& source = (*m_source[lvl])[dit()];
	const Box box           = dbl.get(dit());

	const BaseFab<Real>& state_fab  = state.getSingleValuedFAB();
	const BaseFab<Real>& source_fab = source.getSingleValuedFAB();


	FORT_SOURCE_DT(CHF_REAL(min_dt),
		       CHF_CONST_FRA1(state_fab, comp),
		       CHF_CONST_FRA1(source_fab, comp),
		       CHF_CONST_REAL(a_tolerance),
		       CHF_CONST_REAL(a_max),
		       CHF_BOX(box));

	VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  const Real phi = state(vof, comp);
	  const Real src = source(vof, comp);

	  Real thisdt = 1.E99;
	  if(Abs(phi) > a_tolerance*a_max && src > 0.0){
	    thisdt = Abs(phi/src);
	  }
	  min_dt = Min(min_dt, thisdt);
	}
      }
    }
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("cdr_solver::compute_source_dt() - communication error on norm");
  }
  min_dt = tmp;
#endif
  
  return min_dt;

}

Real cdr_solver::compute_mass(){
  CH_TIME("cdr_solver::compute_mass");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_mass" << endl;
  }

  Real mass = 0.;
  const int base = 0;
  const Real dx = m_amr->get_dx()[base];
  m_amr->average_down(m_state, m_phase);
  
  data_ops::kappa_sum(mass, *m_state[base]);
  mass *= pow(dx, SpaceDim);
  
  return mass;
}

Real cdr_solver::compute_charge(){
  CH_TIME("cdr_solver::compute_charge");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_charge" << endl;
  }

  const Real Q = this->compute_mass()*m_species->get_charge();

  return Q;
}

bool cdr_solver::extrap_source() const {
  return m_extrap_source;
}

bool cdr_solver::is_diffusive(){
  return m_diffusive;
}

bool cdr_solver::is_mobile(){
  return m_mobile;
}

EBAMRCellData& cdr_solver::get_state(){
  return m_state;
}

EBAMRCellData& cdr_solver::get_source(){
  return m_source;
}

EBAMRCellData& cdr_solver::get_velo_cell(){
  return m_velo_cell;
}

EBAMRFluxData& cdr_solver::get_velo_face(){
  return m_velo_face;
}

EBAMRIVData& cdr_solver::get_velo_eb(){
  return m_velo_eb;
}

EBAMRFluxData& cdr_solver::get_diffco_face(){
  return m_diffco;
}

EBAMRIVData& cdr_solver::get_diffco_eb(){
  return m_diffco_eb;
}

EBAMRIVData& cdr_solver::get_ebflux(){
  return m_ebflux;
}

EBAMRIFData& cdr_solver::get_domainflux(){
  return m_domainflux;
}

void cdr_solver::parse_domain_bc(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("domain_bc", str);
  if(str == "kinetic"){
    set_domain_bc(cdr_bc::external);
  }
  else if(str == "outflow"){
    set_domain_bc(cdr_bc::outflow);
  }
  else if(str == "wall"){
    set_domain_bc(cdr_bc::wall);
  }
  else if(str == "extrap"){
    set_domain_bc(cdr_bc::extrap);
  }
  else{
    MayDay::Abort("cdr_solver::parse_domain_bc - unknown BC requested");
  }
}

void cdr_solver::parse_extrap_source(){
  CH_TIME("cdr_solver::parse_extrap_source");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_extrap_source" << endl;
  }
  
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("extrap_source", str);
  m_extrap_source = (str == "true") ? true : false;
}

void cdr_solver::parse_conservation(){
  CH_TIME("cdr_solver::parse_conservation");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_conservation" << endl;
  }

  ParmParse pp(m_class_name.c_str());

  pp.get("redist_mass_weighted", m_redist_mass_weighted);
  pp.get("blend_conservation",   m_blend_conservation);
}

void cdr_solver::parse_plot_vars(){
  ParmParse pp(m_class_name.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  m_plot_phi = false;
  m_plot_vel = false;
  m_plot_dco = false;
  m_plot_src = false;
  m_plot_ebf = false;
  
  for (int i = 0; i < num; i++){
    if(     str[i] == "phi")    m_plot_phi = true;
    else if(str[i] == "vel")    m_plot_vel = true;
    else if(str[i] == "dco")    m_plot_dco = true; 
    else if(str[i] == "src")    m_plot_src = true;
    else if(str[i] == "ebflux") m_plot_ebf = true;
  }
}

void cdr_solver::define_interpolant(){
  CH_TIME("cdr_solver::define_interpolant");
  if(m_verbosity > 5){
    pout() << m_name + "::define_interpolant" << endl;
  }

  const int ncomp = 1;


  for (int dir = 0; dir < SpaceDim; dir++){
    //    Vector<RefCountedPtr<LevelData<BaseIFFAB<Real> > > > interpolant = m_interpolant[SpaceDim];
    //    Vector<RefCountedPtr<LayoutData<BaseIFFAB<Real> > > > grown_sets = m_interp_sets[SpaceDim];

    m_interpolant[dir].resize(1 + m_amr->get_finest_level());
    m_interp_sets[dir].resize(1 + m_amr->get_finest_level());
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const ProblemDomain& domain  = m_amr->get_domains()[lvl];
      
      (m_interpolant[dir])[lvl] = RefCountedPtr<LevelData<BaseIFFAB<Real> > > (new LevelData<BaseIFFAB<Real> >());
      (m_interp_sets[dir])[lvl] = RefCountedPtr<LayoutData<IntVectSet> > (new LayoutData<IntVectSet>());

      EBArith::defineFluxInterpolant(*(m_interpolant[dir])[lvl],
				     *(m_interp_sets[dir])[lvl], 
				     dbl,
				     ebisl,
				     domain,
				     ncomp,
				     dir);
    }
  }
}

void cdr_solver::make_non_negative(EBAMRCellData& a_phi){
  CH_TIME("cdr_solver::make_non_negative");
  if(m_verbosity > 5){
    pout() << m_name + "::make_non_negative" << endl;
  }

  // TLDR: This function injects mass into cut-cells where the state is negative, and removes mass from neighboring
  //       regular cells. This is essentially like flooring, with the exception that total mass is conserved. Very useful
  //       stuff for cut-cells where non-negativity can not be guaranteed. 

  const int comp  = 0;
  const int ncomp = 1;

  data_ops::set_value(m_mass_diff, 0.0);
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    for (DataIterator dit    = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& state       = (*a_phi[lvl])[dit()];
      const EBISBox& ebisbox = state.getEBISBox();
      const Box box          = dbl.get(dit());


      VoFIterator& vofit = (*m_amr->get_vofit(m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa = ebisbox.volFrac(vof);
	const Real phi = state(vof, comp);
	if(phi < 0.0){
	  state(vof,comp) = 0.0;
	  (*m_mass_diff[lvl])[dit()](vof,comp) = kappa*phi;
	}
      }
    }
  }

  this->increment_concentration_redist(m_mass_diff);
  this->concentration_redistribution(a_phi, m_mass_diff);
}

void cdr_solver::GWN_diffusion_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_solver::GWN_diffusion_source");
  if(m_verbosity > 5){
    pout() << m_name + "::GWN_diffusion_source" << endl;
  }

  if(m_diffusive){
  const int comp  = 0;
  const int ncomp = 1;

  this->fill_GWN(m_scratchFluxTwo, 1.0);                         // Gaussian White Noise = W/sqrt(dV)
  this->smooth_heaviside_faces(m_scratchFluxOne, a_cell_states); // m_scratchFluxOne = phis
  data_ops::multiply(m_scratchFluxOne, m_diffco);                // m_scratchFluxOne = D*phis
  data_ops::scale(m_scratchFluxOne, 2.0);                        // m_scratchFluxOne = 2*D*phis
  data_ops::square_root(m_scratchFluxOne);                       // m_scratchFluxOne = sqrt(2*D*phis)

#if 0 // Debug
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    for (int dir = 0; dir <SpaceDim; dir++){
      Real max, min;
      EBLevelDataOps::getMaxMin(max, min, *m_scratchFluxOne[lvl], 0, dir);
      if(min < 0.0 || max < 0.0){
	MayDay::Abort("cdr_solver::GWN_diffusion_source - negative face value");
      }
    }
  }
#endif
  
  data_ops::multiply(m_scratchFluxOne, m_scratchFluxTwo);                     // Holds random, cell-centered flux
  this->conservative_divergence(a_ransource, m_scratchFluxOne, m_eb_zero); 
  
#if 0 // Debug
  m_amr->average_down(a_ransource, m_phase);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(EBLevelDataOps::checkNANINF(*a_ransource[lvl])){
      MayDay::Abort("cdr_solver::GWN_diffusion_source - something is wrong");
    }
  }
#endif
  }
  else{
    data_ops::set_value(a_ransource, 0.0);
  }
}

void cdr_solver::smooth_heaviside_faces(EBAMRFluxData& a_face_states, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_solver::smooth_heaviside_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::smooth_heaviside_faces" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& face_states       = (*a_face_states[lvl])[dit()][dir];
	const EBCellFAB& cell_states = (*a_cell_states[lvl])[dit()];

	BaseFab<Real>& reg_face       = face_states.getSingleValuedFAB();
	const BaseFab<Real>& reg_cell = cell_states.getSingleValuedFAB();

	Box facebox = box;
	//	facebox.grow(dir,1);
	facebox.surroundingNodes(dir);

	// This will also do irregular cells and boundary faces
	FORT_HEAVISIDE_MEAN(CHF_FRA1(reg_face, comp),  
			    CHF_CONST_FRA1(reg_cell, comp),
			    CHF_CONST_INT(dir),
			    CHF_CONST_REAL(dx),
			    CHF_BOX(facebox));

	// Fix up irregular cell faces
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irreg, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();


	  const VolIndex lovof = face.getVoF(Side::Lo);
	  const VolIndex hivof = face.getVoF(Side::Hi);

	  const Real loval = Max(0.0, cell_states(lovof, comp));
	  const Real hival = Max(0.0, cell_states(hivof, comp));


	  Real Hlo;
	  if(loval*vol <= 0.0){
	    Hlo = 0.0;
	  }
	  else if(loval*vol >= 1.0){
	    Hlo = 1.0;
	  }
	  else{
	    Hlo = loval*vol;
	  }

	  Real Hhi;
	  if(hival*vol <= 0.0){
	    Hhi = 0.0;
	  }
	  else if(hival*vol >= 1.0){
	    Hhi = 1.0;
	  }
	  else{
	    Hhi = hival*vol;
	  }


	  face_states(face, comp) = 0.5*(hival + loval)*Hlo*Hhi;
	}

	// No random flux on domain faces. Reset those. 
	for (SideIterator sit; sit.ok(); ++sit){
	  Box sidebox;
	  if(sit() == Side::Lo){
	    sidebox = bdryLo(domain, dir, 1);
	  }
	  else if(sit() == Side::Hi){
	    sidebox = bdryHi(domain, dir, 1);
	  }
	  
	  sidebox &= facebox;

	  Box cellbox = sidebox.enclosedCells(dir);

	  const IntVectSet ivs(cellbox);
	  const FaceStop::WhichFaces stopcrit = FaceStop::AllBoundaryOnly;
	  for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	    face_states(faceit(), comp) = 0.0;
	  }
	}
      }
    }

    // Covered is bogus
    EBLevelDataOps::setCoveredVal(*a_face_states[lvl], 0.0);
  }
}

void cdr_solver::fill_GWN(EBAMRFluxData& a_noise, const Real a_sigma){
  CH_TIME("cdr_solver::fill_GWN");
  if(m_verbosity > 5){
    pout() << m_name + "::fill_GWN" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  std::normal_distribution<double> GWN(0.0, a_sigma);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    const Real ivol              = sqrt(1./vol);
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& noise = (*a_noise[lvl])[dit()][dir];

	Box facebox = box;
	facebox.surroundingNodes(dir);

	noise.setVal(0.0);

	// Regular faces
	BaseFab<Real>& noise_reg = noise.getSingleValuedFAB();
	for (BoxIterator bit(facebox); bit.ok(); ++bit){
	  const IntVect iv = bit();
	  noise_reg(iv, comp) = GWN(m_rng)*ivol;
	}

	// // Irregular faces
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irreg, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  noise(face, comp) = GWN(m_rng)*ivol;
	}
      }
    }
  }
}

void cdr_solver::parse_rng_seed(){
  CH_TIME("cdr_solver::parse_rng_seed");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng_seed" << endl;
  }  
  ParmParse pp(m_class_name.c_str());
  pp.get("seed", m_seed);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();

  m_rng = std::mt19937_64(m_seed);
}

void cdr_solver::parse_plotmode(){
  CH_TIME("cdr_solver::parse_plotmode");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_plotmode" << endl;
  }  
  ParmParse pp(m_class_name.c_str());

  m_plot_numbers = false;
  
  std::string str;
  pp.get("plot_mode", str);
  if(str == "density"){
    m_plot_numbers = false;
  }
  else if(str == "numbers"){
    m_plot_numbers = true;
  }
}
