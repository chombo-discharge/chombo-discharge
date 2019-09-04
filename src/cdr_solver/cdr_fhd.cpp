/*!
  @file   cdr_fhd.cpp
  @brief  Implementation of cdr_fhd.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_fhd.H"
#include "cdr_fhdF_F.H"
#include "data_ops.H"
#include <time.h>
#include <chrono>

#include <ParmParse.H>

cdr_fhd::cdr_fhd() : cdr_gdnv() {
  m_name       = "cdr_fhd";
  m_class_name = "cdr_fhd";
}

cdr_fhd::~cdr_fhd(){
  this->delete_covered();
}

void cdr_fhd::parse_options(){
  CH_TIME("cdr_fhd::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }
  
  parse_diffusion();    // Parses stochastic diffusion
  parse_rng_seed();     // Parses RNG seed
  parse_domain_bc();    // Parses domain BC options
  parse_mass_redist();  // Parses mass redistribution
  parse_hybrid_div();   // Parses options for hybrid divergence
  parse_slopelim();     // Parses slope limiter settings
  parse_plot_vars();    // Parses plot variables
  parse_gmg_settings(); // Parses solver parameters for geometric multigrid
  parse_plotmode();    // Parse plot mode
}

void cdr_fhd::parse_diffusion(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("stochastic_diffusion", str);
  
  m_stochastic_diffusion = (str == "true") ? true : false;
  m_stochastic_advection = false;
}

void cdr_fhd::parse_rng_seed(){
  ParmParse pp(m_class_name.c_str());
  pp.get("seed", m_seed);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();

  m_rng = new std::mt19937_64(m_seed);
}

void cdr_fhd::parse_plotmode(){
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

void cdr_fhd::advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_fhd::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler (no source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_euler(a_new_state, a_old_state, a_dt);
  }
  else{
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    EBAMRCellData ransource;
    m_amr->allocate(ransource, m_phase, 1);
    this->GWN_diffusion_source(ransource, a_old_state);


    // Do the regular aliasing stuff for passing into AMRMultiGrid
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    ransource);
      
    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);
      
    // Multigrid solver
    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
    MayDay::Abort("cdr_fhd::advance_euler - stochastic not implemented");
  }

  data_ops::floor(a_new_state, 0.0);
}

void cdr_fhd::advance_euler(EBAMRCellData&       a_new_state,
			    const EBAMRCellData& a_old_state,
			    const EBAMRCellData& a_source,
			    const Real           a_dt){
  CH_TIME("cdr_fhd::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler (with source)" << endl;
  }

  
  if(m_diffusive){
    if(!m_stochastic_diffusion){
      cdr_tga::advance_euler(a_new_state, a_old_state, a_source, a_dt);
    }
    else{
      bool converged = false;

      const int comp         = 0;
      const int ncomp        = 1;
      const int finest_level = m_amr->get_finest_level();

      EBAMRCellData ransource;
      m_amr->allocate(ransource, m_phase, 1);
      this->GWN_diffusion_source(ransource, a_old_state);

      // Make new state = old_state + dt*random_source
#if 0
      data_ops::copy(a_new_state, a_old_state);
      data_ops::incr(a_new_state, ransource, a_dt);
#else
      data_ops::incr(ransource, a_source, 1.0);
#endif
      
      // Do the regular aliasing stuff for passing into AMRMultiGrid
      Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
      m_amr->alias(new_state, a_new_state);
      m_amr->alias(old_state, a_old_state);
      m_amr->alias(source,    ransource);
      
      const Real alpha = 0.0;
      const Real beta  = 1.0;

      data_ops::set_value(m_diffco_eb, 0.0);
      
      // Multigrid solver
      m_eulersolver->resetAlphaAndBeta(alpha, beta);
      m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
      const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
      if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
	converged = true;
      }
    }
  }
  else{
    data_ops::copy(a_new_state, a_old_state);
  }

  data_ops::floor(a_new_state, 0.0);
}

void cdr_fhd::advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_fhd::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga (no source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_tga(a_new_state, a_old_state, a_dt);
  }
  else{
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    EBAMRCellData ransource;
    m_amr->allocate(ransource, m_phase, 1);
    this->GWN_diffusion_source(ransource, a_old_state);
      
    // Do the regular aliasing stuff for passing into AMRMultiGrid
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    ransource);
      
    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);
      
    // Multigrid solver
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
}

void cdr_fhd::advance_tga(EBAMRCellData&       a_new_state,
			  const EBAMRCellData& a_old_state,
			  const EBAMRCellData& a_source,
			  const Real           a_dt){
  CH_TIME("cdr_fhd::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga (with source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_tga(a_new_state, a_old_state, a_source, a_dt);
  }
  else{
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    EBAMRCellData ransource;
    m_amr->allocate(ransource, m_phase, 1);
    this->GWN_diffusion_source(ransource, a_old_state);
    data_ops::incr(ransource, a_source, 1.0);
      
    // Do the regular aliasing stuff for passing into AMRMultiGrid
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    ransource);
      
    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);
      
    // Multigrid solver
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
}

void cdr_fhd::GWN_diffusion_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_fhd::GWN_diffusion_source");
  if(m_verbosity > 5){
    pout() << m_name + "::GWN_diffusion_source" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  
  EBAMRFluxData ranflux;
  EBAMRFluxData GWN;

  // Putting this in here because it is really, really important that we don't send in any stupid ghost cells or negative
  // values for the diffusion advance
  EBAMRCellData& states = const_cast<EBAMRCellData&> (a_cell_states);
  m_amr->average_down(states, m_phase);
  m_amr->interp_ghost(states, m_phase);
  data_ops::floor(states, 0.0);
  
  m_amr->allocate(ranflux,   m_phase, ncomp);
  m_amr->allocate(GWN,       m_phase, ncomp);
  
  data_ops::set_value(a_ransource, 0.0);
  data_ops::set_value(ranflux,     0.0);
  data_ops::set_value(GWN,         0.0);

  this->fill_GWN(GWN, 1.0);                             // Gaussian White Noise
  this->smooth_heaviside_faces(ranflux, a_cell_states); // ranflux = phis
   data_ops::multiply(ranflux, m_diffco);                // ranflux = D*phis
   data_ops::scale(ranflux, 2.0);                        // ranflux = 2*D*phis
   data_ops::square_root(ranflux);                       // ranflux = sqrt(2*D*phis)

#if 1 // Debug
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    for (int dir = 0; dir <SpaceDim; dir++){
      Real max, min;
      EBLevelDataOps::getMaxMin(max, min, *ranflux[lvl], 0, dir);
      if(min < 0.0 || max < 0.0){
	MayDay::Abort("stop - negative face value");
      }
    }
  }
#endif
  data_ops::multiply(ranflux, GWN);                     // Holds random, cell-centered flux

  // Source term. 
  // I want to re-use conservative_divergence(), but that also computes with the EB fluxes. Back that up
  // first and use the already written and well-tested routine. Then copy back
  EBAMRIVData backup;
  m_amr->allocate(backup, m_phase, 1);
  data_ops::copy(backup, m_ebflux);
  data_ops::set_value(m_ebflux, 0.0);
  data_ops::set_value(a_ransource, 0.0);
  conservative_divergence(a_ransource, ranflux); // Compute the conservative divergence. This also refluxes. 
  data_ops::copy(m_ebflux, backup);

  //  data_ops::set_value(a_ransource, 0.0);

  
#if 1 // Debug
  m_amr->average_down(a_ransource, m_phase);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(EBLevelDataOps::checkNANINF(*a_ransource[lvl])){
      MayDay::Abort("something is wrong");
    }
  }
#endif
}

void cdr_fhd::smooth_heaviside_faces(EBAMRFluxData& a_face_states, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_fhd::smooth_heaviside_faces");
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

void cdr_fhd::fill_GWN(EBAMRFluxData& a_noise, const Real a_sigma){
  CH_TIME("cdr_fhd::fill_GWN");
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
	  noise_reg(iv, comp) = GWN(*m_rng)*ivol;
	}

	// Irregular faces
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irreg, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  noise(face, comp) = GWN(*m_rng)*ivol;
	}
      }
    }
  }
}

void cdr_fhd::compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist){
  CH_TIME("cdr_fhd::compute_divF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF(divF, state)" << endl;
  }

  if(m_mobile){
    cdr_tga::compute_divF(a_divF, a_state, a_extrap_dt, a_redist); // Deterministic flux
    if(m_stochastic_advection){ // Stochastic flux
      EBAMRCellData ransource;
      m_amr->allocate(ransource, m_phase, 1);
      this->GWN_advection_source(ransource, a_state);
      data_ops::incr(a_divF, ransource, 1.0);
    }
  }
  else{
    data_ops::set_value(a_divF, 0.0);
  }
}

void cdr_fhd::GWN_advection_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_fhd::GWN_advection_source");
  if(m_verbosity > 5){
    pout() << m_name + "::GWN_advection_source" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  
  EBAMRFluxData ranflux;
  EBAMRFluxData GWN;
  
  m_amr->allocate(ranflux,   m_phase, ncomp);
  m_amr->allocate(GWN,       m_phase, ncomp);
  
  data_ops::set_value(a_ransource, 0.0);
  data_ops::set_value(ranflux,   0.0);
  data_ops::set_value(GWN,       0.0);
      
  this->fill_GWN(GWN, 1.0);                             // Gaussian White Noise
  this->smooth_heaviside_faces(ranflux, a_cell_states); // ranflux = phis
  data_ops::square_root(ranflux);                       // ranflux = sqrt(phis)
  data_ops::multiply(ranflux, GWN);                     // ranflux = sqrt(phis)*W/sqrt(dV)
  data_ops::multiply(ranflux, m_velo_face);             // ranflux = v*sqrt(phis)*W/sqrt(dV)

  // Source term. 
  // I want to re-use conservative_divergence(), but that also computes with the EB fluxes. Back that up
  // first and use the already written and well-tested routine. Then copy back.
  EBAMRIVData backup;
  m_amr->allocate(backup, m_phase, 1);
  data_ops::copy(backup, m_ebflux);
  data_ops::set_value(m_ebflux, 0.0);
  conservative_divergence(a_ransource, ranflux); // Compute the conservative divergence. This also refluxes. 
  data_ops::copy(m_ebflux, backup);
}

void cdr_fhd::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("cdr_fhd::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  if(m_plot_phi) {
    if(!m_plot_numbers){ // Regular write
      write_data(a_output, a_comp, m_state, true);
    }
    else{ // Scale, write, and scale back
     
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_state[lvl], (pow(m_amr->get_dx()[lvl], 3)));
      }
      write_data(a_output, a_comp, m_state,    false);
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_state[lvl], 1./(pow(m_amr->get_dx()[lvl], 3)));
      }
    }
  }

  if(m_plot_dco && m_diffusive) { // Need to compute the cell-centerd stuff first
    data_ops::set_value(m_scratch, 0.0);
    data_ops::average_face_to_cell(m_scratch, m_diffco, m_amr->get_domains());
    write_data(a_output, a_comp, m_scratch,   false);
  }

  if(m_plot_src) {
    if(!m_plot_numbers){
      write_data(a_output, a_comp, m_source,    false);
    }
    else { // Scale, write, and scale back
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_source[lvl], (pow(m_amr->get_dx()[lvl], 3)));
      }
      write_data(a_output, a_comp, m_source,    false);
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_source[lvl], 1./(pow(m_amr->get_dx()[lvl], 3)));
      }
    }
  }

  if(m_plot_vel && m_mobile) {
    write_data(a_output, a_comp, m_velo_cell, false);
  }
}
