/*!
  @file cdr_solver.cpp
  @brief Implementation of cdr_solver.H
  @author Robert Marskar
  @date Nov. 2017
  @todo See if we can make EBAdvectLevelIntegrator a static in cdr_gdnv
*/

#include "cdr_solver.H"
#include "data_ops.H"

#include <EBArith.H>
#include <EBAMRIO.H>

#if 1 // Move to cdr_gdnv
#include <ExtrapAdvectBC.H>
#endif

cdr_solver::cdr_solver(){

  this->set_phase(phase::gas);
  this->set_time(0, 0., 0.);
  
  m_name = "cdr_solver";
}

cdr_solver::~cdr_solver(){

}

int cdr_solver::query_ghost() const {
  CH_TIME("cdr_solver::query_ghost");
  if(m_verbosity > 5){
    pout() << m_name + "::query_ghost" << endl;
  }

  return 3;
}

void cdr_solver::sanity_check(){
  CH_TIME("cdr_solver::sanity_check");
  if(m_verbosity > 5){
    pout() << m_name + "::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_physdom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_species.isNull());
  CH_assert(!m_ebis.isNull());
}

void cdr_solver::set_species(const RefCountedPtr<species> a_species){
  CH_TIME("cdr_solver::set_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_species" << endl;
  }

  m_species   = a_species;
  m_name      = m_species->get_name();
  m_diffusive = m_species->is_diffusive();
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

void cdr_solver::set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("cdr_solver::set_ebis");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebis" << endl;
  }

  m_ebis = a_ebis;
}

void cdr_solver::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("cdr_solver::set_physical_domain");
  if(m_verbosity > 5){
    pout() << m_name + "::set_physical_domain" << endl;
  }

  m_physdom = a_physdom;
}

void cdr_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("cdr_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;
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
    a_velo[lvl]->copyTo(*m_velo_cell[lvl]);
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
  }

  m_amr->average_down(m_velo_cell, m_phase);
  m_amr->interp_ghost(m_velo_cell, m_phase);
}

void cdr_solver::set_diffco(const EBAMRFluxData& a_diffco, const EBAMRIVData& a_diffco_eb){
  CH_TIME("cdr_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_diffco[lvl]->copyTo(*m_diffco[lvl]);
    a_diffco_eb[lvl]->copyTo(*m_diffco_eb[lvl]);
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
  }

  m_amr->average_down(m_diffco, m_phase);
  m_amr->average_down(m_diffco_eb, m_phase);
}

void cdr_solver::set_source(const EBAMRCellData& a_source){
  CH_TIME("cdr_solver::set_source");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->copyTo(*m_source[lvl]);
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

void cdr_solver::set_phase(const phase::which_phase a_phase){
  CH_TIME("cdr_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void cdr_solver::set_verbosity(const int a_verbosity){
  CH_TIME("cdr_solver::set_verbosity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
  
  m_verbosity = a_verbosity;
}

void cdr_solver::initial_data(){
  CH_TIME("cdr_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }

  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (DataIterator dit = m_state[lvl]->dataIterator(); dit.ok(); ++dit){
      EBCellFAB& state       = (*m_state[lvl])[dit()];
      const Box box          = m_state[lvl]->disjointBoxLayout().get(dit());
      const EBISBox& ebisbox = state.getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, origin, m_amr->get_dx()[lvl]*RealVect::Unit);
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_species->initial_data(pos, m_time);
	}
      }
    }
  }
}

void cdr_solver::allocate_internals(){
  CH_TIME("cdr_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  const int sca = 1;
  const int vec = SpaceDim;

  m_amr->allocate(m_state,      m_phase, sca);
  m_amr->allocate(m_source,     m_phase, sca);
  m_amr->allocate(m_velo_face,  m_phase, sca);
  m_amr->allocate(m_velo_cell,  m_phase, vec); 
  m_amr->allocate(m_ebflux,     m_phase, sca);
  m_amr->allocate(m_diffco,     m_phase, sca);
  m_amr->allocate(m_diffco_eb,  m_phase, sca);

  data_ops::set_value(m_state,      0.0);
  data_ops::set_value(m_source,     0.0);
  data_ops::set_value(m_velo_face,  0.0);
  data_ops::set_value(m_velo_cell,  0.0);
  data_ops::set_value(m_ebflux,     0.0);
  data_ops::set_value(m_diffco,     0.0);
  data_ops::set_value(m_diffco_eb,  0.0);
}

void cdr_solver::advance(const Real& a_dt){
  CH_TIME("cdr_solver::advance(dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(dt)" << endl;
  }
  
  this->advance(m_state, a_dt);
}

void cdr_solver::advance(EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_solver::state(dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::state(dt)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  
  // Compute right-hand-side
  EBAMRCellData rhs;
  m_amr->allocate(rhs, m_phase, ncomp);
  this->compute_rhs(a_state, rhs, a_dt);

  // Increment. This is the explicit Euler method
  data_ops::incr(a_state, rhs, a_dt);

  m_amr->average_down(a_state, m_phase);
  m_amr->interp_ghost(a_state, m_phase);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt){
  CH_TIME("cdr_solver::compute_rhs(rhs, dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rhs(rhs, dt)" << endl;
  }
  
  this->compute_rhs(a_rhs, m_state, a_dt);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_solver::compute_rhs(rhs, state, dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rhs(rhs, state, dt)" << endl;
  }

  data_ops::set_value(a_rhs, 0.0);  

  const int comp  = 0;
  const int ncomp = 1;

  // Advective derivative
  EBAMRCellData advective_term;
  EBAMRCellData diffusion_term;

  m_amr->allocate(advective_term, m_phase, ncomp);
  if(this->is_diffusive()){
    m_amr->allocate(diffusion_term, m_phase, ncomp);
  }

  // Advective term
  this->compute_advective_derivative(advective_term, a_state);
  data_ops::incr(a_rhs, advective_term, -1.0);

  // Diffusion term
  if(this->is_diffusive()){
    this->compute_diffusion_term(diffusion_term, a_state);
    data_ops::incr(a_rhs, diffusion_term, 1.0);
  }

  // Source term
  data_ops::incr(a_rhs, m_source, 1.0);

  m_amr->average_down(a_rhs, m_phase);
  m_amr->interp_ghost(a_rhs, m_phase);
}

void cdr_solver::compute_advective_derivative(EBAMRCellData& a_divF, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::compute_advective_derivative(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_advective_derivative(divF, state)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  EBAMRFluxData face_state;
  EBAMRIVData   div_nc;
  EBAMRIVData   mass_diff;

  m_amr->allocate(face_state, m_phase, ncomp);
  m_amr->allocate(div_nc,     m_phase, ncomp);
  m_amr->allocate(mass_diff,  m_phase, ncomp);

  data_ops::set_value(a_divF,     0.0); 
  data_ops::set_value(face_state, 0.0);
  data_ops::set_value(div_nc,     0.0);
  data_ops::set_value(mass_diff,  0.0);

  // Compute the advective derivative
  this->average_velo_to_faces(m_velo_face, m_velo_cell);            // Average cell-centered velocities to face centers
  this->extrapolate_to_faces(face_state, a_state);                  // Face extrapolation to cell-centered faces
  this->interpolate_to_centroids(face_state);                       // Interpolate to centroids
  this->conservative_divergence(a_divF, face_state, m_velo_face);   // a_divF holds the conservative divergence
  this->nonconservative_divergence(div_nc, a_divF);                 // Compute non-conservative divergence
  this->hybrid_divergence(a_divF, mass_diff, div_nc);               // Make divF = hybrid divergence. Compute mass diff.
  this->increment_flux_register(face_state, m_velo_face);           // Increment flux registers
  this->hyperbolic_redistribution(a_divF, mass_diff, a_state);      // Redistribute mass into hybrid divergence

  const bool ebcf = m_amr->get_ebcf();
  if(ebcf){ 
    this->coarse_fine_increment(mass_diff);     // Increment the coarse-fine redistribution objects
    this->increment_redist_flux();              // Increment flux registers with the redistribution stuff
    this->coarse_fine_redistribution(a_divF);   // Redistribute
  }

  this->reflux(a_divF); // Reflux
}

void cdr_solver::compute_diffusion_term(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::compute_diffusion_term");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_term" << endl;
  }

  MayDay::Abort("cdr_solver::compute_diffusion_term - not implemented");
}

void cdr_solver::average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell){
  CH_TIME("cdr_solver::average_velo_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::average_cell_to_face(*a_velo_face[lvl], *a_velo_cell[lvl], m_amr->get_domains()[lvl]);
  }
}

#if 1 // This should be moved to cdr_gdnv
void cdr_solver::extrapolate_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::extrapolate_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::extrapolate_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  const int comp = 0;
  
  EBAdvectPatchIntegrator::setCurComp(0);
  EBAdvectPatchIntegrator::setDoingVel(0);

  RefCountedPtr<ExtrapAdvectBCFactory> bcfact = RefCountedPtr<ExtrapAdvectBCFactory>
    (new ExtrapAdvectBCFactory());


  for (int lvl = 0; lvl <= finest_level; lvl++){
    const RefCountedPtr<EBAdvectLevelIntegrator>& leveladvect = m_amr->get_level_advect(m_phase)[lvl];
    leveladvect->resetBCs(bcfact);
    CH_assert(!leveladvect.isNull());
    leveladvect->advectToFacesBCG(*a_face_state[lvl],
				  *a_state[lvl],
				  *m_velo_cell[lvl],
				  *m_velo_face[lvl],
				  a_state[lvl],// should be coarse
				  a_state[lvl],// should be coars
				  m_velo_cell[lvl],// should be coarse
				  m_velo_cell[lvl],// should be coarse
				  m_time,
				  m_time,
				  m_time,
				  m_dt,
				  m_source[lvl], 
				  m_source[lvl],// should be coarsse
				  m_source[lvl]);// should be coarsse
  }
						

  MayDay::Abort("cdr_solver::extrapolate_to_faces - not implemented");
}
#endif

void cdr_solver::interpolate_to_centroids(EBAMRFluxData& a_face_state){
  CH_TIME("cdr_solver::interpolate_to_centroids");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_to_centroids" << endl;
  }

  MayDay::Abort("cdr_solver::interpolate_to_centroids - not implemented");
}

void cdr_solver::conservative_divergence(EBAMRCellData&       a_cons_div,
					 const EBAMRFluxData& a_face_vel,
					 const EBAMRFluxData& a_face_state){
  CH_TIME("cdr_solver::conservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::conservative_divergence" << endl;
  }

  MayDay::Abort("cdr_solver::conservative_divergence - not implemented");
}

void cdr_solver::nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divF){
  CH_TIME("cdr_solver::nonconservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::nonconservative_divergence" << endl;
  }

  MayDay::Abort("cdr_solver::nonconservative_divergence - not implemented");
}

void cdr_solver::hybrid_divergence(EBAMRCellData&     a_hybrid_div,
				   EBAMRIVData&       a_mass_diff,
				   const EBAMRIVData& a_noncons_div){
  CH_TIME("cdr_solver::hybrid_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::hybrid_divergence" << endl;
  }

  MayDay::Abort("cdr_solver::hybrid_divergence - not implemented");
}

void cdr_solver::increment_flux_register(const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_velo_face){
  CH_TIME("cdr_solver::increment_flux_register");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_flux_register" << endl;
  }

  MayDay::Abort("cdr_solver::increment_flux_register - not implemented");
}


void cdr_solver::coarse_fine_increment(const EBAMRIVData& m_mass_diff){
  CH_TIME("cdr_solver::coarse_fine_increment");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_increment" << endl;
  }

  MayDay::Abort("cdr_solver::coarse_fine_increment - not implemented");
}

void cdr_solver::hyperbolic_redistribution(EBAMRCellData&       a_del_vel_rho,
					   const EBAMRIVData&   a_mass_diff,
					   const EBAMRCellData& a_redist_weights){
  CH_TIME("cdr_solver::hyberbolic_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::hyperbolic_redistribution" << endl;
  }

  MayDay::Abort("cdr_solver::hyperbolic_redistribution - not implemented");
}

void cdr_solver::increment_redist_flux(){
  CH_TIME("cdr_solver::increment_redist_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist_flux" << endl;
  }

  MayDay::Abort("cdr_solver::increment_redist_flux - not implemented");
}

void cdr_solver::coarse_fine_redistribution(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::coarse_fine_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_redistribution" << endl;
  }

  MayDay::Abort("cdr_solver::coarse_fine_redistribution - not implemented");
}

void cdr_solver::reflux(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::reflux");
  if(m_verbosity > 5){
    pout() << m_name + "::reflux" << endl;
  }

  MayDay::Abort("cdr_solver::reflux - not implemented");
}

#ifdef CH_USE_HDF5
void cdr_solver::write_plot_file(){
  CH_TIME("cdr_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);

  const int ncomps = 2 + SpaceDim;
  Vector<string> names(ncomps);
  names[0] = "density";
  names[1] = "source";
  names[2] = "x-velocity";
  names[3] = "y-velocity";
  if(SpaceDim == 3){
    names[4] = "y-velocity";
  }

  EBAMRCellData output;
  m_amr->allocate(output, m_phase, ncomps);

  for (int lvl = 0; lvl < output.size(); lvl++){
    LevelData<EBCellFAB>& state  = *m_state[lvl];
    LevelData<EBCellFAB>& source = *m_source[lvl];
    LevelData<EBCellFAB>& velo   = *m_velo_cell[lvl];

    state.copyTo(Interval(0,0),           *output[lvl],  Interval(0,0));
    source.copyTo(Interval(0,0),          *output[lvl],  Interval(1,1));
    velo.copyTo(Interval(0,SpaceDim - 1), *output[lvl],  Interval(2, 2 + (SpaceDim-1)));
  }


  // Transform to centroids
  irreg_amr_stencil<centroid_interp>& sten = m_amr->get_centroid_interp_stencils(phase::gas);
  sten.apply(output);

  //
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
	      false,
	      covered_values);
}
#endif

Real cdr_solver::compute_dt(){
  CH_TIME("cdr_solver::compute_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_dt" << endl;
  }

  MayDay::Abort("cdr_solver::compute_dt - not implemented");
}

bool cdr_solver::is_diffusive(){
  CH_TIME("cdr_solver::is_diffusive");
  if(m_verbosity > 5){
    pout() << m_name + "::is_diffusive" << endl;
  }
  
  return m_diffusive;
}

EBAMRCellData& cdr_solver::get_state(){
  return m_state;
}

EBAMRCellData& cdr_solver::get_source(){
  return m_source;
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
