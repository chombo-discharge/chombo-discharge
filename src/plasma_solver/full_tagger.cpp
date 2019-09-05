/*!
  @file full_tagger.cpp
  @brief Implementation of full_tagger.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "full_tagger.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"

#include <EBArith.H>

full_tagger::full_tagger(){
  CH_TIME("full_tagger::full_tagger");
  if(m_verbosity > 5){
    pout() << "full_tagger::full_tagger" << endl;
  }
  
  m_name        = "full_tagger";
}

full_tagger::~full_tagger(){

}

void full_tagger::allocate_storage(){
  CH_TIME("full_tagger::allocate_storage");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_storage" << endl;
  }

  const int sca_ncomp = 1;
  const int vec_ncomp = SpaceDim;


  RefCountedPtr<cdr_layout> cdr = m_timestepper->get_cdr();
  RefCountedPtr<rte_layout> rte = m_timestepper->get_rte();

  m_amr->allocate(m_scratch,  m_phase, sca_ncomp);
  m_amr->allocate(m_E,        m_phase, vec_ncomp);
  m_amr->allocate(m_grad_E,   m_phase, vec_ncomp);
  m_amr->allocate(m_rho,      m_phase, sca_ncomp);
  m_amr->allocate(m_grad_rho, m_phase, vec_ncomp);

  m_cdr_densities.resize(m_plaskin->get_num_species());
  m_cdr_gradients.resize(m_plaskin->get_num_species());
  m_rte_densities.resize(m_plaskin->get_num_photons());

  for(cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_amr->allocate(m_cdr_densities[idx], m_phase, sca_ncomp);
    m_amr->allocate(m_cdr_gradients[idx], m_phase, vec_ncomp);
  }

  for(rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_amr->allocate(m_rte_densities[idx], m_phase, sca_ncomp);
  }
}

void full_tagger::deallocate_storage(){
  CH_TIME("full_tagger::deallocate_storage");
  if(m_verbosity > 5){
    pout() << m_name + "::deallocate_storage" << endl;
  }

  RefCountedPtr<cdr_layout> cdr = m_timestepper->get_cdr();
  RefCountedPtr<rte_layout> rte = m_timestepper->get_rte();

  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_E);
  m_amr->deallocate(m_grad_E);
  m_amr->deallocate(m_rho);
  m_amr->deallocate(m_grad_rho);

  for(cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_amr->deallocate(m_cdr_densities[idx]);
    m_amr->deallocate(m_cdr_gradients[idx]);
  }

  for(rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_amr->deallocate(m_rte_densities[idx]);
  }
}

void full_tagger::compute_tracers(){
  CH_TIME("full_tagger::compute_tracers");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_tracers" << endl;
  }

  this->allocate_storage();
  
  const RealVect origin = m_physdom->get_prob_lo();
  const Real time       = m_timestepper->get_time();
  const int num_species = m_plaskin->get_num_species();
  const int num_photons = m_plaskin->get_num_photons();

  RefCountedPtr<cdr_layout>& cdr = m_timestepper->get_cdr();
  RefCountedPtr<rte_layout>& rte = m_timestepper->get_rte();


  // This is all computes on volumetric centroids
  this->compute_cdr_densities(m_cdr_densities);
  this->compute_cdr_gradients(m_cdr_gradients);
  this->compute_E(m_E, m_grad_E);
  this->compute_rho(m_rho, m_grad_rho);
  this->compute_rte_densities(m_rte_densities);


  // Get maximum and minimum of everything
  Vector<Real> cdr_min(num_species), cdr_max(num_species);
  Vector<Real> rte_min(num_photons), rte_max(num_photons);
  Vector<Real> grad_cdr_max(num_species), grad_cdr_min(num_species);
  Real E_max, E_min;
  Real grad_E_max, grad_E_min;
  Real rho_max, rho_min;
  Real grad_rho_max, grad_rho_min;

  data_ops::get_max_min(cdr_max,           cdr_min,      m_cdr_densities);
  data_ops::get_max_min(rte_max,           rte_min,      m_rte_densities);
  for (int i = 0; i < m_cdr_gradients.size(); i++){
    data_ops::get_max_min_norm(grad_cdr_max[i], grad_cdr_min[i], m_cdr_gradients[i]);
  }
  data_ops::get_max_min_norm(E_max,        E_min,        m_E);
  data_ops::get_max_min_norm(grad_E_max,   grad_E_min,   m_grad_E);
  data_ops::get_max_min(rho_max,           rho_min,      m_rho, 0);
  data_ops::get_max_min_norm(grad_rho_max, grad_rho_min, m_grad_rho);

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box         = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      const EBCellFAB& E_fab    = (*m_E[lvl])[dit()];
      const EBCellFAB& gE_fab   = (*m_grad_E[lvl])[dit()];
      const EBCellFAB& rho_fab  = (*m_rho[lvl])[dit()];
      const EBCellFAB& grho_fab = (*m_grad_rho[lvl])[dit()];

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);


	// Electric field, rho, and their gradients
	const RealVect E        = RealVect(D_DECL(E_fab(vof, 0), E_fab(vof, 1), E_fab(vof, 2)));
	const RealVect grad_E   = RealVect(D_DECL(gE_fab(vof, 0), gE_fab(vof, 1), gE_fab(vof, 2)));
	const Real rho          = rho_fab(vof, 0);
	const RealVect grad_rho = RealVect(D_DECL(grho_fab(vof, 0), grho_fab(vof, 1), grho_fab(vof, 2)));

	Vector<Real>     cdr_densities;
	Vector<RealVect> cdr_gradients;
	Vector<Real>     rte_densities;
	cdr_densities.resize(m_plaskin->get_num_species());
	cdr_gradients.resize(m_plaskin->get_num_species());
	rte_densities.resize(m_plaskin->get_num_photons());

	for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
	  const int idx = solver_it.get_solver();
	  const EBCellFAB& density  = (*(m_cdr_densities[idx])[lvl])[dit()];
	  const EBCellFAB& gradient = (*(m_cdr_gradients[idx])[lvl])[dit()];
	  cdr_densities[idx] = density(vof, 0);
	  cdr_gradients[idx] = RealVect(D_DECL(gradient(vof, 0),
					       gradient(vof, 1),
					       gradient(vof, 2)));
	}

	for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
	  const int idx = solver_it.get_solver();
	  const EBCellFAB& density  = (*(m_rte_densities[idx])[lvl])[dit()];
	  rte_densities[idx] = density(vof, 0);
	}

	Vector<Real> tracers = this->tracer(pos,
					    time,
					    dx,
					    E,
					    E_min,
					    E_max,
					    grad_E,
					    grad_E_min,
					    grad_E_max,
					    rho,
					    rho_min,
					    rho_max,
					    grad_rho,
					    grad_rho_min,
					    grad_rho_max,
					    cdr_densities,
					    cdr_min,
					    cdr_max,
					    cdr_gradients,
					    grad_cdr_min,
					    grad_cdr_max,
					    rte_densities,
					    rte_min,
					    rte_max);
	
	for(int i = 0; i < m_num_tracers; i++){
	  (*m_tracer[i][lvl])[dit()](vof, 0) = tracers[i];
	}
      }
    }
  }


  for (int i = 0; i < m_num_tracers; i++){
    m_amr->average_down(m_tracer[i], m_phase);
    m_amr->interp_ghost(m_tracer[i], m_phase);
  }

  // Compute gradient of tracers
  for (int i = 0; i < m_num_tracers; i++){
    m_amr->compute_gradient(m_grad_tracer[i], m_tracer[i]);
    m_amr->average_down(m_grad_tracer[i], m_phase);
  }

  this->deallocate_storage(); // No reason to keep the extra storage lying around...
}

void full_tagger::compute_cdr_densities(Vector<EBAMRCellData>& a_cdr_densities){
  CH_TIME("full_tagger::compute_cdr_densities");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_cdr_densities" << endl;
  }

  RefCountedPtr<cdr_layout>& cdr = m_timestepper->get_cdr();

  for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();

    // Interpolate density to centroid
    data_ops::set_value(m_scratch, 0.0);
    data_ops::incr(m_scratch, solver->get_state(), 1.0);
    m_amr->interpolate_to_centroids(m_scratch, m_phase);

    // Copy to member data holder
    const int idx = solver_it.get_solver();
    data_ops::set_value(a_cdr_densities[idx], 0.0);
    data_ops::incr(a_cdr_densities[idx], m_scratch, 1.0);
  }
}

void full_tagger::compute_cdr_gradients(Vector<EBAMRCellData>& a_cdr_gradients){
  CH_TIME("full_tagger::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_cdr_gradients" << endl;
  }

  RefCountedPtr<cdr_layout>& cdr = m_timestepper->get_cdr();
    
  for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();

    const int idx = solver_it.get_solver();    
    m_amr->compute_gradient(a_cdr_gradients[idx], solver->get_state());
  }
}

void full_tagger::compute_E(EBAMRCellData& a_E, EBAMRCellData& a_grad_E){
  CH_TIME("full_tagger::compute_E");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_E" << endl;
  }

  m_timestepper->compute_E(a_E, m_phase);
  data_ops::vector_length(m_scratch, a_E);
  m_amr->compute_gradient(a_grad_E, m_scratch);

  m_amr->average_down(a_grad_E, m_phase);
  m_amr->interp_ghost(a_grad_E, m_phase);
  
  // Interpolate to centroids
  m_amr->interpolate_to_centroids(a_E,      m_phase);
  m_amr->interpolate_to_centroids(a_grad_E, m_phase);
}

void full_tagger::compute_rho(EBAMRCellData& a_rho, EBAMRCellData& a_grad_rho){
  CH_TIME("full_tagger::compute_rho");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rho" << endl;
  }

  // Compute cell-centered rho and its gradient
  m_timestepper->compute_rho(a_rho, m_phase);
  m_amr->compute_gradient(a_grad_rho, a_rho);

  // Transform to centroids
  m_amr->interpolate_to_centroids(a_rho,      m_phase);
  m_amr->interpolate_to_centroids(a_grad_rho, m_phase);
}

void full_tagger::compute_rte_densities(Vector<EBAMRCellData>& a_rte_densities){
  CH_TIME("full_tagger::compute_rte_densities");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rte_densities" << endl;
  }

  RefCountedPtr<rte_layout>& rte = m_timestepper->get_rte();

  for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();

    // Interpolate density to centroid
    data_ops::set_value(m_scratch, 0.0);
    data_ops::incr(m_scratch, solver->get_state(), 1.0);
    m_amr->interpolate_to_centroids(m_scratch, m_phase);

    // Copy to member data holder
    const int idx = solver_it.get_solver();
    data_ops::set_value(a_rte_densities[idx], 0.0);
    data_ops::incr(a_rte_densities[idx], m_scratch, 1.0);
  }
}
