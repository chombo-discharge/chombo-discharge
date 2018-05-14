/*!
  @file   ultralw_tagger.cpp
  @brief  Implementation of ultralw_tagger.H
  @author Robert Marskar
  @date   May. 2018
*/

#include "ultralw_tagger.H"
#include "data_ops.H"

#include <EBArith.H>

ultralw_tagger::ultralw_tagger(const int a_num_tracers){
  CH_TIME("ultralw_tagger::ultralw_tagger");
  this->set_verbosity(-1);
  if(m_verbosity > 5){
    pout() << "ultralw_tagger::ultralw_tagger" << endl;
  }
  
  m_num_tracers = a_num_tracers;
  m_name        = "cell_tagger";


  this->set_phase(phase::gas);
}

ultralw_tagger::~ultralw_tagger(){

}

void ultralw_tagger::allocate_storage(){
  CH_TIME("ultralw_tagger::allocate_storage");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_storage" << endl;
  }

  m_amr->allocate(m_scratch,  m_phase, 1);
  m_amr->allocate(m_E,        m_phase, SpaceDim);
  m_amr->allocate(m_grad_E,   m_phase, SpaceDim);
}

void ultralw_tagger::deallocate_storage(){
  CH_TIME("ultralw_tagger::deallocate_storage");
  if(m_verbosity > 5){
    pout() << m_name + "::deallocate_storage" << endl;
  }

  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_E);
  m_amr->deallocate(m_grad_E);
}

void ultralw_tagger::compute_E(EBAMRCellData& a_E, EBAMRCellData& a_grad_E){
  CH_TIME("ultralw_tagger::compute_E");
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

void ultralw_tagger::compute_tracers(){
  CH_TIME("ultralw_tagger::compute_tracers");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_tracers" << endl;
  }

  this->allocate_storage();
  
  const RealVect origin = m_physdom->get_prob_lo();
  const Real time       = m_timestepper->get_time();

  // Compute electric field on volumetric centroids
  this->compute_E(m_E, m_grad_E);

  // Get maximum and minimum of everything
  Real E_max, E_min;
  Real grad_E_max, grad_E_min;

  data_ops::get_max_min_norm(E_max,        E_min,        m_E);
  data_ops::get_max_min_norm(grad_E_max,   grad_E_min,   m_grad_E);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
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

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);

	// Electric field and grad(|E|)
	const RealVect E        = RealVect(D_DECL(E_fab(vof, 0), E_fab(vof, 1), E_fab(vof, 2)));
	const RealVect grad_E   = RealVect(D_DECL(gE_fab(vof, 0), gE_fab(vof, 1), gE_fab(vof, 2)));

	Vector<Real> tracers = this->tracer(pos,
					    time,
					    dx,
					    E,
					    E_min,
					    E_max,
					    grad_E,
					    grad_E_min,
					    grad_E_max);
	
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
