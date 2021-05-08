/*!
  @file   cdr_plasma_field_tagger.cpp
  @brief  Implementation of cdr_plasma_field_tagger.H
  @author Robert Marskar
  @date   May. 2018
*/

#include "cdr_plasma_field_tagger.H"
#include "data_ops.H"

#include <EBArith.H>

#include "CD_NamespaceHeader.H"
using namespace physics::cdr_plasma;

cdr_plasma_field_tagger::cdr_plasma_field_tagger(){
  CH_TIME("cdr_plasma_field_tagger::cdr_plasma_field_tagger");
  if(m_verbosity > 5){
    pout() << "cdr_plasma_field_tagger::cdr_plasma_field_tagger" << endl;
  }

  m_name = "cdr_plasma_field_tagger";
}

cdr_plasma_field_tagger::~cdr_plasma_field_tagger(){

}

void cdr_plasma_field_tagger::allocate_storage(){
  CH_TIME("cdr_plasma_field_tagger::allocate_storage");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_storage" << endl;
  }

  m_amr->allocate(m_scratch,  m_realm, m_phase, 1);
  m_amr->allocate(m_E,        m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_grad_E,   m_realm, m_phase, SpaceDim);
}

void cdr_plasma_field_tagger::deallocate_storage(){
  CH_TIME("cdr_plasma_field_tagger::deallocate_storage");
  if(m_verbosity > 5){
    pout() << m_name + "::deallocate_storage" << endl;
  }

  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_E);
  m_amr->deallocate(m_grad_E);
}

void cdr_plasma_field_tagger::compute_E(EBAMRCellData& a_E, EBAMRCellData& a_grad_E){
  CH_TIME("cdr_plasma_field_tagger::compute_E");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_E" << endl;
  }

  m_timestepper->compute_E(a_E, m_phase);
  data_ops::vector_length(m_scratch, a_E);
  m_amr->compute_gradient(a_grad_E, m_scratch, m_realm, phase::gas);

  m_amr->average_down(a_grad_E, m_realm, m_phase);
  m_amr->interp_ghost(a_grad_E, m_realm, m_phase);
  
  // Interpolate to centroids
  m_amr->interpolate_to_centroids(a_E,      m_realm, m_phase);
  m_amr->interpolate_to_centroids(a_grad_E, m_realm, m_phase);
}

void cdr_plasma_field_tagger::compute_tracers(){
  CH_TIME("cdr_plasma_field_tagger::compute_tracers");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_tracers" << endl;
  }

  this->allocate_storage();
  
  const RealVect origin = m_amr->get_prob_lo();
  const Real time       = m_timestepper->get_time();

  // Compute electric field on volumetric centroids
  this->compute_E(m_E, m_grad_E);

  // Get maximum and minimum of everything
  Real E_max, E_min;
  Real grad_E_max, grad_E_min;

  data_ops::get_max_min_norm(E_max,        E_min,        m_E);
  data_ops::get_max_min_norm(grad_E_max,   grad_E_min,   m_grad_E);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box         = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      const EBCellFAB& E_fab    = (*m_E[lvl])[dit()];
      const EBCellFAB& gE_fab   = (*m_grad_E[lvl])[dit()];

      const BaseFab<Real>& E_reg  = E_fab.getSingleValuedFAB();
      const BaseFab<Real>& gE_reg = gE_fab.getSingleValuedFAB();

      // Avoid the extra point lookups by getting these before the point loops
      Vector<EBCellFAB*> tr;
      Vector<BaseFab<Real>* > tr_fab;
      for (int i = 0; i < m_num_tracers; i++){
	tr.push_back(&((*m_tracer[i][lvl])[dit()]));
	tr_fab.push_back(&(tr[i]->getSingleValuedFAB()));
      }

      // Regular box loop
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv   = bit();
	const RealVect pos = origin + RealVect(iv)*dx;

	const RealVect E        = RealVect(D_DECL(E_reg(iv, 0),  E_reg(iv, 1),  E_reg(iv, 2)));
	const RealVect grad_E   = RealVect(D_DECL(gE_reg(iv, 0), gE_reg(iv, 1), gE_reg(iv, 2)));

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
	  (*tr_fab[i])(iv, 0) = tracers[i];
	}
      }

      // Irregular box loop
      const IntVectSet& irreg = ebisbox.getIrregIVS(box);
      for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);

	// Electric field and grad(|E|)
	const RealVect E        = RealVect(D_DECL(E_fab(vof, 0),  E_fab(vof, 1), E_fab(vof, 2)));
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
	  (*tr[i])(vof, 0);
	}
      }
    }
  }


  for (int i = 0; i < m_num_tracers; i++){
    m_amr->average_down(m_tracer[i], m_realm, m_phase);
    m_amr->interp_ghost(m_tracer[i], m_realm, m_phase);
  }

  // Compute gradient of tracers
  for (int i = 0; i < m_num_tracers; i++){
    m_amr->compute_gradient(m_grad_tracer[i], m_tracer[i], m_realm, phase::gas);
    m_amr->average_down(m_grad_tracer[i], m_realm, m_phase);
  }

  this->deallocate_storage(); // No reason to keep the extra storage lying around...
}
#include "CD_NamespaceFooter.H"
