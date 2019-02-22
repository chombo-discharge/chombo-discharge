/*!
  @file   sisdc_tagger.cpp
  @brief  Implementation of sisdc_tagger.H
  @author Robert Marskar
  @date   Oct. 2018
*/

#include "sisdc_tagger.H"
#include "data_ops.H"
#include "sisdc.H"

#include <ParmParse.H>
  
sisdc_tagger::sisdc_tagger(){
  m_num_tracers = 1;
  m_name        = "sisdc_tagger";

  m_electron_idx = 0;
  m_thresh       = 0.1;

  ParmParse pp("sisdc_tagger");
  pp.query("electron_index",  m_electron_idx);
  pp.query("err_thresh",      m_thresh);

  this->set_phase(phase::gas);

  std::cout << m_thresh << std::endl;
}

sisdc_tagger::~sisdc_tagger(){
  
}

void sisdc_tagger::compute_tracers(){
  CH_TIME("sisdc_tagger::compute_tracers");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_tracers" << endl;
  }

  const int comp = 0;
  const int max_amr_depth = m_amr->get_max_amr_depth();
  const int finest_level  = m_amr->get_finest_level();

  // Get electron error
  sisdc* stepper = (sisdc*) (&(*m_timestepper));
  EBAMRCellData& ne_err  = *(stepper->get_cdr_errors()[m_electron_idx]);

  // Get maximum and minimum ne and Se, and the electric field
  Real err_max,  err_min;
  data_ops::get_max_min(err_max,  err_min,  ne_err,   comp);

  err_max = Max(Abs(err_max), Abs(err_min));

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box         = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      const EBCellFAB& ne_err_fab = (*ne_err[lvl])[dit()];
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa = ebisbox.volFrac(vof);
	const Real err_vof  = ne_err_fab(vof, comp);
	(*m_tracer[0][lvl])[dit()](vof, 0) = (err_max > 0.0) ? Abs(kappa*err_vof/err_max) : 0.0;
      }
    }
  }

  // Compute gradient of the tracer
  for (int i = 0; i < m_num_tracers; i++){
    m_amr->average_down(m_tracer[i], m_phase);
    m_amr->interp_ghost(m_tracer[i], m_phase);
  }
}


bool sisdc_tagger::coarsen_cell(const RealVect&         a_pos,
				   const Real&             a_time,
				   const Real&             a_dx,
				   const int&              a_lvl,
				   const Vector<Real>&     a_tracer,
				   const Vector<RealVect>& a_grad_tracer){

  return true;
}

bool sisdc_tagger::refine_cell(const RealVect&         a_pos,
				  const Real&             a_time,
				  const Real&             a_dx,
				  const int&              a_lvl,
				  const Vector<Real>&     a_tracer,
				  const Vector<RealVect>& a_grad_tracer){
  const bool refine = (a_tracer[0] > m_thresh) ? true : false;
  
  return refine;
}
