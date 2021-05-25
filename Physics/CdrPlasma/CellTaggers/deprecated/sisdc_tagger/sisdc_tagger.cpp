/*!
  @file   sisdc_tagger.cpp
  @brief  Implementation of sisdc_tagger.H
  @author Robert Marskar
  @date   Oct. 2018
*/

#include "sisdc_tagger.H"
#include <CD_DataOps.H>
#include "sisdc.H"

#include <ParmParse.H>
  
sisdc_tagger::sisdc_tagger(){
  m_num_tracers = 2;
  m_name        = "sisdc_tagger";

  m_cdr_idx     = 0;
  m_err_thresh  = 0.01;
  m_curv_thresh = 0.1;
  m_mag_thresh  = 0.9;

  ParmParse pp("sisdc_tagger");
  pp.query("cdr_index",        m_cdr_idx);
  pp.query("err_thresh",       m_err_thresh);
  pp.query("curv_thresh",      m_curv_thresh);
  pp.query("magnitude_thresh", m_mag_thresh);

  this->setPhase(phase::gas);
}

sisdc_tagger::~sisdc_tagger(){
  
}

void sisdc_tagger::computeTracers(){
  CH_TIME("sisdc_tagger::computeTracers");
  if(m_verbosity > 5){
    pout() << m_name + "::computeTracers" << endl;
  }

  const int comp = 0;
  const int max_amr_depth = m_amr->getMaxAmrDepth();
  const int finest_level  = m_amr->getFinestLevel();

  // Get electron density and error
  sisdc* stepper = (sisdc*) (&(*m_timeStepper));
  EBAMRCellData& ne_err  = *(stepper->getCdrErrors()[m_cdr_idx]);
  EBAMRCellData& ne      = *stepper->getCdrSolvers()->getPhis()[m_cdr_idx];

  // Get maximum and minimum stuff
  Real err_max,  err_min, Emax, Emin;
  Real ne_min, ne_max;
  DataOps::getMaxMin(err_max,  err_min,  ne_err,   comp);
  DataOps::getMaxMin(ne_max,  ne_min,  ne,   comp);
  err_max = Max(Abs(err_max), Abs(err_min));

  // Compute the electric field and take the magnitude onto tracer1. Also scale it. 
  EBAMRCellData E;
  m_amr->allocate(E, phase::gas, SpaceDim);
  m_timeStepper->compute_E(E, phase::gas);

  DataOps::vectorLength(m_tracer[1], E);
  m_amr->interpToCentroids(m_tracer[1], m_phase);
  DataOps::getMaxMin(Emax, Emin, m_tracer[1], 0);
  DataOps::scale(m_tracer[1], 1./Emax);


  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box         = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      const EBCellFAB& ne_fab     = (*ne[lvl])[dit()];
      const EBCellFAB& ne_err_fab = (*ne_err[lvl])[dit()];
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa = ebisbox.volFrac(vof);
	const Real ne_vof   = ne_fab(vof, comp);
	const Real err_vof  = ne_err_fab(vof, comp);
	(*m_tracer[0][lvl])[dit()](vof, 0) = (err_max > 0.0) ? Abs(kappa*err_vof/err_max) : 0.0;
      }
    }
  }

  DataOps::getMaxMin(err_max, err_min, m_tracer[0], 0);
  err_max = Max(Abs(err_max), Abs(err_min));
  DataOps::scale(m_tracer[0], 1./err_max);
}

bool sisdc_tagger::coarsenCell(const RealVect&         a_pos,
				const Real&             a_time,
				const Real&             a_dx,
				const int&              a_lvl,
				const Vector<Real>&     a_tracer,
				const Vector<RealVect>& a_grad_tracer){
  return true;
}

bool sisdc_tagger::refineCell(const RealVect&         a_pos,
			       const Real&             a_time,
			       const Real&             a_dx,
			       const int&              a_lvl,
			       const Vector<Real>&     a_tracer,
			       const Vector<RealVect>& a_grad_tracer){
  const bool refine_err  = (Abs(a_tracer[0]) > m_err_thresh) ? true : false;
  const bool refine_curv = (a_grad_tracer[1].vectorLength()*a_dx)/a_tracer[1] > m_curv_thresh ? true : false;
  const bool refine_magn = a_tracer[1] > m_mag_thresh;
  
  return refine_err || refine_magn || refine_curv;
#include "CD_NamespaceFooter.H"
