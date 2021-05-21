/*!
  @file   CdrPlasmaStreamerTagger.cpp
  @brief  Implementation of CdrPlasmaStreamerTagger.H
  @author Robert Marskar
  @date   Oct. 2018
*/

#include <CD_CdrPlasmaStreamerTagger.H>
#include <CD_DataOps.H>

#include <ParmParse.H>
  
CdrPlasmaStreamerTagger::CdrPlasmaStreamerTagger(){
  m_num_tracers = 2;
  m_name        = "CdrPlasmaStreamerTagger";

  m_electron_idx   = 0;
  m_fudge          = 1.E-4;
  m_thresh1.resize(1, 0.5);
  m_thresh2.resize(1, 0.1);




  {
    ParmParse pp("CdrPlasmaStreamerTagger");
    
    const int num_src = pp.countval("refine_source");
    const int num_crv = pp.countval("refine_E_curv");

    pp.query("relative_factor", m_fudge);
    pp.query("electron_index",  m_electron_idx);
    if(num_src > 0){
      m_thresh1.resize(num_src);
      pp.getarr("refine_source", m_thresh1, 0, num_src);
    }
    if(num_crv > 0){
      m_thresh2.resize(num_crv);
      pp.getarr("refine_E_curv", m_thresh2, 0, num_crv);
    }
  }

  this->setPhase(phase::gas);
}

CdrPlasmaStreamerTagger::~CdrPlasmaStreamerTagger(){
  
}

void CdrPlasmaStreamerTagger::computeTracers(){
  CH_TIME("CdrPlasmaStreamerTagger::computeTracers");
  if(m_verbosity > 5){
    pout() << m_name + "::computeTracers" << endl;
  }

  const int comp = 0;
  const int max_amr_depth = m_amr->getMaxAmrDepth();
  const int finest_level  = m_amr->getFinestLevel();
  
  while (m_thresh1.size() <= max_amr_depth){
    m_thresh1.push_back(m_thresh1.back());
  }
  while (m_thresh2.size() < max_amr_depth){
    m_thresh2.push_back(m_thresh2.back());
  }



  // Get electron stuff
  EBAMRCellData& ne  = (m_timeStepper->get_cdr())->getSolvers()[m_electron_idx]->getPhi();
  EBAMRCellData& Se  = (m_timeStepper->get_cdr())->getSolvers()[m_electron_idx]->getSource();
  EBAMRCellData& ve  = (m_timeStepper->get_cdr())->getSolvers()[m_electron_idx]->getCellCenteredVelocity();

  // Compute the electric field
  EBAMRCellData rho;
  EBAMRCellData Efield;
  m_amr->allocate(Efield,  m_phase, SpaceDim);
  m_timeStepper->compute_E(Efield, m_phase);
  m_amr->allocatePointer(rho);
  m_amr->alias(rho, phase::gas, m_timeStepper->get_poisson()->getSource());

  // Compute the electric field magnitude
  EBAMRCellData Emag;
  m_amr->allocate(Emag, m_phase, 1);
  DataOps::vectorLength(Emag, Efield);
  m_amr->averageDown(Emag, m_phase);
  m_amr->interpGhost(Emag, m_phase);



  // Get maximum and minimum ne and Se, and the electric field
  Real ne_max,  ne_min;
  Real Se_max,  Se_min;
  Real E_max,   E_min;
  Real rho_max, rho_min;
  DataOps::getMaxMin(ne_max,  ne_min,  ne,   comp);
  DataOps::getMaxMin(Se_max,  Se_min,  Se,   comp);
  DataOps::getMaxMin(E_max,   E_min,   Emag, comp);
  DataOps::getMaxMin(rho_max, rho_min, rho,  comp);

  // Compute the FLASH code error
  const Real FLASH_eps = 1.E-2;
  if(ne_max > 1.E-2*ne_min && Abs(ne_max) > 0.0){
    DataOps::flashError(m_tracer[0], ne, FLASH_eps);
  }
  else{
    DataOps::setValue(m_tracer[0], 0.0);
  }

  //  if(Abs(rho_max) > 1.E-2*Abs(rho_min)){
  //    DataOps::flashError(m_tracer[1], Emag, FLASH_eps);
  // }
  // else{
  //   DataOps::setValue(m_tracer[1], 0.0);
  // }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box         = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      const EBCellFAB& ne_fab  = (*ne[lvl])[dit()];
      const EBCellFAB& Se_fab  = (*Se[lvl])[dit()];
      const EBCellFAB& rho_fab = (*rho[lvl])[dit()];
      const EBCellFAB& E_fab   = (*Emag[lvl])[dit()];

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real S        = Se_fab(vof, comp);
	const Real n        = ne_fab(vof, comp);
	const Real ro       = rho_fab(vof, comp);
	const Real e        = E_fab(vof,comp);

	Real& tracer0 = (*m_tracer[0][lvl])[dit()](vof, 0);
	Real& tracer1 = (*m_tracer[1][lvl])[dit()](vof, 0);
	
	if(S > 1.E-3*Se_max && n > 1.E-3*ne_max){
	  tracer0 = Abs(tracer0);
#if 1
	  tracer0 = S/Se_max;
#endif
	}
	else{
	  tracer0 = 0.0;
	}

#if 1
	if(e > 1.E-1*E_max && e > 0.0){
	  tracer1 = Abs(e);
	}
	else{
	  tracer1 = 0.0;
	}
#endif
      }
    }
  }

  // Compute gradient of the tracer
  for (int i = 0; i < m_num_tracers; i++){
    m_amr->averageDown(m_tracer[i], m_phase);
    m_amr->interpGhost(m_tracer[i], m_phase);
    m_amr->computeGradient(m_grad_tracer[i], m_tracer[i]);
    m_amr->averageDown(m_grad_tracer[i], m_phase);
  }
}


bool CdrPlasmaStreamerTagger::coarsenCell(const RealVect&         a_pos,
				   const Real&             a_time,
				   const Real&             a_dx,
				   const int&              a_lvl,
				   const Vector<Real>&     a_tracer,
				   const Vector<RealVect>& a_grad_tracer){

  const bool coarsen1 = a_tracer[0] < 0.125*m_thresh1[a_lvl] ? true : false;
  const bool coarsen2 = a_tracer[1] < 0.25*m_thresh2[a_lvl] ? true : false;
  return true;
  return coarsen1;// && coarsen2; 
}

bool CdrPlasmaStreamerTagger::refineCell(const RealVect&         a_pos,
				  const Real&             a_time,
				  const Real&             a_dx,
				  const int&              a_lvl,
				  const Vector<Real>&     a_tracer,
				  const Vector<RealVect>& a_grad_tracer){
  const bool refine1 = (a_tracer[0] > m_thresh1[a_lvl]) ? true : false;
  //  const bool refine2 = (a_tracer[1] > m_thresh2[a_lvl]) ? true : false;
  const bool refine2 = a_grad_tracer[1].vectorLength()*a_dx/a_tracer[1] > m_thresh2[a_lvl] ? true : false;

  return refine1 || refine2;
#include "CD_NamespaceFooter.H"
