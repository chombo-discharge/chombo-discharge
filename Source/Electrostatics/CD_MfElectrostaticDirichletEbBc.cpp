/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MfElectrostaticDirichletEbBc.cpp
  @brief  Implementation of CD_MfElectrostaticDirichletEbBc.H
  @author Robert Marskar
*/

#include <CD_LeastSquares.H>
#include <CD_MfElectrostaticDirichletEbBc.H>
#include <ParmParse.H>

#include <CD_NamespaceHeader.H>

bool MfElectrostaticDirichletEbBc::s_areaFracWeighted = false;
bool MfElectrostaticDirichletEbBc::s_quadrant_based   = true;
int  MfElectrostaticDirichletEbBc::s_lsq_radius       = 1;

#define DEBUG 0

MfElectrostaticDirichletEbBc::MfElectrostaticDirichletEbBc(const ProblemDomain& a_domain,
							   const EBISLayout&    a_ebisl,
							   const RealVect&      a_dx,
							   const IntVect*       a_phig,
							   const IntVect*       a_rhsg,
							   const int            a_phase) {
  m_domain = a_domain;
  m_ebisl  = a_ebisl;
  m_dx     = a_dx;
  m_phase  = a_phase;
  m_order  = 1;
}

MfElectrostaticDirichletEbBc::~MfElectrostaticDirichletEbBc(){

}

LayoutData<BaseIVFAB<VoFStencil> >* MfElectrostaticDirichletEbBc::getFluxStencil(int ivar){
  return &m_IrregStencils;
}

void MfElectrostaticDirichletEbBc::setJumpObject(const RefCountedPtr<JumpBc> a_jumpbc){
  m_jumpbc = a_jumpbc;
}

void MfElectrostaticDirichletEbBc::setOrder(int a_order){
  m_order = a_order;
}

void MfElectrostaticDirichletEbBc::defineIVS(const MFLevelGrid& a_mflg){
    
  m_ivs.define(a_mflg.getGrids());

  const DisjointBoxLayout& dbl = a_mflg.getGrids();

  int num = 0;
  m_ivs.define(dbl);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    m_ivs[dit()] = a_mflg.interfaceRegion(dbl.get(dit()), dit());
    num += m_ivs[dit()].numPts();
  }
 
  m_definedmf = true;
}

void MfElectrostaticDirichletEbBc::define(const LayoutData<IntVectSet>& a_cfivs, const Real& a_factor){

  if(!m_definedmf){
    MayDay::Error("MfElectrostaticDirichletEbBc::define - must set multifluid cells first!");
  }
  else if(!m_coefSet){
    MayDay::Error("MfElectrostaticDirichletEbBc::define - must call setCoef BEFORE calling define");
  }
  
  const int comp      = 0;
  const int num_comps = 1;

  const int otherphase = m_phase == 0 ? 1 : 0;

  const DisjointBoxLayout& dbl = m_ebisl.getDisjointLayout();

  m_irreg_weights.define(dbl);
  m_matching_weights.define(dbl);
  m_IrregStencils.define(dbl);
  m_matching_stencils.define(dbl);
  m_vofit_irreg.define(dbl);
  m_vofit_diri.define(dbl);
  m_vofit_matching.define(dbl);

  const LayoutData<MFInterfaceFAB<VoFStencil> >& jumpstens = m_jumpbc->getStencils();
  const LayoutData<MFInterfaceFAB<Real> >& jumpweights     = m_jumpbc->getWeights();
  const LayoutData<MFInterfaceFAB<Real> >& jump_coef       = m_jumpbc->getAvgBco();

  // Create stencils on irregular cells that are not part of m_ivs
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box          = dbl[dit()];
    const EBISBox& ebisbox  = m_ebisl[dit()];
    const EBGraph& ebgraph  = ebisbox.getEBGraph();
    const IntVectSet& cfivs = a_cfivs[dit()];

    VoFIterator& vofit_irreg    = m_vofit_irreg[dit()];
    VoFIterator& vofit_diri     = m_vofit_diri[dit()];
    VoFIterator& vofit_matching = m_vofit_matching[dit()];
    
    IntVectSet irreg_ivs;  // All irregular cells
    IntVectSet diri_ivs;   // Pure dirichlet cells (i.e. non-matched irregular cells)
    IntVectSet match_ivs;  // Quasi-dirichlet cells (i.e. matched irregular cells)
    
    irreg_ivs  = ebisbox.getIrregIVS(box);
    irreg_ivs |= ebisbox.getMultiCells(box);

    diri_ivs = irreg_ivs - m_ivs[dit()];
    match_ivs = m_ivs[dit()];

    m_irreg_weights[dit()].define(irreg_ivs,  ebgraph, num_comps);
    m_IrregStencils[dit()].define(irreg_ivs, ebgraph, num_comps);

    vofit_irreg.define(irreg_ivs, ebgraph);
    vofit_diri.define(diri_ivs, ebgraph);
    vofit_matching.define(match_ivs, ebgraph);

    // Build stencils for pure Dirichlet type irreg cells
    for (vofit_irreg.reset(); vofit_irreg.ok(); ++vofit_irreg){
      const VolIndex& vof   = vofit_irreg();

      Real& cur_weight        = m_irreg_weights[dit()](vof, comp);
      VoFStencil& cur_stencil = m_IrregStencils[dit()](vof, comp);

      bool drop_order = false;

      if(m_order == 2){
	drop_order = this->getSecondOrderStencil(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	if(drop_order){
#if DEBUG
	  pout() << "dropping order on domain = " << m_domain << " and cell = " << vof.gridIndex() << endl;
#endif
	  this->getFirstOrderStencil(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}
      }
      else if(m_order == 1){
	this->getFirstOrderStencil(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }
    }

    // Adjust stencils for matching cells
    const MFInterfaceFAB<VoFStencil>& avgStens = m_jumpbc->getAvgStencils()[dit()];
    const MFInterfaceFAB<Real>& avgWeights     = m_jumpbc->getAvgWeights()[dit()];
    const MFInterfaceFAB<Real>& avgBco         = m_jumpbc->get_avgBco()[dit()];

    for (vofit_matching.reset(); vofit_matching.ok(); ++vofit_matching){
      const VolIndex& vof = vofit_matching();
      const VolIndex vof0(vof.gridIndex(), 0); // Multi-cell averages stored on first component

      // Everything that comes from JumpBc is from vof0, as it should!!!
      VoFStencil& cur_stencil     = m_IrregStencils[dit()](vof, comp);
      const VoFStencil& jump_sten = avgStens.getIVFAB(m_phase)(vof0, comp);
      const Real cur_weight       = m_irreg_weights[dit()](vof, comp); 
      const Real wp               = avgWeights.getIVFAB(m_phase)(vof0, comp);
      const Real wq               = avgWeights.getIVFAB(otherphase)(vof0, comp);
      const Real bp               = avgBco.getIVFAB(m_phase)(vof0, comp);
      const Real bq               = avgBco.getIVFAB(otherphase)(vof0, comp);
      const Real factor           = cur_weight*bp/(bp*wp + bq*wq);

      
#if DEBUG // Debug. Check if anything is NaN
      for (int i = 0; i < cur_stencil.size(); i++){
	const Real weight = cur_stencil.weight(i);
	if(weight != weight) MayDay::Abort("MfElectrostaticDirichletEbBc::define - got NaN weight in BC stencil");
      }
      for (int i = 0; i < jump_sten.size(); i++){
	const Real weight = jump_sten.weight(i);
	if(weight != weight) MayDay::Abort("MfElectrostaticDirichletEbBc::define - got NaN weight in matching stencil");
      }
      if(cur_weight != cur_weight) MayDay::Abort("MfElectrostaticDirichletEbBc::define - cur_weight is NaN");
      if(bp         != bp)         MayDay::Abort("MfElectrostaticDirichletEbBc::define - bp is NaN");
      if(bq         != bq)         MayDay::Abort("MfElectrostaticDirichletEbBc::define - bq is NaN");
      if(wp         != wp)         MayDay::Abort("MfElectrostaticDirichletEbBc::define - wp is NaN");
      if(wq         != wq)         MayDay::Abort("MfElectrostaticDirichletEbBc::define - wq is NaN");
      if(factor     != factor)     MayDay::Abort("MfElectrostaticDirichletEbBc::define - factor is NaN");
#endif

      // Adjust stencil for BC jump
      VoFStencil addsten(jump_sten);
      addsten     *= -1.0*factor;
      cur_stencil += addsten;


      // Debug after adjust. 
#if DEBUG
      for (int i = 0; i < cur_stencil.size(); i++){
	const Real weight = cur_stencil.weight(i);
	if(weight != weight) MayDay::Abort("MfElectrostaticDirichletEbBc::define - got NaN weight");
      }
#endif
    }

    // Scale stencils appropriately. They should be scaled by beta*bco*area_frac*a_factor
    for (vofit_irreg.reset(); vofit_irreg.ok(); ++vofit_irreg){
      const VolIndex& vof   = vofit_irreg();
      const Real& area_frac = ebisbox.bndryArea(vof);
      const Real factor     = m_beta*(*m_bcoe)[dit()](vof, comp)*area_frac*a_factor;
      
      VoFStencil& cur_stencil = m_IrregStencils[dit()](vof, comp);
      cur_stencil *= factor;

      if(MfElectrostaticDirichletEbBc::s_areaFracWeighted){
	cur_stencil *= ebisbox.areaFracScaling(vof);
      }
    }
  }
}

bool MfElectrostaticDirichletEbBc::getSecondOrderStencil(Real&             a_weight,
							 VoFStencil&       a_stencil,
							 const VolIndex&   a_vof,
							 const EBISBox&    a_ebisbox,
							 const IntVectSet& a_cfivs){
  CH_TIME("MfElectrostaticDirichletEbBc::getSecondOrderStencil");

  a_stencil.clear();
  bool drop_order = false;

  Vector<VoFStencil> point_stencils;
  Vector<Real> distance_along_lines;
  
  EBArith::johanStencil(drop_order, point_stencils, distance_along_lines, a_vof, a_ebisbox, m_dx, a_cfivs);
  if(drop_order){
    return true;
  }

  // If we got this far we have a stencil
  CH_assert(distance_along_lines.size() >= 2);
  CH_assert(point_stencils.size() >= 2);

  const Real& x1   = distance_along_lines[0];
  const Real& x2   = distance_along_lines[1];
  const Real denom = x2*x2*x1 - x1*x1*x2;

  VoFStencil& phi1Sten = point_stencils[0];
  VoFStencil& phi2Sten = point_stencils[1];
  
  phi1Sten *= -x2*x2/denom;
  phi2Sten *=  x1*x1/denom;

  a_weight   = -x1*x1/denom + x2*x2/denom;
  a_stencil +=  phi1Sten;
  a_stencil +=  phi2Sten;

  return false;
}

void MfElectrostaticDirichletEbBc::getFirstOrderStencil(Real&             a_weight,
							VoFStencil&       a_stencil,
							const VolIndex&   a_vof,
							const EBISBox&    a_ebisbox,
							const IntVectSet& a_cfivs){
  CH_TIME("MfElectrostaticDirichletEbBc::getFirstOrderStencil");

  const RealVect normal   = a_ebisbox.normal(a_vof);
  const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);

  if(s_quadrant_based){
    EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisbox, m_dx, m_domain, 0);
  }
  else{
    EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil,
					       a_weight,
					       normal,
					       centroid,
					       a_vof,
					       a_ebisbox,
					       m_dx*RealVect::Unit,
					       m_domain,
					       0,
					       s_lsq_radius);
  }

  // Damn, couldn't find a stencil. Try including more cells
  if(s_quadrant_based && a_stencil.size() == 0){
    EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil,
					       a_weight,
					       normal,
					       centroid,
					       a_vof,
					       a_ebisbox,
					       m_dx*RealVect::Unit,
					       m_domain,
					       0,
					       s_lsq_radius);
  }

  // Still haven't found a stencil. Approximate using one side only. If both sides do this then we can't solve. 
  if(a_stencil.size() == 0){
#if DEBUG
    pout() << "MfElectrostaticDirichletEbBc::get_first_order sten - no sten on domain = "
	   << m_domain << "\t vof = "
	   << a_vof.gridIndex() << endl;
    MayDay::Warning("MfElectrostaticDirichletEbBc::getFirstOrderStencil - could not find a stencil. ");
#endif

    // Do an approximation for the gradient at the boundary by taking the gradient at the cell center and dotting it with the normal vector
    for (int dir = 0; dir < SpaceDim; dir++){
      VoFStencil sten; 
      EBArith::getFirstDerivStencilWidthOne(sten, a_vof, a_ebisbox, dir, m_dx[dir], (IntVectSet*) &a_cfivs, 0);

      if(sten.size() > 0){
	sten      *= normal[dir];
	a_stencil += sten;
      }
    }
    a_weight = 0.0;
  }

#if 0 // Test LeastSquares
  ParmParse pp("lsq");

  int p = 0;
  int order = 1;
  bool useLSQ = false;
  pp.query("use", useLSQ);
  pp.query("p", p);
  pp.query("order", order);
  

  if(useLSQ){
    VoFStencil mySten = LeastSquares::getBndryGradSten(a_vof, a_ebisbox, m_dx[0], p, order);

    if(mySten.size() == 0){
      mySten = LeastSquares::getBndryGradSten(a_vof, a_ebisbox, m_dx[0], p, 1);

      if(mySten.size() == 0){
    	std::cout << "did not find lsq gradient stencil" << std::endl;
      }
    }

    // This is how you do the weights for the boundary point
    if(mySten.size() > 0){
      a_stencil =  LeastSquares::projectGradSten(mySten, -normal);
      a_weight  = -LeastSquares::sumAllWeights(a_stencil);
    }
  }
#endif

}

void MfElectrostaticDirichletEbBc::applyEBFlux(EBCellFAB&                    a_lphi,
					       const EBCellFAB&              a_phi,
					       VoFIterator&                  a_vofit,
					       const LayoutData<IntVectSet>& a_cfivs,
					       const DataIndex&              a_dit,
					       const RealVect&               a_probLo,
					       const RealVect&               a_dx,
					       const Real&                   a_factor,
					       const bool&                   a_useHomogeneous,
					       const Real&                   a_time){

  const int comp         = 0;
  const int otherphase   = m_phase == 0 ? 1 : 0;
  const EBISBox& ebisbox = a_phi.getEBISBox();

  const MFInterfaceFAB<Real>& inhomo = m_jumpbc->getInhomogeneous()[a_dit];
  const MFInterfaceFAB<Real>& homog  = m_jumpbc->getHomogeneous()[a_dit];

  // Multi-fluid cells
  VoFIterator& vofit_matching = m_vofit_matching[a_dit];
  for (vofit_matching.reset(); vofit_matching.ok(); ++vofit_matching){
    const VolIndex& vof   = vofit_matching();
    const Real weight     = m_irreg_weights[a_dit](vof, comp);
    const Real beta       = m_beta;
    const Real bco        = (*m_bcoe)[a_dit](vof, comp);
    const Real area_frac  = ebisbox.bndryArea(vof);


    Real value = inhomo.getIVFAB(m_phase)(vof, comp); // This is the flux from the other side. Always use it. 
    if(!a_useHomogeneous){
      value += homog.getIVFAB(m_phase)(vof,comp);
    }

    // Flux - this contains the contribution from the surface charge (inhomogeneous bc) and the
    // other side of the interface (quasi-homogeneous bc). 
    Real flux = beta*weight*value*bco*area_frac*a_factor;

    if(MfElectrostaticDirichletEbBc::s_areaFracWeighted){
      flux *= ebisbox.areaFracScaling(vof);
    }

    // Increment
    a_lphi(vof, comp) += flux;
  }

  // Pure Dirichlet irregular cells
  if(!a_useHomogeneous){
    const EBISBox& ebisbox = a_lphi.getEBISBox();

    VoFIterator& vofit_diri = m_vofit_diri[a_dit];
    for (vofit_diri.reset(); vofit_diri.ok(); ++vofit_diri){
      const VolIndex& vof = vofit_diri();

      const Real& value      = (*m_data)[a_dit](vof, comp);
      const Real& beta       = m_beta;
      const Real& bco        = (*m_bcoe)[a_dit](vof, comp);
      const Real& area_frac  = ebisbox.bndryArea(vof);
      const Real& weight     = m_irreg_weights[a_dit](vof, comp);
      Real flux              = weight*value*beta*bco*area_frac*a_factor;

      if(MfElectrostaticDirichletEbBc::s_areaFracWeighted){
	flux *= ebisbox.areaFracScaling(vof);
      }

      a_lphi(vof, comp) += flux;
    }
  }
}

#include <CD_NamespaceFooter.H>
