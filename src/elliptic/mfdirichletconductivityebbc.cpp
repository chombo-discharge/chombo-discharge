/*!
  @file mfdirichletconductivityebbc.cpp
  @brief Implementation of mfdiricheltconductivityebbc.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "mfdirichletconductivityebbc.H"

bool mfdirichletconductivityebbc::s_areaFracWeighted = false;
bool mfdirichletconductivityebbc::s_quadrant_based   = true;
int  mfdirichletconductivityebbc::s_lsq_radius       = 1;

#define ADJUST_STENCILS 1
#define DEBUG 0

mfdirichletconductivityebbc::mfdirichletconductivityebbc(const ProblemDomain& a_domain,
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

mfdirichletconductivityebbc::~mfdirichletconductivityebbc(){

}

LayoutData<BaseIVFAB<VoFStencil> >* mfdirichletconductivityebbc::getFluxStencil(int ivar){
  return &m_irreg_stencils;
}

void mfdirichletconductivityebbc::set_jump_object(const RefCountedPtr<jump_bc> a_jumpbc){
  m_jumpbc = a_jumpbc;
}

void mfdirichletconductivityebbc::setOrder(int a_order){
  m_order = a_order;
}

void mfdirichletconductivityebbc::define_ivs(const MFLevelGrid& a_mflg){
    
  m_ivs.define(a_mflg.get_grids());

  const DisjointBoxLayout& dbl = a_mflg.get_grids();

  int num = 0;
  m_ivs.define(dbl);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    m_ivs[dit()] = a_mflg.interface_region(dbl.get(dit()), dit());
    num += m_ivs[dit()].numPts();
  }
 
  m_definedmf = true;
}

void mfdirichletconductivityebbc::define(const LayoutData<IntVectSet>& a_cfivs, const Real& a_factor){

  if(!m_definedmf){
    MayDay::Error("mfdirichletconductivityebbc::define - must set multifluid cells first!");
  }
  else if(!m_coefSet){
    MayDay::Error("mfdirichletconductivityebbc::define - must call setCoef BEFORE calling define");
  }
  
  const int comp      = 0;
  const int num_comps = 1;

  const int otherphase = m_phase == 0 ? 1 : 0;

  const DisjointBoxLayout& dbl = m_ebisl.getDisjointLayout();

  m_irreg_weights.define(dbl);
  m_matching_weights.define(dbl);
  m_irreg_stencils.define(dbl);
  m_matching_stencils.define(dbl);

  const LayoutData<MFInterfaceFAB<VoFStencil> >& jumpstens = m_jumpbc->get_stencils();
  const LayoutData<MFInterfaceFAB<Real> >& jumpweights     = m_jumpbc->get_weights();
  const LayoutData<MFInterfaceFAB<Real> >& jump_coef       = m_jumpbc->get_bco();

  // Create stencils on irregular cells that are not part of m_ivs
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box          = dbl[dit()];
    const EBISBox& ebisbox  = m_ebisl[dit()];
    const EBGraph& ebgraph  = ebisbox.getEBGraph();
    const IntVectSet& cfivs = a_cfivs[dit()];



    IntVectSet irreg_ivs;  // All irregular cells
    IntVectSet diri_ivs;   // Pure dirichlet cells (i.e. non-matched irregular cells)
    IntVectSet match_ivs;  // Quasi-dirichlet cells (i.e. matched irregular cells)
    
    irreg_ivs  = ebisbox.getIrregIVS(box);
    irreg_ivs |= ebisbox.getMultiCells(box);

    diri_ivs = irreg_ivs - m_ivs[dit()];
    match_ivs = m_ivs[dit()];

    m_irreg_weights[dit()].define(irreg_ivs,  ebgraph, num_comps);
    m_irreg_stencils[dit()].define(irreg_ivs, ebgraph, num_comps);

    // Build stencils for pure Dirichlet type irreg cells
    for (VoFIterator vofit(irreg_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();

      Real& cur_weight        = m_irreg_weights[dit()](vof, comp);
      VoFStencil& cur_stencil = m_irreg_stencils[dit()](vof, comp);

      bool drop_order = false;

      if(m_order == 2){
	drop_order = this->get_second_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	if(drop_order){
#if DEBUG
	  pout() << "dropping order on domain = " << m_domain << " and cell = " << vof.gridIndex() << endl;
#endif
	  this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}
      }
      else if(m_order == 1){
	this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }
    }

    // Adjust stencils for matching cells
    const MFInterfaceFAB<VoFStencil>& avgStens = m_jumpbc->get_avgStencils()[dit()];
    const MFInterfaceFAB<Real>& avgWeights     = m_jumpbc->get_avgWeights()[dit()];
    const MFInterfaceFAB<Real>& avgBco         = m_jumpbc->get_avgBco()[dit()];
    for (VoFIterator vofit(match_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const VolIndex vof0(vof.gridIndex(), 0); // Multi-cell averages stored on first component

      // Everything that comes from jump_bc is from vof0, as it should!!!
      VoFStencil& cur_stencil     = m_irreg_stencils[dit()](vof, comp);
      const VoFStencil& jump_sten = avgStens.get_ivfab(m_phase)(vof0, comp);
      const Real cur_weight       = m_irreg_weights[dit()](vof, comp); 
      const Real wp               = avgWeights.get_ivfab(m_phase)(vof0, comp);
      const Real wq               = avgWeights.get_ivfab(otherphase)(vof0, comp);
      const Real bp               = avgBco.get_ivfab(m_phase)(vof0, comp);
      const Real bq               = avgBco.get_ivfab(otherphase)(vof0, comp);
      const Real factor           = cur_weight*bp/(bp*wp + bq*wq);


#if ADJUST_STENCILS
      VoFStencil addsten(jump_sten);
      addsten *= -1.0*factor;

      cur_stencil += addsten;
#endif
#if DEBUG
      const Real denom = 1./(bp*wp + bq*wq);
      if(denom      != denom     ) MayDay::Abort("mfdirichletconductivityebbc::define - got NaN denom");
      if(cur_weight != cur_weight) MayDay::Abort("mfdirichletconductivityebbc::define - got NaN cur_weight");
      if(bp         != bp        ) MayDay::Abort("mfdirichletconductivityebbc::define - got NaN bp");
      for (int i = 0; i < cur_stencil.size(); i++){
	const Real weight = cur_stencil.weight(i);
	if(weight != weight) MayDay::Abort("mfdirichletconductivityebbc::define - got NaN weight");
      }
#endif
    }

    // Scale stencils appropriately. They should be scaled by beta*bco*area_frac*a_factor
    for (VoFIterator vofit(irreg_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();
      const Real& area_frac = ebisbox.bndryArea(vof);
      const Real factor     = m_beta*(*m_bcoe)[dit()](vof, comp)*area_frac*a_factor;
      
      VoFStencil& cur_stencil = m_irreg_stencils[dit()](vof, comp);
      cur_stencil *= factor;

      if(mfdirichletconductivityebbc::s_areaFracWeighted){
	cur_stencil *= ebisbox.areaFracScaling(vof);
      }
    }
  }
}

bool mfdirichletconductivityebbc::get_second_order_sten(Real&             a_weight,
							VoFStencil&       a_stencil,
							const VolIndex&   a_vof,
							const EBISBox&    a_ebisbox,
							const IntVectSet& a_cfivs){
  CH_TIME("mfdirichletconductivityebbc::get_second_order_sten");

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

void mfdirichletconductivityebbc::get_first_order_sten(Real&             a_weight,
						       VoFStencil&       a_stencil,
						       const VolIndex&   a_vof,
						       const EBISBox&    a_ebisbox,
						       const IntVectSet& a_cfivs){
  CH_TIME("mfdirichletconductivityebbc::get_first_order_sten");

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


  if(a_stencil.size() == 0){
#if DEBUG
    pout() << "mfdirichletconductivityebbc::get_first_order sten - no sten on domain = "
	   << m_domain << "\t vof = "
	   << a_vof.gridIndex() << endl;
    MayDay::Warning("mfdirichletconductivityebbc::get_first_order_sten - could not find a stencil. ");
#endif

#if 1 // Original code. Do an approximation for the gradient at the boundary
    // Get an approximation for the cell-centered gradient
    for (int dir = 0; dir < SpaceDim; dir++){
      VoFStencil sten; 
      EBArith::getFirstDerivStencilWidthOne(sten, a_vof, a_ebisbox, dir, m_dx[dir], (IntVectSet*) &a_cfivs, 0);

      if(sten.size() > 0){
	sten *= normal[dir];
	a_stencil += sten;
      }
    }
    a_weight = 0.0;
#else
    a_stencil.clear();
    a_weight = 0.0;
#endif
  }

}

void mfdirichletconductivityebbc::applyEBFlux(EBCellFAB&                    a_lphi,
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

  const MFInterfaceFAB<Real>& inhomo = m_jumpbc->get_inhomo()[a_dit];
  const MFInterfaceFAB<Real>& homog  = m_jumpbc->get_homog()[a_dit];

  // Multi-fluid cells
  for (VoFIterator vofit(m_ivs[a_dit], ebisbox.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof   = vofit();
    const Real weight     = m_irreg_weights[a_dit](vof, comp);
    const Real beta       = m_beta;
    const Real bco        = (*m_bcoe)[a_dit](vof, comp);
    const Real area_frac  = ebisbox.bndryArea(vof);


#if ADJUST_STENCILS
    Real value = inhomo.get_ivfab(m_phase)(vof, comp); // This is the flux from the other side. Always use it. 
    if(!a_useHomogeneous){
      value += homog.get_ivfab(m_phase)(vof,comp);
    }

    // Flux - this contains the contribution from the surface charge (inhomogeneous bc) and the
    // other side of the interface (quasi-homogeneous bc). 
    Real flux = beta*weight*value*bco*area_frac*a_factor;
#else // This shouldn't work...
    Real flux = inhomo.get_ivfab(m_phase)(vof, comp);
    flux *= beta*weight*bco*area_frac*a_factor;
#endif

    if(mfdirichletconductivityebbc::s_areaFracWeighted){
      flux *= ebisbox.areaFracScaling(vof);
    }

    // Increment
    a_lphi(vof, comp) += flux;
  }

  // Pure Dirichlet irregular cells
  if(!a_useHomogeneous){
    const EBISBox& ebisbox = a_lphi.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box& box         = m_ebisl.getDisjointLayout()[a_dit];
    IntVectSet ivs   = ebisbox.getIrregIVS(box) - m_ivs[a_dit];

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      const Real& value      = (*m_data)[a_dit](vof, comp);
      const Real& beta       = m_beta;
      const Real& bco        = (*m_bcoe)[a_dit](vof, comp);
      const Real& area_frac  = ebisbox.bndryArea(vof);
      const Real& weight     = m_irreg_weights[a_dit](vof, comp);
      Real flux              = weight*value*beta*bco*area_frac*a_factor;

      if(mfdirichletconductivityebbc::s_areaFracWeighted){
      	flux *= ebisbox.areaFracScaling(vof);
      }

      a_lphi(vof, comp) += flux;
    }
  }
}
