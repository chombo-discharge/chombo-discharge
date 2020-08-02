/*!
  @file robinconductivityebbc.cpp
  @brief Implementation of robinconductivityebbc.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "robinconductivityebbc.H"

bool robinconductivityebbc::s_quadrant_based = false;
int  robinconductivityebbc::s_lsq_radius     = 1;

robinconductivityebbc::robinconductivityebbc(const RealVect a_dx, const RealVect a_origin){
  m_dx     = a_dx;
  m_origin = a_origin;
}

robinconductivityebbc:: ~robinconductivityebbc(){
}

void robinconductivityebbc::define(const LayoutData<IntVectSet>& a_cfivs, const Real& a_factor){

  const int ncomp = 1;
  const int comp  = 0;

  const EBISLayout& ebisl      = m_eblg.getEBISL();
  const ProblemDomain& domain  = m_eblg.getDomain();
  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  m_bcstencils.define(dbl);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box             = dbl[dit()];
    const EBISBox& ebisbox     = ebisl[dit()];
    const EBGraph& ebgraph     = ebisbox.getEBGraph();
    const IntVectSet& irregIVS = ebisbox.getIrregIVS(box);
    const IntVectSet& multiIVS = ebisbox.getMultiCells(box);
    const IntVectSet& cfivs    = a_cfivs[dit()];


    BaseIVFAB<VoFStencil>& stencils = m_bcstencils[dit()];
    stencils.define(irregIVS, ebgraph, ncomp);
    
    for (VoFIterator vofit(irregIVS, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const RealVect pos  = EBArith::getVofLocation(vof, m_dx, m_origin) + ebisbox.bndryCentroid(vof)*m_dx;
      VoFStencil& stencil = stencils(vof, comp);
      
      bool found_stencil = false;
      
      if(m_type == stencil_type::taylor){
	found_stencil = this->get_taylor_sten(stencil, vof, ebisbox, cfivs);
      }
      else if(m_type == stencil_type::lsq){
	found_stencil = this->get_lsq_sten(stencil, vof, ebisbox, domain);
      }

      // Second last resort is least squares, but if we've already tried this, skip it. 
      if(!found_stencil && m_type != stencil_type::lsq){
	found_stencil = this->get_lsq_sten(stencil, vof, ebisbox, domain);
      }


      // Last resort is cell-centered value
      if(!found_stencil){
	stencil.clear();
	stencil.add(vof, 1.0);
      }


      // dphi/dn = g/b - a*phi/b. So we must scale. We only want 
      Real aco, bco;
      if(m_const_coeff){
	aco = m_aco;
	bco = m_bco;
      }
      else if(m_func_coeff){
	aco = m_robinco->aco(pos);
	bco = m_robinco->bco(pos);
      }
      else if(m_data_coeff){
	aco = (*m_acodata)[dit](vof, comp);
	bco = (*m_bcodata)[dit](vof, comp);
      }
      else{
	MayDay::Abort("robinconductivityebbc::define - must set coefficients before calling define");
      }

      stencil *= aco/bco; // Here, stencil => dphi/dn

      const Real area = ebisbox.bndryArea(vof);      
      const Real bcoe = (*m_bcoe)[dit()](vof, comp); 
      const Real fact = m_beta*bcoe*area*a_factor;   

      stencil *= fact; 
    }
  }
}

bool robinconductivityebbc::get_taylor_sten(VoFStencil&       a_stencil,
					    const VolIndex&   a_vof,
					    const EBISBox&    a_ebisbox,
					    const IntVectSet& a_cfivs){

  a_stencil.clear();
  const RealVect& dist = a_ebisbox.bndryCentroid(a_vof)*m_dx;
  IntVectSet* cfivs    = const_cast<IntVectSet*> (&a_cfivs);

  int order = EBArith::getExtrapolationStencil(a_stencil, dist, m_dx, a_vof, a_ebisbox, -1, cfivs, 0);

  return order > 0;
}

bool robinconductivityebbc::get_lsq_sten(VoFStencil&          a_stencil,
					 const VolIndex&      a_vof,
					 const EBISBox&       a_ebisbox,
					 const ProblemDomain& a_domain){
#if   CH_SPACEDIM == 2
  const int minStenSize  = 3;
#elif CH_SPACEDIM == 3
  const int minStenSize  = 7;
#endif

  bool found_stencil = false;

  Real weight;
  const RealVect dist = a_ebisbox.bndryCentroid(a_vof);
  if(s_quadrant_based){
    EBArith::getLeastSquaresGradSten(a_stencil, weight, a_vof, a_ebisbox, m_dx, a_domain, 0);
  }
  else{
    EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil,
					       weight,
					       dist,
					       RealVect::Zero,
					       a_vof,
					       a_ebisbox,
					       m_dx,
					       a_domain,
					       0,
					       s_lsq_radius);
  }
  
  if(a_stencil.size() < minStenSize){
    found_stencil = false;
  }
  else{
    found_stencil = true;
  }

  if(found_stencil){
    a_stencil.add(a_vof, weight);
    a_stencil *= m_dx[0];
    a_stencil.add(a_vof, 1.0);
  }

  return found_stencil;
}

void robinconductivityebbc::set_type(const stencil_type a_type){
  if(a_type == stencil_type::taylor || a_type == stencil_type::lsq){
    m_type = a_type;
  }
  else{
    MayDay::Abort("robinconductivityebbc::set_type - only taylor and least squares supported for now");
  }
}

void robinconductivityebbc::set_coefs(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aco = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void robinconductivityebbc::set_coefs(const RefCountedPtr<robin_coef> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void robinconductivityebbc::set_coefs(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_aco,
				      const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_bco,
				      const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_rhs){
  m_acodata = a_aco;
  m_bcodata = a_bco;
  m_rhsdata = a_rhs;

  m_const_coeff = false;
  m_func_coeff  = false;
  m_data_coeff  = true;
}

void robinconductivityebbc::applyEBFlux(EBCellFAB&                    a_lphi, 
					const EBCellFAB&              a_phi, 
					VoFIterator&                  a_vofit, 
					const LayoutData<IntVectSet>& a_cfivs, 
					const DataIndex&              a_dit, 
					const RealVect&               a_probLo, 
					const RealVect&               a_dx, 
					const Real&                   a_factor, 
					const bool&                   a_useHomogeneous, 
					const Real&                   a_time){

  const int ncomp        = 1;
  const int comp         = 0;
  
  const EBISBox& ebisbox = a_lphi.getEBISBox();

  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit){
    const VolIndex& vof = a_vofit();
    const RealVect pos  = EBArith::getVofLocation(vof, m_dx, m_origin) + ebisbox.bndryCentroid(vof)*m_dx;

    Real aco, bco, rhs;
    if(m_const_coeff){
      aco = m_aco;
      bco = m_bco;
      rhs = m_rhs;
    }
    else if(m_func_coeff){
      aco = m_robinco->aco(pos);
      bco = m_robinco->bco(pos);
      rhs = m_robinco->rhs(pos);
    }
    else if(m_data_coeff){
      aco = (*m_acodata)[a_dit](vof, comp);
      bco = (*m_bcodata)[a_dit](vof, comp);
      rhs = (*m_rhsdata)[a_dit](vof, comp);
    }

    Real flux = rhs/bco;

    const Real area = ebisbox.bndryArea(vof);
    const Real bcoe = (*m_bcoe)[a_dit](vof, comp); // bco from Helmholtz equation, not the one from Robin bc

    flux *= m_beta*bcoe*area;

    a_lphi(vof, comp) += flux*a_factor;
  }
}

LayoutData<BaseIVFAB<VoFStencil> >* robinconductivityebbc::getFluxStencil(int ivar){
  return &m_bcstencils;
}
