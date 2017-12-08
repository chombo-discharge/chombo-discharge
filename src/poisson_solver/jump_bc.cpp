/*!
  @file jump_bc.cpp
  @brief Implementation of jump_bc.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "jump_bc.H"

#include <EBArith.H>

bool jump_bc::s_quadrant_based = true;
int  jump_bc::s_lsq_radius     = 2;

jump_bc::jump_bc(){
  CH_TIME("jump_bc::jump_bc(weak)");
}

jump_bc::jump_bc(const MFLevelGrid& a_mflg, const Real& a_dx, const int a_order, const LayoutData<IntVectSet>* a_cfivs){
  CH_TIME("jump_bc::jump_bc(full)");

  this->define(a_mflg, a_dx, a_order, a_cfivs);
}

jump_bc::~jump_bc(){
  CH_TIME("jump_bc::~jump_bc(full)");
}

void jump_bc::define(const MFLevelGrid& a_mflg, const Real& a_dx, const int a_order, const LayoutData<IntVectSet>* a_cfivs){
  CH_TIME("jump_bc::define");
  m_mflg   = a_mflg;
  m_dx     = a_dx;
  m_domain = m_mflg.get_domain();
  m_mfis   = m_mflg.get_mfis();
  m_grids  = m_mflg.get_grids();
  m_order  = a_order;
  m_cfivs  = a_cfivs;

  m_bco.define(m_grids);
  m_soln.define(m_grids);
  m_weights.define(m_grids);
  m_stencils.define(m_grids);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    MFInterfaceFAB<Real>& bco         = m_bco[dit()];
    MFInterfaceFAB<Real>& soln        = m_soln[dit()];
    MFInterfaceFAB<Real>& weights     = m_weights[dit()];
    MFInterfaceFAB<VoFStencil>& stens = m_stencils[dit()];

    bco.define(m_mflg,     dit());
    soln.define(m_mflg,    dit());
    weights.define(m_mflg, dit());
    stens.define(m_mflg,   dit());
  }

  this->build_stencils();
}

void jump_bc::build_stencils(){
  CH_TIME("jump_bc::build_stencils");

  const int comp = 0;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    for (int iphase = 0; iphase < m_mfis->num_phases(); iphase ++){

      BaseIVFAB<Real>& weights        = m_weights[dit()].get_ivfab(iphase);
      BaseIVFAB<VoFStencil>& stencils = m_stencils[dit()].get_ivfab(iphase);

      const EBLevelGrid& eblg = m_mflg.get_eblg(iphase);
      const EBISLayout ebisl  = eblg.getEBISL();
      const EBISBox& ebisbox  = ebisl[dit()];
      const IntVectSet& cfivs = (*m_cfivs)[dit()];

      for (VoFIterator vofit(weights.getIVS(), weights.getEBGraph()); vofit.ok(); ++vofit){
	const VolIndex& vof     = vofit();
	Real& cur_weight        = weights(vof, comp);
	VoFStencil& cur_stencil = stencils(vof, comp);


	bool drop_order = false;
	
	if(m_order == 2){
	  drop_order = this->get_second_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}
	else if(m_order == 1 || drop_order){
	  this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}
      }
    }
  }
}


bool jump_bc::get_second_order_sten(Real&             a_weight,
				    VoFStencil&       a_stencil,
				    const VolIndex&   a_vof,
				    const EBISBox&    a_ebisbox,
				    const IntVectSet& a_cfivs){
  CH_TIME("jump_bc::get_second_order_sten");

  a_stencil.clear();
  bool drop_order = false;

  Vector<VoFStencil> point_stencils;
  Vector<Real> distance_along_lines;
  
  EBArith::johanStencil(drop_order, point_stencils, distance_along_lines, a_vof, a_ebisbox, m_dx*RealVect::Unit, a_cfivs);
  if(drop_order){
#if 0 // Debug
    pout() << "could not get stencil for vof = " << a_vof << "\t dx = " << m_dx << endl;
    MayDay::Abort("jump_bc::get_second_order_sten - could not get stencil");
#endif
    return true;
  }
#if 0 // Debug
  else{
    pout() << "getting stencil" << endl;
  }
#endif

  // If we got this far we have a stnecil
  CH_assert(distance_along_lines.size() >= 2);
  CH_assert(point_stencils.size() >= 2);

  Real& x1    = distance_along_lines[0];
  Real& x2    = distance_along_lines[0];
  Real denom  = x2*x2*x1 - x1*x1*x2;

  VoFStencil& phi1Sten = point_stencils[0];
  VoFStencil& phi2Sten = point_stencils[1];
  
  phi1Sten *=-x2*x2/denom;
  phi2Sten *= x1*x1/denom;

  a_weight =-x1*x1/denom + x2*x2/denom;
  a_stencil += phi1Sten;
  a_stencil += phi2Sten;
}

void jump_bc::get_first_order_sten(Real&             a_weight,
				   VoFStencil&       a_stencil,
				   const VolIndex&   a_vof,
				   const EBISBox&    a_ebisbox,
				   const IntVectSet& a_cfivs){
  CH_TIME("jump_bc::get_first_order_sten");

  MayDay::Abort("jump_bc::get_first_order_sten - not implemented");

  const RealVect& normal   = a_ebisbox.normal(a_vof);
  const RealVect& centroid = a_ebisbox.bndryCentroid(a_vof);

  if(s_quadrant_based){
    EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisbox, m_dx*RealVect::Unit, m_domain, 0);
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
}

void jump_bc::match_bc(MFInterfaceFAB<Real>& m_phibc, const MFCellFAB& m_phi){

}
