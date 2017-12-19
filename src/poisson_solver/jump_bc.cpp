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

jump_bc::jump_bc(const MFLevelGrid&            a_mflg,
		 const LevelData<MFBaseIVFAB>& a_bco,
		 const Real&                   a_dx,
		 const int                     a_order,
		 const LayoutData<IntVectSet>* a_cfivs){
  CH_TIME("jump_bc::jump_bc(full)");

  this->define(a_mflg, a_bco, a_dx, a_order, a_cfivs);
}

jump_bc::~jump_bc(){
  CH_TIME("jump_bc::~jump_bc(full)");
}

void jump_bc::define(const MFLevelGrid&            a_mflg,
		     const LevelData<MFBaseIVFAB>& a_bco,
		     const Real&                   a_dx,
		     const int                     a_order,
		     const LayoutData<IntVectSet>* a_cfivs){
  CH_TIME("jump_bc::define");
  m_mflg   = a_mflg;
  m_dx     = a_dx;
  m_domain = m_mflg.get_domain();
  m_mfis   = m_mflg.get_mfis();
  m_grids  = m_mflg.get_grids();
  m_order  = a_order;
  m_cfivs  = a_cfivs;

  CH_assert(m_mfis->num_phases() == 2); // Only two-phase

  m_bco.define(m_grids);
  m_weights.define(m_grids);
  m_stencils.define(m_grids);

  int num = 0;
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    MFInterfaceFAB<Real>& bco         = m_bco[dit()];
    MFInterfaceFAB<Real>& weights     = m_weights[dit()];
    MFInterfaceFAB<VoFStencil>& stens = m_stencils[dit()];

    bco.define(m_mflg,     dit());
    weights.define(m_mflg, dit());
    stens.define(m_mflg,   dit());
  }

  this->set_bco(a_bco);
  this->build_stencils();

  m_defined = true;
}

void jump_bc::set_bco(const LevelData<MFBaseIVFAB>& a_bco){
  CH_TIME("jump_bc::build_stencils");

  const int comp = 0;
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    for (int iphase = 0; iphase < m_mfis->num_phases(); iphase ++){

      BaseIVFAB<Real>& bco             = m_bco[dit()].get_ivfab(iphase);
      const BaseIVFAB<Real>& bco_irreg = a_bco[dit()].get_ivfab(iphase);

      for (VoFIterator vofit(bco.getIVS(), bco.getEBGraph()); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	bco(vof, comp) = bco_irreg(vof, comp);
      }
    }
  }
}

void jump_bc::build_stencils(){
  CH_TIME("jump_bc::build_stencils");

  const int comp = 0;

  CH_assert(!m_mfis.isNull());
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
	else if(m_order == 1){
	  this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}

	if(m_order == 2 && drop_order){
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
    return true;
  }

  // If we got this far we have a stnecil
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
}

void jump_bc::get_first_order_sten(Real&             a_weight,
				   VoFStencil&       a_stencil,
				   const VolIndex&   a_vof,
				   const EBISBox&    a_ebisbox,
				   const IntVectSet& a_cfivs){
  CH_TIME("jump_bc::get_first_order_sten");

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

void jump_bc::match_bc(LevelData<BaseIVFAB<Real> >&       a_phibc,
		       const LevelData<BaseIVFAB<Real> >& a_jump,
		       const LevelData<MFCellFAB>&        a_phi,
		       const bool                         a_homogeneous){
  CH_TIME("jump_bc::match_bc(1)");

  for (DataIterator dit = a_phibc.dataIterator(); dit.ok(); ++dit){
    this->match_bc(a_phibc[dit()], a_jump[dit()], a_phi[dit()], m_bco[dit()], m_weights[dit()], m_stencils[dit()], a_homogeneous);
  }
}

void jump_bc::match_bc(LevelData<BaseIVFAB<Real> >&       a_phibc,
		       const LevelData<MFCellFAB>&        a_phi,
		       const bool                         a_homogeneous){
  CH_TIME("jump_bc::match_bc(1)");

  const int ncomp = 1;

  for (DataIterator dit = a_phibc.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real> zero(a_phibc[dit()].getIVS(), a_phibc[dit()].getEBGraph(), ncomp);
    zero.setVal(0.0);
    this->match_bc(a_phibc[dit()], zero, a_phi[dit()], m_bco[dit()], m_weights[dit()], m_stencils[dit()], a_homogeneous);
  }
}

void jump_bc::match_bc(BaseIVFAB<Real>&                  a_phibc,
		       const BaseIVFAB<Real>&            a_jump,
		       const MFCellFAB&                  a_phi,
		       const MFInterfaceFAB<Real>&       a_bco,
		       const MFInterfaceFAB<Real>&       a_weights,
		       const MFInterfaceFAB<VoFStencil>& a_stencils,
		       const bool                        a_homogeneous){
  CH_TIME("jump_bc::match_bc(2)");

#if 1 // Don't modify bc value for now. I'm working on something else
  return;
#endif


  const int comp   = 0;
  const int phase1 = 0;
  const int phase2 = 1;
  
  const IntVectSet& ivs = a_bco.get_ivs();
  
  const BaseIVFAB<Real>& bco1        = a_bco.get_ivfab(phase1);
  const BaseIVFAB<Real>& bco2        = a_bco.get_ivfab(phase2);
  const BaseIVFAB<Real>& w1          = a_weights.get_ivfab(phase1);
  const BaseIVFAB<Real>& w2          = a_weights.get_ivfab(phase2);
  
  const BaseIVFAB<VoFStencil>& sten1 = a_stencils.get_ivfab(phase1);
  const BaseIVFAB<VoFStencil>& sten2 = a_stencils.get_ivfab(phase2);
  
  const EBCellFAB& phi1              = a_phi.getPhase(phase1);
  const EBCellFAB& phi2              = a_phi.getPhase(phase2);

  const EBGraph& graph1              = bco1.getEBGraph();
  const EBGraph& graph2              = bco2.getEBGraph();

  
  // Set phibc = a_jump
  for (VoFIterator vofit(ivs, a_phibc.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit(); 
    a_phibc(vof, comp) = a_jump(vof, comp);
  }

  // First phase loop. Add first stencil stuff
  for (VoFIterator vofit(ivs, graph1); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    const VoFStencil& sten = sten1(vof, comp);
    for (int i = 0; i < sten.size(); i++){
      const VolIndex& ivof = sten.vof(i);
      const Real& iweight  = sten.weight(i);
      a_phibc(vof, comp) -= bco1(vof, comp)*phi1(ivof,comp)*iweight;
    }
  }
  
  // Second phase loop. Add second stencil stuff
  for (VoFIterator vofit(ivs, graph2); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    
    const VoFStencil& sten = sten2(vof, comp);
    for (int i = 0; i < sten.size(); i++){
      const VolIndex& ivof = sten.vof(i);
      const Real& iweight  = sten.weight(i);
      a_phibc(vof, comp) -= bco2(vof, comp)*phi2(ivof,comp)*iweight;
    }
  }
  
  // Divide by weights
  for (VoFIterator vofit(ivs, a_phibc.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const Real denom    = bco1(vof, comp)*w1(vof,comp) + bco2(vof, comp)*w2(vof, comp);
    a_phibc(vof, comp) *= 1./denom;
  }
}
