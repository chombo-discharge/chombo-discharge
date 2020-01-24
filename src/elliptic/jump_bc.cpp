/*!
  @file jump_bc.cpp
  @brief Implementation of jump_bc.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "jump_bc.H"

#include <EBArith.H>

#define USE_NEW_MATCHING 1
#define DEBUG 0

bool jump_bc::s_quadrant_based = true;
int  jump_bc::s_lsq_radius     = 1;

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
  CH_TIME("jump_bc::~jump_bc");
}

void jump_bc::get_first_order_sten(Real&             a_weight,
				   VoFStencil&       a_stencil,
				   const VolIndex&   a_vof,
				   const EBISBox&    a_ebisbox,
				   const IntVectSet& a_cfivs){
  const RealVect normal   = a_ebisbox.normal(a_vof);
  const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);

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

  if(a_stencil.size() == 0){
#if DEBUG
    pout() << "jump_bc::get_first_order sten - no sten on domain = "
	   << m_domain << "\t vof = "
	   << a_vof.gridIndex() << endl;
    MayDay::Warning("jump_bc::get_first_order_sten - could not find a stencil.");
#endif

    // Make an approximation to the cell-centered gradient
#if 1 // Original code
    for (int dir = 0; dir < SpaceDim; dir++){
      VoFStencil sten; 
      EBArith::getFirstDerivStencilWidthOne(sten, a_vof, a_ebisbox, dir, m_dx, (IntVectSet*) &a_cfivs, 0);

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

bool jump_bc::get_second_order_sten(Real&             a_weight,
				    VoFStencil&       a_stencil,
				    const VolIndex&   a_vof,
				    const EBISBox&    a_ebisbox,
				    const IntVectSet& a_cfivs){
  a_stencil.clear();
  bool drop_order = false;

  Vector<VoFStencil> point_stencils;
  Vector<Real> distance_along_lines;
  
  EBArith::johanStencil(drop_order, point_stencils, distance_along_lines, a_vof, a_ebisbox, m_dx*RealVect::Unit, a_cfivs);
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

  m_ivs.define(m_grids);
  m_bco.define(m_grids);
  m_weights.define(m_grids);
  m_stencils.define(m_grids);
  m_inhomo.define(m_grids);
  m_homog.define(m_grids);
  m_vofit_gas.define(m_grids);
  m_vofit_sol.define(m_grids);
  m_avgWeights.define(m_grids);
  m_avgStencils.define(m_grids);
  m_avgBco.define(m_grids);
  m_avgFactor.define(m_grids);
  m_avgJump.define(m_grids);

  int num = 0;

  // If we only have one phase, there are no jump cells and we don't need to do anything. 
  if(m_mfis->num_phases() > 1){
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
      MFInterfaceFAB<Real>& bco            = m_bco[dit()];
      MFInterfaceFAB<Real>& weights        = m_weights[dit()];
      MFInterfaceFAB<Real>& inhomo         = m_inhomo[dit()];
      MFInterfaceFAB<Real>& homog          = m_homog[dit()];
      MFInterfaceFAB<Real>& avgWeights     = m_avgWeights[dit()];
      MFInterfaceFAB<Real>& avgFactor      = m_avgFactor[dit()];
      MFInterfaceFAB<Real>& avgBco         = m_avgBco[dit()];
      MFInterfaceFAB<Real>& avgJump        = m_avgJump[dit()];
      MFInterfaceFAB<VoFStencil>& stens    = m_stencils[dit()];
      MFInterfaceFAB<VoFStencil>& avgStens = m_avgStencils[dit()];

      bco.define(m_mflg,        dit());
      homog.define(m_mflg,      dit());
      stens.define(m_mflg,      dit());
      inhomo.define(m_mflg,     dit());
      avgBco.define(m_mflg,     dit());
      avgJump.define(m_mflg,    dit());
      weights.define(m_mflg,    dit());
      avgStens.define(m_mflg,   dit());
      avgWeights.define(m_mflg, dit());

      m_ivs[dit()] = bco.get_ivs();
    }

    this->define_vofiter();
    this->set_bco(a_bco);
    this->build_stencils();
  }

  m_defined = true;
}

void jump_bc::define_vofiter(){
  CH_TIME("jump_bc::define_vofiter");

  const int phase1 = 0;
  const int phase2 = 1;
  
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    VoFIterator& vofit_g = m_vofit_gas[dit()];
    VoFIterator& vofit_s = m_vofit_sol[dit()];

    const IntVectSet& ivs = m_bco[dit()].get_ivs();
    const EBGraph& graph1 = m_bco[dit()].get_ivfab(phase1).getEBGraph();
    const EBGraph& graph2 = m_bco[dit()].get_ivfab(phase2).getEBGraph();

    vofit_g.define(ivs, graph1);
    vofit_s.define(ivs, graph2);
  }
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

      BaseIVFAB<Real>& bco               = m_bco[dit()].get_ivfab(iphase);
      BaseIVFAB<Real>& avgBco            = m_avgBco[dit()].get_ivfab(iphase);
      BaseIVFAB<Real>& weights           = m_weights[dit()].get_ivfab(iphase);
      BaseIVFAB<Real>& avgWeights        = m_avgWeights[dit()].get_ivfab(iphase);
      BaseIVFAB<VoFStencil>& stencils    = m_stencils[dit()].get_ivfab(iphase);
      BaseIVFAB<VoFStencil>& avgStencils = m_avgStencils[dit()].get_ivfab(iphase);

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

      // Now build the average stencils, weights, and coefficients
      const IntVectSet& ivs = avgWeights.getIVS();
      for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit){
	const IntVect iv = ivsit();
	const Vector<VolIndex> vofs = ebisbox.getVoFs(iv);

	Real totalArea = 0.0;
	Real curWeight = 0.0;
	Real curBco    = 0.0;
	int num        = 0;
	VoFStencil curStencil;
	for (int v = 0; v < vofs.size(); v++){
	  num++;
	  const VolIndex vof = vofs[v];

	  const Real area = ebisbox.bndryArea(vof);

	  totalArea += area;
	  curWeight += area*weights(vof, comp);
	  curBco    += area*bco(vof,comp);

	  VoFStencil sten = stencils(vof,comp);
	  sten *= area;
	  curStencil += sten;
	}

#if 0 // Test scale
	totalArea = totalArea/num;
#endif

	curStencil *= 1./totalArea;
	curWeight  *= 1./totalArea;
	curBco     *= 1./totalArea;



	// Put average stuff on the first vof component. No multicell storage for the average stencils
	const VolIndex vof(iv,0);
	avgWeights(vof, comp) = curWeight;
	avgStencils(vof,comp) = curStencil;
	avgBco(vof,comp)      = curBco;

#if DEBUG
	if(curWeight != curWeight) MayDay::Abort("jumpc_bc::build_stencils - got NaN weight");
	if(curBco    != curBco)    MayDay::Abort("jumpc_bc::build_stencils - got NaN bco");
#endif
      }
    }
  }
}

void jump_bc::match_bc(LevelData<BaseIVFAB<Real> >&       a_phibc,
		       const LevelData<BaseIVFAB<Real> >& a_jump,
		       const LevelData<MFCellFAB>&        a_phi,
		       const bool                         a_homogeneous){
  CH_TIME("jump_bc::match_bc(1)");

  for (DataIterator dit = a_phibc.dataIterator(); dit.ok(); ++dit){
#if USE_NEW_MATCHING
    this->new_match_bc(a_phibc[dit()], a_jump[dit()], a_phi[dit()], dit());
#else
    this->match_bc(a_phibc[dit()],
		   m_inhomo[dit()],
		   m_homog[dit()],
		   a_jump[dit()],
		   a_phi[dit()],
		   m_bco[dit()],
		   m_weights[dit()],
		   m_stencils[dit()],
		   a_homogeneous,
		   dit());
#endif
  }
}

void jump_bc::match_bc(BaseIVFAB<Real>&                  a_phibc,
		       MFInterfaceFAB<Real>&             a_inhomo,
		       MFInterfaceFAB<Real>&             a_homog,
		       const BaseIVFAB<Real>&            a_jump,
		       const MFCellFAB&                  a_phi,
		       const MFInterfaceFAB<Real>&       a_bco,
		       const MFInterfaceFAB<Real>&       a_weights,
		       const MFInterfaceFAB<VoFStencil>& a_stencils,
		       const bool                        a_homogeneous,
		       const DataIndex&                  a_dit){
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

  const EBISBox& ebisbox1            = phi1.getEBISBox();
  const EBISBox& ebisbox2            = phi2.getEBISBox();

  const EBGraph& graph1              = bco1.getEBGraph();
  const EBGraph& graph2              = bco2.getEBGraph();

  BaseIVFAB<Real>& inhomo1 = a_inhomo.get_ivfab(phase1);
  BaseIVFAB<Real>& inhomo2 = a_inhomo.get_ivfab(phase2);
  BaseIVFAB<Real>& homog1  = a_homog.get_ivfab(phase1);
  BaseIVFAB<Real>& homog2  = a_homog.get_ivfab(phase2);

  inhomo1.setVal(0.0);
  inhomo2.setVal(0.0);

  // Set phibc = 
  for (VoFIterator vofit(ivs, a_phibc.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit(); 
    a_phibc(vof, comp) = a_jump(vof, comp);
    homog1(vof, comp)  = a_jump(vof, comp);
    homog2(vof, comp)  = a_jump(vof, comp);
  }


  // First phase loop. Add first stencil stuff
  {
    VoFIterator& vofit = m_vofit_gas[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      const VoFStencil& sten = sten1(vof, comp);
      for (int i = 0; i < sten.size(); i++){
	const VolIndex& ivof = sten.vof(i);
	const Real& iweight  = sten.weight(i);
	a_phibc(vof, comp) -= bco1(vof, comp)*phi1(ivof,comp)*iweight;
	inhomo2(vof, comp) -= bco1(vof, comp)*phi1(ivof,comp)*iweight;
      }
    }
  }
  
  // Second phase loop. Add second stencil stuff
  {
    VoFIterator& vofit = m_vofit_sol[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
    
      const VoFStencil& sten = sten2(vof, comp);
      for (int i = 0; i < sten.size(); i++){
	const VolIndex& ivof = sten.vof(i);
	const Real& iweight  = sten.weight(i);
	a_phibc(vof, comp) -= bco2(vof, comp)*phi2(ivof,comp)*iweight;
	inhomo1(vof, comp) -= bco2(vof, comp)*phi2(ivof,comp)*iweight;
      }
    }
  }

  // Divide by weights
  for (VoFIterator vofit(ivs, a_phibc.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    const Real factor   = 1.0/(bco1(vof, comp)*w1(vof,comp) + bco2(vof, comp)*w2(vof, comp));
    
    a_phibc(vof, comp) *= factor;
    inhomo1(vof, comp) *= factor;
    inhomo2(vof, comp) *= factor;
    homog1(vof, comp)  *= factor;
    homog2(vof, comp)  *= factor;
  }
}


void jump_bc::new_match_bc(BaseIVFAB<Real>&                  a_phibc,
			   const BaseIVFAB<Real>&            a_jump,
			   const MFCellFAB&                  a_phi,
			   const DataIndex&                  a_dit){
  const int comp   = 0;
  const int phase1 = 0;
  const int phase2 = 1;

  // First, compute the area-weighted jump coefficient
  jump_bc::compute_avg_jump(a_jump, a_phi, a_dit);
  
  const BaseIVFAB<Real>& avg_bco1    = m_avgBco[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& avg_bco2    = m_avgBco[a_dit].get_ivfab(phase2);
  
  const BaseIVFAB<Real>& avg_w1      = m_avgWeights[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& avg_w2      = m_avgWeights[a_dit].get_ivfab(phase2);
  
  const BaseIVFAB<VoFStencil>& avg_sten1 = m_avgStencils[a_dit].get_ivfab(phase1);
  const BaseIVFAB<VoFStencil>& avg_sten2 = m_avgStencils[a_dit].get_ivfab(phase2);

  const BaseIVFAB<Real>& avg_jump1   = m_avgJump[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& avg_jump2   = m_avgJump[a_dit].get_ivfab(phase2);
  
  const EBCellFAB& phi1              = a_phi.getPhase(phase1);
  const EBCellFAB& phi2              = a_phi.getPhase(phase2);

  const EBISBox& ebisbox1            = phi1.getEBISBox();
  const EBISBox& ebisbox2            = phi2.getEBISBox();
  
  BaseIVFAB<Real>& inhomo1 = m_inhomo[a_dit].get_ivfab(phase1);
  BaseIVFAB<Real>& inhomo2 = m_inhomo[a_dit].get_ivfab(phase2);

  BaseIVFAB<Real>& homog1  = m_homog[a_dit].get_ivfab(phase1);
  BaseIVFAB<Real>& homog2  = m_homog[a_dit].get_ivfab(phase2);

  homog1.setVal(0.0); // This is the contribution of the surface charge for phase1
  homog2.setVal(0.0); // This is the contribution of the surface charge for phase2

  inhomo1.setVal(0.0); // This is the contribution from the other side for phase1
  inhomo2.setVal(0.0); // This is the contribution from the other side for phase2
  
  // Do average matching in multivalued cells.
  const IntVectSet& ivs = m_bco[a_dit].get_ivs();
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit){
    const IntVect iv = ivsit();

    const VolIndex vof0 = VolIndex(iv, 0);

    // Apply the stencils
    Real apply_sten1 = 0.0; // This result will be an EB flux contribution to phase 2
    Real apply_sten2 = 0.0; // This result will be an EB flux contribution to phase 1

    const VoFStencil& avgStencil1 = avg_sten1(vof0, 0);
    const VoFStencil& avgStencil2 = avg_sten2(vof0, 0);

    // Apply stencil on the gas side
    for (int i = 0; i < avgStencil1.size(); i++){
      const VolIndex ivof    = avgStencil1.vof(i);
      const Real     iweight = avgStencil1.weight(i);
      apply_sten1 += iweight*phi1(ivof,comp);
#if DEBUG
      const IntVect iiv = ivof.gridIndex();
      if(iweight != iweight)      MayDay::Abort("got weight2 NaN");
      if(ebisbox1.isCovered(iiv)) MayDay::Abort("reaching into covered cell for phase1");
#endif
    }

    // Apply stencil on the other side
    for (int i = 0; i < avgStencil2.size(); i++){
      const VolIndex ivof    = avgStencil2.vof(i);
      const Real     iweight = avgStencil2.weight(i);
      apply_sten2 += iweight*phi2(ivof,comp);
#if DEBUG
      const IntVect iiv = ivof.gridIndex();
      if(iweight != iweight)      MayDay::Abort("got weight2 NaN");
      if(ebisbox2.isCovered(iiv)) MayDay::Abort("reaching into covered cell for phase2");
#endif
    }


    // Compute the factor 1/(bp*wb + bq*wq)
    const Real factor = 1.0/(avg_bco1(vof0,0)*avg_w1(vof0,0) + avg_bco2(vof0,0)*avg_w2(vof0,0));

#if DEBUG 
    if(factor      != factor)      MayDay::Abort("got factor NaN");
    if(apply_sten1 != apply_sten1) MayDay::Abort("got sten1 NaN");
    if(apply_sten2 != apply_sten2) MayDay::Abort("got sten2 NaN");
#endif

    // Now make the "inhomogeneous" and "homogeneous" contributions
    const Vector<VolIndex> vofs1 = ebisbox1.getVoFs(iv);
    const Vector<VolIndex> vofs2 = ebisbox2.getVoFs(iv);

    // Phase 1 -- should the contributions be weighted somehow...?
    for (int v = 0; v < vofs1.size(); v++){
      const VolIndex vof = vofs1[v];

      homog1(vof,comp)  = avg_jump1(vof0, 0);
      inhomo1(vof,comp) = -avg_bco2(vof0,0)*apply_sten2;

      homog1(vof,comp)  *= factor;
      inhomo1(vof,comp) *= factor;

#if DEBUG
      if(homog1(vof,comp)  != homog1(vof,comp))  MayDay::Abort("got homo1 NaN");
      if(inhomo1(vof,comp) != inhomo1(vof,comp)) MayDay::Abort("got inhomo1 NaN");
#endif
    }

    // Phase 2 -- should the contributions be weighted somehow...?
    for (int v = 0; v < vofs2.size(); v++){
      const VolIndex vof = vofs2[v];

      homog2(vof,comp)  = avg_jump2(vof0, 0);
      inhomo2(vof,comp) = -avg_bco1(vof0,0)*apply_sten1;

      homog2(vof,comp)  *= factor;
      inhomo2(vof,comp) *= factor;

#if DEBUG
      if(homog2(vof,comp)  != homog2(vof,comp))  MayDay::Abort("got homo2 NaN");
      if(inhomo2(vof,comp) != inhomo2(vof,comp)) MayDay::Abort("got inhomo2 NaN");
#endif
    }
  }
}

void jump_bc::compute_avg_jump(const BaseIVFAB<Real>& a_jump, const MFCellFAB& a_phi, const DataIndex& a_dit){

  const IntVectSet ivs = m_avgJump[a_dit].get_ivs();

  const int phase1 = 0;
  const int phase2 = 1;
  
  const BaseIVFAB<Real>& bco1 = m_bco[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& bco2 = m_bco[a_dit].get_ivfab(phase2);

  BaseIVFAB<Real>& avg1 = m_avgJump[a_dit].get_ivfab(phase1);
  BaseIVFAB<Real>& avg2 = m_avgJump[a_dit].get_ivfab(phase2);
  
  const EBISBox& ebisbox1 = a_phi.getPhase(phase1).getEBISBox();
  const EBISBox& ebisbox2 = a_phi.getPhase(phase2).getEBISBox();

  
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit){
    const IntVect iv = ivsit();
    
    const Vector<VolIndex> vofs1 = ebisbox1.getVoFs(iv);
    const Vector<VolIndex> vofs2 = ebisbox2.getVoFs(iv);

    // Compute average jump condition
    Real totalArea = 0.0;
    Real totalJump = 0.0;
    for (int v = 0; v < vofs1.size(); v++){
      const VolIndex vof = vofs1[v];

      const Real area = ebisbox1.bndryArea(vof);

      totalArea += area;
      totalJump += area*a_jump(vof,0);
    }
    totalJump = totalJump/totalArea;


    // Set the jump
    const VolIndex vof0(iv, 0);
    avg1(vof0,0) = totalJump;
    avg2(vof0,0) = totalJump;
  }
}

void jump_bc::compute_dphidn(Vector<LevelData<BaseIVFAB<Real> > >&       a_dphidn,
			     const Vector<LevelData<BaseIVFAB<Real> > >& a_phibc,
			     const LevelData<MFCellFAB>&                 a_phi){
  for (int iphase = 0; iphase < m_mfis->num_phases(); iphase++){
    this->compute_dphidn(a_dphidn[iphase], a_phibc[iphase], a_phi, iphase);
  }
}

void jump_bc::compute_dphidn(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			     const LevelData<BaseIVFAB<Real> >& a_phibc,
			     const LevelData<MFCellFAB>&        a_phi,
			     const int                          a_phase){

  for (DataIterator dit = a_dphidn.dataIterator(); dit.ok(); ++dit){
    this->compute_dphidn(a_dphidn[dit()], a_phibc[dit()], a_phi[dit()].getPhase(a_phase), dit(), a_phase);
  }
}

void jump_bc::compute_dphidn(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			     const LevelData<BaseIVFAB<Real> >& a_phibc,
			     const LevelData<EBCellFAB>&        a_phi,
			     const int                          a_phase){
  
  for (DataIterator dit = a_dphidn.dataIterator(); dit.ok(); ++dit){
    this->compute_dphidn(a_dphidn[dit()], a_phibc[dit()], a_phi[dit()], dit(), a_phase);
  }
}

void jump_bc::compute_dphidn(BaseIVFAB<Real>&       a_dphidn,
			     const BaseIVFAB<Real>& a_phibc,
			     const EBCellFAB&       a_phi,
			     const DataIndex&       a_dit,
			     const int              a_phase){
  const int comp = 0.;

  const BaseIVFAB<Real>& wb             = m_weights[a_dit].get_ivfab(a_phase);
  const BaseIVFAB<VoFStencil>& stencils = m_stencils[a_dit].get_ivfab(a_phase);

  for (VoFIterator vofit(m_ivs[a_dit], a_dphidn.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    // Apply stencil wb*phi_B + sum(wi*phi_i)
    Real dphidn = wb(vof, comp)*a_phibc(vof,comp);

    const VoFStencil& sten = stencils(vof, comp);
    for (int i = 0; i < sten.size(); i++){
      const VolIndex& ivof = sten.vof(i);
      const Real& weight   = sten.weight(i);
      dphidn += weight*a_phi(ivof, comp);
    }

    a_dphidn(vof, comp) = dphidn;
  }
}

LayoutData<MFInterfaceFAB<VoFStencil> >& jump_bc::get_stencils(){
  return m_stencils;
}

LayoutData<MFInterfaceFAB<Real> >& jump_bc::get_weights(){
  return m_weights;
}

LayoutData<MFInterfaceFAB<Real> >& jump_bc::get_bco(){
  return m_bco;
}

LayoutData<MFInterfaceFAB<VoFStencil> >& jump_bc::get_avgStencils(){
#if USE_NEW_MATCHING
  return m_avgStencils;
#else
  return m_stencils;
#endif
}

LayoutData<MFInterfaceFAB<Real> >& jump_bc::get_avgWeights(){
#if USE_NEW_MATCHING
  return m_avgWeights;
#else
  return m_weights;
#endif
}

LayoutData<MFInterfaceFAB<Real> >& jump_bc::get_avgBco(){
#if USE_NEW_MATCHING
  return m_avgBco;
#else
  return m_bco;
#endif
}

LayoutData<MFInterfaceFAB<Real> >& jump_bc::get_inhomo(){
  return m_inhomo;
}

LayoutData<MFInterfaceFAB<Real> >& jump_bc::get_homog(){
  return m_homog;
}
