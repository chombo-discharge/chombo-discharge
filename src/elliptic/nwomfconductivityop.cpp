/*!
  @file nwomfconductivityop.cpp
  @brief Implementation of nwomfconductivityop.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "mfdirichletconductivityebbc.H"
#include "nwomfconductivityop.H"
#include "mfalias.H"

#include <EBArith.H>
#include <DirichletConductivityDomainBC.H>
#include <MFLevelDataOps.H>
#include <BaseIVFactory.H>
#include <EBAMRDataOps.H>

#define verb 0

nwomfconductivityop::nwomfconductivityop(){
  
}
  

nwomfconductivityop::~nwomfconductivityop(){
  for (int i = 0; i < m_alias.size(); i++){
    m_alias[i] = NULL;
  }
}

void nwomfconductivityop::define(const RefCountedPtr<mfis>&                    a_mfis,
			      const RefCountedPtr<BaseDomainBCFactory>&     a_dombc,
			      const RefCountedPtr<LevelData<MFCellFAB> >&   a_aco,
			      const RefCountedPtr<LevelData<MFFluxFAB> >&   a_bco,
			      const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_bco_irreg,
			      const NWOMFQuadCFInterp&                      a_quadcfi,
			      const MFLevelGrid&                            a_mflg_fine,
			      const MFLevelGrid&                            a_mflg,
			      const MFLevelGrid&                            a_mflg_coar,
			      const MFLevelGrid&                            a_mflg_coar_mg,
			      const ProblemDomain&                          a_domain,
			      const bool&                                   a_layout_changed,
			      const bool&                                   a_has_mg,
			      const bool&                                   a_has_fine,
			      const bool&                                   a_has_coar,
			      const int&                                    a_ref_to_fine,
			      const int&                                    a_ref_to_coar,
			      const int&                                    a_relax_type,
			      const int&                                    a_order_ebbc,
			      const IntVect&                                a_ghost_phi,
			      const IntVect&                                a_ghost_rhs,
			      const Real&                                   a_dx,
			      const Real&                                   a_dx_coar,
			      const Real&                                   a_alpha,
			      const Real&                                   a_beta,
			      const RealVect&                               a_origin){


  const int num_phases = a_mfis->num_phases();
  
  const int num_alias  = 6;

  m_mfis = a_mfis;
  m_ncomp = 1;
  m_relax = a_relax_type;
  m_domain = a_domain;
  m_phases = num_phases;
  m_ebops.resize(num_phases);
  m_ebbc.resize(num_phases);
  m_acoeffs.resize(num_phases);
  m_bcoeffs.resize(num_phases);
  m_bcoeffs_irr.resize(num_phases);
  m_alias.resize(num_alias);
  m_ref_to_coarser = a_ref_to_coar;
  m_ghost_phi = a_ghost_phi;
  m_ghost_rhs = a_ghost_rhs;
  m_origin = a_origin;
  m_dx = a_dx;
  m_dirival.resize(num_phases);
  m_multifluid = m_mfis->num_phases() > 1;

  if(a_has_mg){
    m_mflg_coar_mg = a_mflg_coar_mg;
  }

  for (int i = 0; i < num_alias; i++){
    m_alias[i] = new LevelData<EBCellFAB>();
  }

  // Object for matching boundary conditions. Native EBBC is data-based Dirichlet
  m_jumpbc = RefCountedPtr<jump_bc> (new jump_bc(a_mflg, *a_bco_irreg, a_dx, a_order_ebbc, (a_mflg.get_eblg(0)).getCFIVS()));

  Vector<EBISLayout> layouts(num_phases);
  Vector<int> comps(num_phases);;
  
  for (int iphase = 0; iphase < num_phases; iphase++){

    const EBLevelGrid& eblg      = a_mflg.get_eblg(iphase);
    const EBISLayout&  ebisl     = eblg.getEBISL();

    layouts[iphase] = ebisl;
    comps[iphase]   = m_ncomp;


    EBLevelGrid eblg_fine;
    EBLevelGrid eblg_coar;
    EBLevelGrid eblg_mg;

    if(a_has_fine){
      eblg_fine = a_mflg_fine.get_eblg(iphase);
    }
    if(a_has_coar){
      eblg_coar = a_mflg_coar.get_eblg(iphase);
    }
    if(a_has_mg){
      eblg_mg   = a_mflg_coar_mg.get_eblg(iphase);
    }


    RefCountedPtr<NWOEBQuadCFInterp> quadcfi;
    if(a_has_coar){
      quadcfi = a_quadcfi.get_quadcfi_ptr(iphase);
      CH_assert(!quadcfi.isNull());
    }

    // Coefficients
    m_acoeffs[iphase]     = RefCountedPtr<LevelData<EBCellFAB> >        (new LevelData<EBCellFAB>());
    m_bcoeffs[iphase]     = RefCountedPtr<LevelData<EBFluxFAB> >        (new LevelData<EBFluxFAB>());
    m_bcoeffs_irr[iphase] = RefCountedPtr<LevelData<BaseIVFAB<Real> > > (new LevelData<BaseIVFAB<Real> >());

    
    // Domain BC
    ConductivityBaseDomainBC* bc = (ConductivityBaseDomainBC*) a_dombc->create(a_domain, ebisl, a_dx*RealVect::Unit);
    RefCountedPtr<ConductivityBaseDomainBC> dbc(bc);

    // EB BC
    m_ebbc[iphase] = RefCountedPtr<mfdirichletconductivityebbc> (new mfdirichletconductivityebbc(a_domain,
												 ebisl,
												 a_dx*RealVect::Unit,
												 &a_ghost_phi,
												 &a_ghost_rhs,
												 iphase)),
    m_ebbc[iphase]->set_jump_object(m_jumpbc);
    m_ebbc[iphase]->setOrder(a_order_ebbc);
    m_ebbc[iphase]->define_ivs(a_mflg);


    // Create storage for data-based dirichlet boundary conditions
    LayoutData<IntVectSet> ivs(eblg.getDBL());
    for (DataIterator dit = ivs.dataIterator(); dit.ok(); ++dit){
      ivs[dit()] = ebisl[dit()].getIrregIVS(eblg.getDBL().get(dit()));
    }

    BaseIVFactory<Real> ivfact(ebisl, ivs);
    m_dirival[iphase] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> > (eblg.getDBL(), 1, IntVect::Zero, ivfact));

    EBLevelDataOps::setVal(*m_dirival[iphase], 0.0);
    m_ebbc[iphase]->setData(m_dirival[iphase]);

    mfalias::aliasMF(*m_acoeffs[iphase],     iphase, *a_aco);
    mfalias::aliasMF(*m_bcoeffs[iphase],     iphase, *a_bco);
    mfalias::aliasMF(*m_bcoeffs_irr[iphase], iphase, *a_bco_irreg);

    const RefCountedPtr<LevelData<EBCellFAB> >&        aco     = m_acoeffs[iphase];
    const RefCountedPtr<LevelData<EBFluxFAB> >&        bco     = m_bcoeffs[iphase];
    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& bco_irr = m_bcoeffs_irr[iphase];
    const RefCountedPtr<mfdirichletconductivityebbc>&  ebbc    = m_ebbc[iphase];

    const Real alpha = a_alpha;
    const Real beta  = a_beta;

#if verb
    pout() << "nwomfconductivityop::creating oper" << endl;
#endif
#if 0
    m_ebops[iphase] = RefCountedPtr<nwoebconductivityop> (new nwoebconductivityop(eblg_fine,
									    eblg,
									    eblg_coar,
									    eblg_mg,
									    quadcfi,
									    dbc,
									    ebbc,
									    a_dx,
									    a_dx_coar,
									    a_ref_to_fine,
									    a_ref_to_coar,
									    a_has_fine,
									    a_has_coar,
									    a_has_mg,
									    a_layout_changed,
									    alpha,
									    beta,
									    aco,
									    bco,
									    bco_irr,
									    a_ghost_phi,
									    a_ghost_rhs,
									    a_relax_type));
#endif
#if verb
    pout() << "nwomfconductivityop::done creating oper" << endl;
#endif
  }

  MFCellFactory* factory = new MFCellFactory(layouts, comps);
  RefCountedPtr<DataFactory<MFCellFAB> > fac(factory);
  m_ops.define(fac);

  EBArith::getMultiColors(m_colors);
}

void nwomfconductivityop::set_jump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump){
  m_jump = a_jump;
}

void nwomfconductivityop::set_electrodes(const Vector<electrode>& a_electrodes, const RefCountedPtr<BaseBCFuncEval> a_potential){
#if verb
  pout() << "nwomfconductivityop::set_electrodes"<< endl;
#endif

  m_electrodes = a_electrodes;
  m_potential  = a_potential;
  this->set_bc_from_levelset();
}

void nwomfconductivityop::update_bc(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
  CH_TIME("nwomfconductivityop::update_bc");
#if verb
  pout() << "nwomfconductivityop::update_bc"<< endl;
#endif

  LevelData<MFCellFAB>* phi = const_cast<LevelData<MFCellFAB>* > (&a_phi);
  phi->exchange();

  this->set_bc_from_levelset();
  this->set_bc_from_matching(a_phi, a_homogeneous);

#if verb
  pout() << "nwomfconductivityop::update_bc - done" << endl;
#endif
}

void nwomfconductivityop::set_bc_from_levelset(){
#if verb
  pout() << "nwomfconductivityop::set_bc_from_levelset"<< endl;
#endif

  const int comp = 0;
  
  for (int iphase = 0; iphase < m_phases; iphase++){
    LevelData<BaseIVFAB<Real> >& val = *m_dirival[iphase];

    for (DataIterator dit = val.dataIterator(); dit.ok(); ++dit){
      const IntVectSet& ivs  = val[dit()].getIVS();
      const EBGraph& ebgraph = val[dit()].getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect& pos = EBArith::getVofLocation(vof, m_dx, m_origin);

	if(m_electrodes.size() > 0){
	  int  func = 0;
	  Real dist = m_electrodes[0].get_function()->value(pos);

	  for (int i = 1; i < m_electrodes.size(); i++){
	    Real cur_val = (m_electrodes[i].get_function())->value(pos);
	    if(Abs(cur_val) < dist){
	      func = i;
	      dist  = cur_val;
	    }
	  }

	  Real potential;
	  if(m_electrodes[func].is_live()){
	    potential = m_potential->value(pos, comp);
	  }
	  else{
	    potential = 0.;
	  }
	  val[dit()](vof, comp) = potential;
	}
      }
    }
  }
}

void nwomfconductivityop::set_bc_from_matching(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
  CH_TIME("nwomfconductivityop::set_bc_from_matching");
#if verb
  pout() << "nwomfconductivityop::set_bc_from_matching"<< endl;
#endif

  if(m_multifluid){
    for (int iphase = 0; iphase < m_phases; iphase++){
      m_jumpbc->match_bc(*m_dirival[iphase], *m_jump, a_phi, a_homogeneous); 
    }
  }

#if verb
  pout() << "nwomfconductivityop::set_bc_from_matching - done"<< endl;
#endif
}

void nwomfconductivityop::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta){
  CH_TIME("nwomfconductivityop::setAlphaAndBeta");
#if verb
  pout() << "nwomfconductivityop::setalaphandbeta"<< endl;
#endif
  for (int iphase = 0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setAlphaAndBeta(a_alpha, a_beta);
  }
}

void nwomfconductivityop::diagonalScale(LevelData<MFCellFAB>& a_rhs){
  CH_TIME("nwomfconductivityop::diagonalScale");
#if verb
  pout() << "nwomfconductivityop::setdiagonalscale"<< endl;
#endif
  MayDay::Abort("nwomfconductivityop::DiagonalScale - where did i get called?");

  // // Operator diagonal scale
  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->diagonalScale(*m_alias[0], true);
  }

  MayDay::Abort("nwomfconductivityop::diagonalScale - add scaling for area fractions");
  //  MFLevelDataOps::kappaWeight(a_rhs); // This came from MFPoissonOp
}

void nwomfconductivityop::divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "nwomfconductivityop::dividebyidenticyoef"<< endl;
#endif

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->divideByIdentityCoef(*m_alias[0]);
  }
}

void nwomfconductivityop::applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
					 const LevelData<MFCellFAB>& a_phi){
#if verb
  pout() << "nwomfconductivityop::applyopnoboundary"<< endl;
#endif
  CH_TIME("nwomfconductivityop::applyOpNoBoundary");

  this->update_bc(a_phi, true);
  for (int iphase=0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_opPhi);
    mfalias::aliasMF(*m_alias[1], iphase, a_phi);
    m_ebops[iphase]->applyOpNoBoundary(*m_alias[0], *m_alias[1]);
  }
}

void nwomfconductivityop::setTime(Real a_oldTime, Real a_mu, Real a_dt){
  CH_TIME("nwomfconductivityop::setTime");
#if verb
  pout() << "nwomfconductivityop::settime"<< endl;
#endif
  for (int iphase=0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setTime(a_oldTime, a_mu, a_dt);
  }
}

void nwomfconductivityop::residual(LevelData<MFCellFAB>&        a_lhs,
				const LevelData<MFCellFAB>&  a_phi,
				const LevelData<MFCellFAB>&  a_rhs,
				bool                         a_homogeneous){
  CH_TIME("nwomfconductivityop::residual");
#if verb
  pout() << "nwomfconductivityop::residual"<< endl;
#endif
  
  this->applyOp(a_lhs, a_phi, a_homogeneous);

  // Make residual rhs - L(phi)
  this->scale(a_lhs,      -1.0);
  this->incr(a_lhs, a_rhs, 1.0);

}

void nwomfconductivityop::preCond(LevelData<MFCellFAB>&       a_correction,
			       const LevelData<MFCellFAB>& a_residual){
  CH_TIME("nwomfconductivityop::preCond");
#if verb
  pout() << "nwomfconductivityop::precond"<< endl;
#endif
  this->relax(a_correction, a_residual, 40);
}

void nwomfconductivityop::applyOp(LevelData<MFCellFAB>&        a_lhs,
			       const LevelData<MFCellFAB>&  a_phi,
			       bool                         a_homogeneous){
  CH_TIME("nwomfconductivityop::applyOp");
#if verb
  pout() << "nwomfconductivityop::applyop"<< endl;
#endif

  this->update_bc(a_phi, a_homogeneous);

  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_lhs);
    mfalias::aliasMF(*m_alias[1], i, a_phi);
    
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneous, true);
  }
}

void nwomfconductivityop::create(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("nwomfconductivityop::create");
#if verb
  pout() << "nwomfconductivityop::create"<< endl;
#endif
  m_ops.create(a_lhs, a_rhs);
}

void nwomfconductivityop::createCoarsened(LevelData<MFCellFAB>&       a_lhs,
				       const LevelData<MFCellFAB>& a_rhs,
				       const int&                  a_refRat) {
  CH_TIME("nwomfconductivityop::createCoarsened");
#if verb
  pout() << "nwomfconductivityop::createCoarsened"<< endl;
#endif

  const IntVect ghostVec = a_rhs.ghostVect();

  Vector<EBLevelGrid> eblg(m_phases);
  Vector<EBISLayout> ebisl(m_phases);
  for (int i = 0; i < m_phases; i++){
    eblg[i] = m_ebops[i]->getEBLG();
  }

  DisjointBoxLayout dbl_coar;
  ProblemDomain     dom_coar = coarsen(eblg[0].getDomain(), a_refRat);
  coarsen(dbl_coar, eblg[0].getDBL(), a_refRat);

  for (int i = 0; i < m_phases; i++){
    eblg[i].getEBIS()->fillEBISLayout(ebisl[i], dbl_coar, dom_coar, a_rhs.ghostVect()[0]);
    if(a_refRat > 2){
      ebisl[i].setMaxRefinementRatio(a_refRat, eblg[i].getEBIS());
    }
  }

  
  Vector<int> vncomp(m_phases, a_rhs.nComp());
  MFCellFactory fact(ebisl, vncomp);
  a_lhs.define(dbl_coar, a_rhs.nComp(), a_rhs.ghostVect(), fact);

#if verb
  pout() << "nwomfconductivityop::createCoarsened - done"<< endl;
#endif
}

void nwomfconductivityop::assign(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "nwomfconductivityop::assign"<< endl;
#endif
  m_ops.assign(a_lhs, a_rhs);
#if verb
  pout() << "nwomfconductivityop::assign - doen"<< endl;
#endif
}


Real nwomfconductivityop::dotProduct(const LevelData<MFCellFAB>& a_data1,
				  const LevelData<MFCellFAB>& a_data2){
#if verb
  pout() << "nwomfconductivityop::dotproduct"<< endl;
#endif
  Real accum = 0.0;

  Real volume = 0.0;

  for (DataIterator dit = a_data1.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      Real phaseVolume;
      for (int i=0; i<m_phases; i++)
	{
	  const EBCellFAB& data1 = a_data1[d].getPhase(i);
	  const EBCellFAB& data2 = a_data2[d].getPhase(i);
	  const Box& box = a_data1.getBoxes().get(d);


	  accum += EBLevelDataOps::sumKappaDotProduct(phaseVolume,data1,data2,box,
						      EBLEVELDATAOPS_ALLVOFS,m_domain);

	  volume += phaseVolume;
	}
    }



#ifdef CH_MPI
  Real recv;
  int result;

  result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL,
			 MPI_SUM, Chombo_MPI::comm);
  accum = recv;

  result = MPI_Allreduce(&volume, &recv, 1, MPI_CH_REAL,
			 MPI_SUM, Chombo_MPI::comm);
  volume = recv;
#endif

  if (volume > 0.0)
    {
      accum = accum / volume;
    }
  else
    {
      accum = 0.0;
    }

  return accum;
}

void nwomfconductivityop::incr(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_rhs,
			    Real                        a_scale){
#if verb
  pout() << "nwomfconductivityop::incr"<< endl;
#endif
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
#if verb
  pout() << "nwomfconductivityop::incr - done"<< endl;
#endif
}

void nwomfconductivityop::axby(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_x,
			    const LevelData<MFCellFAB>& a_y,
			    Real a,
			    Real b){
#if verb
  pout() << "nwomfconductivityop::axby"<< endl;
#endif

  m_ops.axby(a_lhs, a_x, a_y, a, b);
}

void nwomfconductivityop::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale){
  CH_TIME("nwomfconductivityop::scale");
#if verb
  pout() << "nwomfconductivityop::scale"<< endl;
#endif
  m_ops.scale(a_lhs, a_scale);
#if verb
  pout() << "nwomfconductivityop::scale - doen"<< endl;
#endif
}

Real nwomfconductivityop::norm(const LevelData<MFCellFAB>& a_x, int a_ord){
  CH_TIME("nwomfconductivityop::norm");
#if verb
  pout() << "nwomfconductivityop::norm"<< endl;
#endif
  this->update_bc(a_x, true);
  Real volume;
  Real rtn = this->kappaNorm(volume, a_x, a_ord);
#if verb
  pout() << "nwomfconductivityop::norm - done"<< endl;
#endif
  return rtn;
}

Real nwomfconductivityop::kappaNorm(Real&                       a_volume,
				 const LevelData<MFCellFAB>& a_data,
				 int                         a_p) const {

#if verb
  pout() << "nwomfconductivityop::kappaNorm"<< endl;
#endif
  Real accum = 0.0;

  a_volume = 0.0;
  int ncomp=a_data.nComp();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      Vector<Real> cur(ncomp,0.0);
      Real curVolume = 0.0;
      Vector<Real> phaseCur(ncomp,0.0);
      Real phaseVolume;
      for (int i=0; i<m_phases; i++)
	{
	  const EBCellFAB& data = a_data[d].getPhase(i);
	  const Box& box = a_data.getBoxes().get(d);
	  phaseCur = EBLevelDataOps::vectorSumKappaPow(phaseVolume,data,box,
						       EBLEVELDATAOPS_ALLVOFS,m_domain,a_p);
	  for (int i=0; i<ncomp; i++){
	    if (a_p == 0){
	      if ( phaseCur[i] > cur[i]){
		cur[i] = phaseCur[i];
	      }
	    }
	    else{
	      cur[i]+= phaseCur[i];
	    }
	  }
	  curVolume += phaseVolume;
	}
      a_volume += curVolume;

      for (int i=0; i<ncomp; i++) {
	if (a_p == 0){
	  if (cur[i] > accum){
	    accum = cur[i];
	  }
	}
	else {
	  accum += cur[i];
	}
      }
    }
#ifdef CH_MPI
  if (a_p == 0)
    {
      Real recv;
      int result;

      result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL,
			     MPI_MAX, Chombo_MPI::comm);
      accum = recv;
    }
  else
    {
      Real recv;
      int result;

      result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL,
			     MPI_SUM, Chombo_MPI::comm);
      accum = recv;

      result = MPI_Allreduce(&a_volume, &recv, 1, MPI_CH_REAL,
			     MPI_SUM, Chombo_MPI::comm);
      a_volume = recv;
    }
#endif
  if (a_p != 0)
    {
      if (a_volume > 0.0)
	{
	  accum = accum / a_volume;
	}
      else
	{
	  accum = 0;
	}
      Real invPow = 1.0/a_p;
      accum = pow(accum,invPow);
    }

  return accum;
}

void nwomfconductivityop::setToZero(LevelData<MFCellFAB>& a_x){
  CH_TIME("nwomfconductivityop::setToZero");
#if verb
  pout() << "nwomfconductivityop::setToZero"<< endl;
#endif
  m_ops.setToZero(a_x);
#if verb
  pout() << "nwomfconductivityop::setToZero - done"<< endl;
#endif
}

void nwomfconductivityop::relax(LevelData<MFCellFAB>&       a_e,
			     const LevelData<MFCellFAB>& a_residual,
			     int                         a_iterations){
  CH_TIME("nwomfconductivityop::relax");
#if verb
  pout() << "nwomfconductivityop::relax"<< endl;
#endif
#if 0
// #if 0
//   for (int i=0; i < a_iterations; i++){
//     if (m_relax == 0){
//       this->levelJacobi(a_e, a_residual);
//     }
//     else if (m_relax == 1){
//       this->levelMulticolorGS(a_e, a_residual);
//     }
//     else if (m_relax == 2){
//       this->levelGSRB(a_e, a_residual);
//     }
//     else{
//       MayDay::Error("nwomfconductivityop: Invalid relaxation type");
//     }
//   }
// #else
//   const bool homogeneous = true;
//   for (int i = 0; i < a_iterations; i++){
//     this->update_bc(a_e, homogeneous);

//     for (int iphase = 0; iphase < m_phases; iphase++){
//       mfalias::aliasMF(*m_alias[0], iphase, a_e);
//       mfalias::aliasMF(*m_alias[1], iphase, a_residual);

//       m_ebops[iphase]->lazyGauSai(*m_alias[0], *m_alias[1]);
//     }
//   }
// #endif

  const bool homogeneous = true;
  
  if(!m_multifluid){ // Single-fluid code
    for (int iphase = 0; iphase < m_phases; iphase++){
      mfalias::aliasMF(*m_alias[0], iphase, a_e);
      mfalias::aliasMF(*m_alias[1], iphase, a_residual);

      for (int i = 0; i < a_iterations; i++){
	this->update_bc(a_e, homogeneous);
	m_ebops[iphase]->relax(*m_alias[0], *m_alias[1], 1);
      }
    }
  }
  else { // Multifluid code
#if 0 // Original code
    for (int i = 0; i < a_iterations; i++){
      this->update_bc(a_e, homogeneous);

      for (int iphase = 0; iphase < m_phases; iphase++){
	mfalias::aliasMF(*m_alias[0], iphase, a_e);
	mfalias::aliasMF(*m_alias[1], iphase, a_residual);

	m_ebops[iphase]->lazyGauSai(*m_alias[0], *m_alias[1]);
      }
    }
#else // Place to put optimized code
    const DisjointBoxLayout& dbl = a_e.disjointBoxLayout();
    
    LevelData<MFCellFAB> lphi;
    this->create(lphi, a_residual);
    
    const bool homogeneous = true;

    mfalias::aliasMF(*m_alias[0], 0, a_e);        
    mfalias::aliasMF(*m_alias[1], 0, a_residual);
    mfalias::aliasMF(*m_alias[2], 0, lphi);
    mfalias::aliasMF(*m_alias[3], 1, a_e);
    mfalias::aliasMF(*m_alias[4], 1, a_residual);
    mfalias::aliasMF(*m_alias[5], 1, lphi);

    
    for (int i = 0; i < a_iterations; i++){
      this->update_bc(a_e, homogeneous);

      for (int icolor = 0; icolor < m_colors.size(); icolor++){

	// Apply homogeneous CFBCs

	  // Get coarse-fine boundary conditions
	  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	    for (int dir = 0; dir < SpaceDim; dir++){
	      for (SideIterator sit; sit.ok(); sit.next()){
		m_ebops[0]->applyHomogeneousCFBCs((*m_alias[0])[dit()], dit(), dir, sit());
		m_ebops[1]->applyHomogeneousCFBCs((*m_alias[3])[dit()], dit(), dir, sit());
	      }
	    }
	  }

	  m_ebops[0]->applyOp(*m_alias[2], *m_alias[0], NULL, true, true, false);
	  m_ebops[1]->applyOp(*m_alias[5], *m_alias[3], NULL, true, true, false);


	  // Don't relax covered boxes. This implies that only intersection boxes will be relaxed twice. Can't avoid this. 
	  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	    const EBISBox& ebisbox0 = (*m_alias[0])[dit()].getEBISBox();
	    if(!ebisbox0.isAllCovered()){
	      m_ebops[0]->gsrbColor((*m_alias[0])[dit()],
				    (*m_alias[2])[dit()],
				    (*m_alias[1])[dit()],
				    dbl.get(dit()),
				    dit(),
				    m_colors[icolor]);
	    }

	    // Check if box is covered. Relax if its not
	    const EBISBox& ebisbox1 = (*m_alias[3])[dit()].getEBISBox();
	    if(!ebisbox1.isAllCovered()){
	      m_ebops[1]->gsrbColor((*m_alias[3])[dit()],
				    (*m_alias[5])[dit()],
				    (*m_alias[4])[dit()],
				    dbl.get(dit()),
				    dit(),
				    m_colors[icolor]);
	    }
	  }

	  if((icolor-1) % 2 == 0 && icolor - 1 < m_colors.size()){
	    m_alias[0]->exchange();
	    m_alias[3]->exchange();
	  }
	}
      }
#endif
  }
#endif
}

void nwomfconductivityop::levelJacobi(LevelData<MFCellFAB>&       a_phi,
				   const LevelData<MFCellFAB>& a_rhs,
				   const int                   a_iterations){
  CH_TIME("nwomfconductivityop::levelJacobi");
  bool homogeneous = true;

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_phi);
    mfalias::aliasMF(*m_alias[1], iphase, a_rhs);

#if 1 // This is the only way we can make it converge for now
    // m_ebops[iphase]->relaxPoiJac(*m_alias[0], *m_alias[1], 1);
    //    m_ebops[iphase]->relaxGauSai(*m_alias[0], *m_alias[1], 1);
    //m_ebops[iphase]->lazyGauSai(*m_alias[0], *m_alias[1]);
    /// m_ebops[iphase]->relaxGSRBFast(*m_alias[0], *m_alias[1], 1);
#else
    m_ebops[iphase]->relax_mf(*m_alias[0], *m_alias[1], 1);
#endif
  }
}

void nwomfconductivityop::createCoarser(LevelData<MFCellFAB>&       a_coarse,
				     const LevelData<MFCellFAB>& a_fine,
				     bool                        ghosted){
  CH_TIME("nwomfconductivityop::createCoarser");
#if verb
  pout() << "nwomfconductivityop::createCoarser"<< endl;
#endif

  DisjointBoxLayout dbl = m_mflg_coar_mg.get_grids();

  Vector<int> comps(m_phases, m_ncomp);
  Vector<EBISLayout> ebisl(m_phases);
  for (int iphase = 0; iphase < m_phases; iphase++){
    ebisl[iphase] = m_mflg_coar_mg.get_eblg(iphase).getEBISL();
  }

  MFCellFactory factory(ebisl, comps);
  a_coarse.define(dbl, m_ncomp, a_fine.ghostVect(), factory);
}

void nwomfconductivityop::restrictResidual(LevelData<MFCellFAB>&       a_resCoarse,
					LevelData<MFCellFAB>&       a_phiFine,
					const LevelData<MFCellFAB>& a_rhsFine){
  CH_TIME("nwomfconductivityop::restrictResidual");
#if verb
  pout() << "nwomfconductivityop::restrictResidual"<< endl;
#endif
  for (int i=0; i < m_phases; i++) {
    mfalias::aliasMF(*m_alias[0], i, a_resCoarse);
    mfalias::aliasMF(*m_alias[1], i, a_phiFine);
    mfalias::aliasMF(*m_alias[2], i, a_rhsFine);
    
    m_ebops[i]->restrictResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void nwomfconductivityop::prolongIncrement(LevelData<MFCellFAB>&       a_phiThisLevel,
					const LevelData<MFCellFAB>& a_correctCoarse){
  CH_TIME("nwomfconductivityop::prolongIncrement");
#if verb
  pout() << "nwomfconductivityop::prolongIncrement"<< endl;
#endif
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_phiThisLevel);
    mfalias::aliasMF(*m_alias[1], i, a_correctCoarse);
    
    m_ebops[i]->prolongIncrement(*m_alias[0], *m_alias[1]);
  }
}

void nwomfconductivityop::AMRResidual(LevelData<MFCellFAB>&       a_residual,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   const LevelData<MFCellFAB>& a_rhs,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("nwomfconductivityop::AMRResidual");
#if verb
  pout() << "nwomfconductivityop::amrresidual"<< endl;
#endif

  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void nwomfconductivityop::AMROperator(LevelData<MFCellFAB>&       a_LofPhi,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("nwomfconductivityop::AMROperator");
#if verb
  pout() << "nwomfconductivityop::amroperator" << endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_LofPhi);
    mfalias::aliasMF(*m_alias[1], iphase, a_phiFine);
    mfalias::aliasMF(*m_alias[2], iphase, a_phi);
    mfalias::aliasMF(*m_alias[3], iphase, a_phiCoarse);

    nwomfconductivityop* finerOp = (nwomfconductivityop*) a_finerOp;
#if verb
    pout() << "nwomfconductivityop::AMROperator - apply nwoebconductivityops" << endl;
#endif
    m_ebops[iphase]->AMROperator(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_homogeneousBC, finerOp->m_ebops[iphase]);
    // m_ebops[iphase]->applyOp(*m_alias[0], *m_alias[1], m_alias[2], a_homogeneousBC, false);
    // m_ebops[iphase]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[iphase]);
#if verb
    pout() << "nwomfconductivityop::AMROperator - apply nwoebconductivityops - done" << endl;
#endif
  }

										 


#if verb
  pout() << "nwomfconductivityop::amroperator - done" << endl;
#endif
}

void nwomfconductivityop::AMRResidualNC(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("nwomfconductivityop::AMRResidualNC");
#if verb
  pout() << "nwomfconductivityop::amrresidualnc"<< endl;
#endif

  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);

#if verb
  pout() << "nwomfconductivityop::amrresidualnc - done"<< endl;
#endif
}

void nwomfconductivityop::AMROperatorNC(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("nwomfconductivityop::AMROperatorNC");
#if verb
  pout() << "nwomfconductivityop::amroperatornc"<< endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_LofPhi);
    mfalias::aliasMF(*m_alias[1], i, a_phiFine);
    mfalias::aliasMF(*m_alias[2], i, a_phi);

    nwomfconductivityop* finerOp = (nwomfconductivityop*) a_finerOp;
    m_ebops[i]->AMROperatorNC(*m_alias[0], *m_alias[1], *m_alias[2], a_homogeneousBC, finerOp->m_ebops[i]);

    // m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneousBC, true);

    // nwomfconductivityop* finerOp = (nwomfconductivityop*) a_finerOp;
    // m_ebops[i]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[i]);
  }
#if verb
  pout() << "nwomfconductivityop::amroperatornc - done"<< endl;
#endif
}


void nwomfconductivityop::AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC){
  CH_TIME("nwomfconductivityop::AMRResidualNF");
#if verb
  pout() << "nwomfconductivityop::amrresidualnf"<< endl;
#endif

  this->AMROperatorNF(a_residual, a_phi, a_phiCoarse, a_homogeneousBC);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void nwomfconductivityop::AMROperatorNF(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     bool                        a_homogeneousBC){
  CH_TIME("nwomfconductivityop::AMROperatorNF");
#if verb
  pout() << "nwomfconductivityop::amroperatornf"<< endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);
  
  for (int i=0; i<m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_LofPhi);
    mfalias::aliasMF(*m_alias[1], i, a_phi);
    mfalias::aliasMF(*m_alias[2], i, a_phiCoarse);

    m_ebops[i]->AMROperatorNF(*m_alias[0], *m_alias[1], *m_alias[2], a_homogeneousBC);
    //    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], m_alias[2], a_homogeneousBC, false);
  }
}



void nwomfconductivityop::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
					 const LevelData<MFCellFAB>& a_correction,
					 const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("nwomfconductivityop::AMRUpdateResidual");
#if verb
  pout() << "nwomfconductivityop::amrupdateresidual"<< endl;
#endif

  this->update_bc(a_correction, true);

  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_residual);
    mfalias::aliasMF(*m_alias[1], i, a_correction);
    mfalias::aliasMF(*m_alias[2], i, a_coarseCorrection);
    
    m_ebops[i]->AMRUpdateResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void nwomfconductivityop::AMRRestrict(LevelData<MFCellFAB>&       a_resCoarse,
				   const LevelData<MFCellFAB>& a_residual,
				   const LevelData<MFCellFAB>& a_correction,
				   const LevelData<MFCellFAB>& a_coarseCorrection,
				   bool                        a_skip_res){
  CH_TIME("nwomfconductivityop::AMRRestrict");
#if verb
  pout() << "nwomfconductivityop::amrrestrict"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_resCoarse);
    mfalias::aliasMF(*m_alias[1], i, a_residual);
    mfalias::aliasMF(*m_alias[2], i, a_correction);
    mfalias::aliasMF(*m_alias[3], i, a_coarseCorrection);
    
    m_ebops[i]->AMRRestrict(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_skip_res);
  }
}

void nwomfconductivityop::AMRProlong(LevelData<MFCellFAB>&       a_correction,
				  const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("nwomfconductivityop::AMRProlong");
#if verb
  pout() << "nwomfconductivityop::amrprolong"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_correction);
    mfalias::aliasMF(*m_alias[1], i, a_coarseCorrection);
    
    m_ebops[i]->AMRProlong(*m_alias[0], *m_alias[1]);
  }
}



Real nwomfconductivityop::AMRNorm(const LevelData<MFCellFAB>& a_coar_resid,
			       const LevelData<MFCellFAB>& a_fine_resid,
			       const int&                  a_ref_rat,
			       const int&                  a_ord){
  CH_TIME("mf_helmholtzop::AMRNorm");
#if verb
  pout() << "nwomfconductivityop::amrnorm"<< endl;
#endif
  Real m = 0;

  for (int i = 0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_coar_resid);
    mfalias::aliasMF(*m_alias[1], i, a_fine_resid);

    Real norm = m_ebops[i]->AMRNorm(*m_alias[0], *m_alias[1], a_ref_rat, a_ord);

    m = Max(m, norm);
  }

  return m;
}

int nwomfconductivityop::refToCoarser(){
#if verb
  pout() << "nwomfconductivityop::reftocoarser"<< endl;
#endif
  return m_ref_to_coarser;
}
