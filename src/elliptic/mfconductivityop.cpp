/*!
  @file mfconductivityop.cpp
  @brief Implementation of mfconductivityop.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "mfdirichletconductivityebbc.H"
#include "mfconductivityop.H"
#include "mfalias.H"

#include <DirichletConductivityDomainBC.H>
#include <MFLevelDataOps.H>
#include <BaseIVFactory.H>
#include <EBAMRDataOps.H>

#define verb 0

mfconductivityop::mfconductivityop(){

}
  

mfconductivityop::~mfconductivityop(){
  
}

void mfconductivityop::define(const RefCountedPtr<mfis>&                    a_mfis,
			      const RefCountedPtr<BaseDomainBC>&            a_dombc,
			      const RefCountedPtr<LevelData<MFCellFAB> >&   a_aco,
			      const RefCountedPtr<LevelData<MFFluxFAB> >&   a_bco,
			      const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_bco_irreg,
			      const MFQuadCFInterp&                         a_quadcfi,
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



#if verb
  MayDay::Warning("mfconductivityop::mfconductivityop - remember to check how many aliasing holders we need");
#endif
  const int num_alias  = 4;

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


    RefCountedPtr<EBQuadCFInterp> quadcfi;
    if(a_has_coar){
      quadcfi = a_quadcfi.get_quadcfi_ptr(iphase);
      CH_assert(!quadcfi.isNull());
    }


    m_acoeffs[iphase]     = RefCountedPtr<LevelData<EBCellFAB> >        (new LevelData<EBCellFAB>());
    m_bcoeffs[iphase]     = RefCountedPtr<LevelData<EBFluxFAB> >        (new LevelData<EBFluxFAB>());
    m_bcoeffs_irr[iphase] = RefCountedPtr<LevelData<BaseIVFAB<Real> > > (new LevelData<BaseIVFAB<Real> >());

    
#if 1 // In the long run, we will replace these with more general boundary conditions.
#if verb
    MayDay::Warning("mfconductivityop::mfconductivityop - fix up boundary conditions!");
#endif
    DirichletConductivityDomainBCFactory bcfact;
    bcfact.setValue(0.0);
    
    RefCountedPtr<DirichletConductivityDomainBC> dombc = RefCountedPtr<DirichletConductivityDomainBC>
      (bcfact.create(a_domain, ebisl, a_dx*RealVect::Unit));



#if 1 // Testing
    m_ebbc[iphase] = RefCountedPtr<mfdirichletconductivityebbc> (new mfdirichletconductivityebbc(a_domain,
												 ebisl,
												 a_dx*RealVect::Unit,
												 &a_ghost_phi,
												 &a_ghost_rhs,
												 iphase)),
#else
      
      m_ebbc[iphase]        = RefCountedPtr<DirichletConductivityEBBC>    (new DirichletConductivityEBBC(a_domain,
													 ebisl,
													 a_dx*RealVect::Unit,
													 &a_ghost_phi,
													 &a_ghost_rhs));
#endif

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


#endif

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
    pout() << "mfconductivityop::creating factory" << endl;
#endif
    m_ebops[iphase] = RefCountedPtr<ebconductivityop> (new ebconductivityop(eblg_fine,
									    eblg,
									    eblg_coar,
									    eblg_mg,
									    quadcfi,
									    dombc,
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
#if verb
    pout() << "mfconductivityop::done creating factory" << endl;
#endif
  }



  MFCellFactory* factory = new MFCellFactory(layouts, comps);
  RefCountedPtr<DataFactory<MFCellFAB> > fac(factory);
  m_tmp.define(a_mflg.get_grids(), m_ncomp, m_ghost_rhs, *factory);
  m_weights.define(a_mflg.get_grids(), m_ncomp, m_ghost_phi, *factory);
	       
  m_ops.define(fac);
}

void mfconductivityop::set_jump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump){
  m_jump = a_jump;
}

void mfconductivityop::set_electrodes(const Vector<electrode>& a_electrodes){
#if verb
  pout() << "mfconductivityop::set_electrodes"<< endl;
#endif

  m_electrodes = a_electrodes;
  this->set_bc_from_levelset();
}

void mfconductivityop::update_bc(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
#if 0
  pout() << "mfconductivityop::update_bc"<< endl;
#endif

  this->set_bc_from_levelset();
  this->set_bc_from_matching(a_phi, a_homogeneous);
}

void mfconductivityop::set_bc_from_levelset(){
  CH_TIME("mfconductivityop::set_bc_from_levelset");
  
  for (int iphase = 0; iphase < m_phases; iphase++){
    LevelData<BaseIVFAB<Real> >& val = *m_dirival[iphase];

    for (DataIterator dit = val.dataIterator(); dit.ok(); ++dit){
      const IntVectSet& ivs  = val[dit()].getIVS();
      const EBGraph& ebgraph = val[dit()].getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect& pos = EBArith::getVofLocation(vof, m_dx, m_origin);

	int  func = 0;
	Real dist = m_electrodes[0].get_function()->value(pos);

	for (int i = 1; i < m_electrodes.size(); i++){
	  Real cur_val = (m_electrodes[i].get_function())->value(pos);
	  if(Abs(cur_val) < dist){
	    func = i;
	    dist  = cur_val;
	  }
	}

	val[dit()](vof, 0) = m_electrodes[func].is_live() ? 1 : 0;
      }
    }
  }
}

void mfconductivityop::set_bc_from_matching(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
  for (int iphase = 0; iphase < m_phases; iphase++){
    m_jumpbc->match_bc(*m_dirival[iphase], a_phi, a_homogeneous);
  }
}

void mfconductivityop::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta){
  CH_TIME("mfconductivityop::setAlphaAndBeta");
#if verb
  pout() << "mfconductivityop::setalaphandbeta"<< endl;
#endif
  for (int iphase = 0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setAlphaAndBeta(a_alpha, a_beta);
#if 0
    MayDay::Warning("mfconductivityop::setAlphaAndBeta - MFPoisson multiplies by acoef, why?");
#endif
  }
}

void mfconductivityop::diagonalScale(LevelData<MFCellFAB>& a_rhs){
  CH_TIME("mfconductivityop::diagonalScale");
#if verb
  pout() << "mfconductivityop::setdiagonalscale"<< endl;
#endif
  MayDay::Abort("where did i get called?");

  // // Operator diagonal scale
  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->diagonalScale(*m_alias[0], true);
  }

  MayDay::Abort("mfconductivityop::diagonalScale - add scaling for area fractions");
  //  MFLevelDataOps::kappaWeight(a_rhs); // This came from MFPoissonOp
}

void mfconductivityop::divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "mfconductivityop::dividebyidenticyoef"<< endl;
#endif

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->divideByIdentityCoef(*m_alias[0]);
  }
}

void mfconductivityop::applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
					 const LevelData<MFCellFAB>& a_phi){
#if verb
  pout() << "mfconductivityop::applyopnoboundary"<< endl;
#endif
  CH_TIME("mfconductivityop::applyOpNoBoundary");

  this->update_bc(a_phi, true);
  for (int iphase=0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_opPhi);
    mfalias::aliasMF(*m_alias[1], iphase, a_phi);
    m_ebops[iphase]->applyOpNoBoundary(*m_alias[0], *m_alias[1]);
  }
}

void mfconductivityop::setTime(Real a_oldTime, Real a_mu, Real a_dt){
  CH_TIME("mfconductivityop::setTime");
#if verb
  pout() << "mfconductivityop::settime"<< endl;
#endif
  for (int iphase=0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setTime(a_oldTime, a_mu, a_dt);
  }

#if verb
  MayDay::Warning("mfconductivityop::setTime - jump bc should also be time dependent");
#endif
}

void mfconductivityop::residual(LevelData<MFCellFAB>&        a_lhs,
				const LevelData<MFCellFAB>&  a_phi,
				const LevelData<MFCellFAB>&  a_rhs,
				bool                         a_homogeneous){
  CH_TIME("mfconductivityop::residual");
#if verb
  pout() << "mfconductivityop::residual"<< endl;
#endif
  
  this->applyOp(a_lhs, a_phi, a_homogeneous);

  // Make residual rhs - L(phi)
  this->scale(a_lhs,      -1.0);
  this->incr(a_lhs, a_rhs, 1.0);

}

void mfconductivityop::preCond(LevelData<MFCellFAB>&       a_correction,
			       const LevelData<MFCellFAB>& a_residual){
  CH_TIME("mfconductivityop::preCond");
#if verb
  pout() << "mfconductivityop::precond"<< endl;
#endif
#if 0

  this->relax(a_correction, a_residual, 40);
#endif
}

void mfconductivityop::applyOp(LevelData<MFCellFAB>&        a_lhs,
			       const LevelData<MFCellFAB>&  a_phi,
			       bool                         a_homogeneous){
  CH_TIME("mfconductivityop::applyOp");
#if verb
  pout() << "mfconductivityop::applyop"<< endl;
#endif
#if verb
  MayDay::Warning("mfconductivityop::applyOp - the matching condition should be updated first");
#endif

  this->update_bc(a_phi, a_homogeneous);

  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_lhs);
    mfalias::aliasMF(*m_alias[1], i, a_phi);
    
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneous, true);
  }
}

void mfconductivityop::create(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("mfconductivityop::create");
#if verb
  pout() << "mfconductivityop::create"<< endl;
#endif
  m_ops.create(a_lhs, a_rhs);
}

void mfconductivityop::createCoarsened(LevelData<MFCellFAB>&       a_lhs,
				       const LevelData<MFCellFAB>& a_rhs,
				       const int&                  a_refRat) {
  CH_TIME("mfconductivityop::createCoarsened");
#if verb
  pout() << "mfconductivityop::createCoarsened"<< endl;
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
  pout() << "mfconductivityop::createCoarsened - done"<< endl;
#endif
}

void mfconductivityop::assign(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "mfconductivityop::assign"<< endl;
#endif
  m_ops.assign(a_lhs, a_rhs);
#if verb
  pout() << "mfconductivityop::assign - doen"<< endl;
#endif
}


Real mfconductivityop::dotProduct(const LevelData<MFCellFAB>& a_data1,
				  const LevelData<MFCellFAB>& a_data2){
#if verb
  pout() << "mfconductivityop::dotproduct"<< endl;
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

void mfconductivityop::incr(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_rhs,
			    Real                        a_scale){
#if verb
  pout() << "mfconductivityop::incr"<< endl;
#endif
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
#if verb
  pout() << "mfconductivityop::incr - done"<< endl;
#endif
}

void mfconductivityop::axby(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_x,
			    const LevelData<MFCellFAB>& a_y,
			    Real a,
			    Real b){
#if verb
  pout() << "mfconductivityop::axby"<< endl;
#endif

  m_ops.axby(a_lhs, a_x, a_y, a, b);
}

void mfconductivityop::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale){
  CH_TIME("mfconductivityop::scale");
#if verb
  pout() << "mfconductivityop::scale"<< endl;
#endif
  m_ops.scale(a_lhs, a_scale);
#if verb
  pout() << "mfconductivityop::scale - doen"<< endl;
#endif
}

Real mfconductivityop::norm(const LevelData<MFCellFAB>& a_x, int a_ord){
  CH_TIME("mfconductivityop::norm");
#if verb
  pout() << "mfconductivityop::norm"<< endl;
#endif
  Real volume;
  Real rtn = this->kappaNorm(volume, a_x, a_ord);
#if verb
  pout() << "mfconductivityop::norm - done"<< endl;
#endif
  return rtn;
}

Real mfconductivityop::kappaNorm(Real&                       a_volume,
				 const LevelData<MFCellFAB>& a_data,
				 int                         a_p) const {

#if verb
  pout() << "mfconductivityop::kappaNorm"<< endl;
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

void mfconductivityop::setToZero(LevelData<MFCellFAB>& a_x){
  CH_TIME("mfconductivityop::setToZero");
#if verb
  pout() << "mfconductivityop::setToZero"<< endl;
#endif
  m_ops.setToZero(a_x);
#if verb
  pout() << "mfconductivityop::setToZero - done"<< endl;
#endif
}

void mfconductivityop::relax(LevelData<MFCellFAB>&       a_e,
			     const LevelData<MFCellFAB>& a_residual,
			     int                         iterations){
  CH_TIME("mfconductivityop::relax");
#if verb
  pout() << "mfconductivityop::relax"<< endl;
#endif

#if 0
  for (int i=0; i < iterations; i++){
    if (m_relax == 0){
      this->levelJacobi(a_e, a_residual);
    }
    else if (m_relax == 1){
      this->levelMulticolorGS(a_e, a_residual);
    }
    else if (m_relax == 2){
      this->levelGSRB(a_e, a_residual);
    }
    else{
      MayDay::Error("mfconductivityop: Invalid relaxation type");
    }
  }
#else
  const bool homogeneous = true;
  for (int i = 0; i < iterations; i++){
    this->update_bc(a_e, homogeneous);
    this->levelJacobi(a_e, a_residual, iterations);
  }
#endif
}

void mfconductivityop::levelJacobi(LevelData<MFCellFAB>&       a_phi,
				   const LevelData<MFCellFAB>& a_rhs,
				   const int                   a_iterations){
  CH_TIME("mfconductivityop::levelJacobi");
  LevelData<MFCellFAB>& resid = m_tmp;
  bool homogeneous = true;

  //  this->residual(resid, a_phi, a_rhs, true);

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_phi);
    mfalias::aliasMF(*m_alias[1], iphase, a_rhs);
    //    mfalias::aliasMF(*m_alias[2], iphase, resid);
    
    //    m_ebops[iphase]->relax(*m_alias[0], *m_alias[1], 1);
    m_ebops[iphase]->relaxGauSai(*m_alias[0], *m_alias[1], 1);
    //m_ebops[iphase]->relaxGSRBFast(*m_alias[0], *m_alias[1], 1);
  }
}

void mfconductivityop::createCoarser(LevelData<MFCellFAB>&       a_coarse,
				     const LevelData<MFCellFAB>& a_fine,
				     bool                        ghosted){
  CH_TIME("mfconductivityop::createCoarser");
#if verb
  pout() << "mfconductivityop::createCoarser"<< endl;
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

void mfconductivityop::restrictResidual(LevelData<MFCellFAB>&       a_resCoarse,
					LevelData<MFCellFAB>&       a_phiFine,
					const LevelData<MFCellFAB>& a_rhsFine){
  CH_TIME("mfconductivityop::restrictResidual");
#if verb
  pout() << "mfconductivityop::restrictResidual"<< endl;
#endif
  for (int i=0; i < m_phases; i++) {
    mfalias::aliasMF(*m_alias[0], i, a_resCoarse);
    mfalias::aliasMF(*m_alias[1], i, a_phiFine);
    mfalias::aliasMF(*m_alias[2], i, a_rhsFine);
    
    m_ebops[i]->restrictResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void mfconductivityop::prolongIncrement(LevelData<MFCellFAB>&       a_phiThisLevel,
					const LevelData<MFCellFAB>& a_correctCoarse){
  CH_TIME("mfconductivityop::prolongIncrement");
#if verb
  pout() << "mfconductivityop::prolongIncrement"<< endl;
#endif
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_phiThisLevel);
    mfalias::aliasMF(*m_alias[1], i, a_correctCoarse);
    
    m_ebops[i]->prolongIncrement(*m_alias[0], *m_alias[1]);
  }
}

void mfconductivityop::AMRResidual(LevelData<MFCellFAB>&       a_residual,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   const LevelData<MFCellFAB>& a_rhs,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("mfconductivityop::AMRResidual");
#if verb
  pout() << "mfconductivityop::amrresidual"<< endl;
#endif

  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void mfconductivityop::AMROperator(LevelData<MFCellFAB>&       a_LofPhi,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("mfconductivityop::AMROperator");
#if verb
  pout() << "mfconductivityop::amroperator" << endl;
#endif

#if verb
  MayDay::Warning("mfconductivityop::AMROperator - must update BCfirst");
#endif
  this->update_bc(a_phi, a_homogeneousBC);

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_LofPhi);
    mfalias::aliasMF(*m_alias[1], iphase, a_phiFine);
    mfalias::aliasMF(*m_alias[2], iphase, a_phi);
    mfalias::aliasMF(*m_alias[3], iphase, a_phiCoarse);

    mfconductivityop* finerOp = (mfconductivityop*) a_finerOp;

    //    pout() << "amroperator, phase = " << iphase << endl;
    m_ebops[iphase]->AMROperator(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_homogeneousBC, finerOp->m_ebops[iphase]);
    //    pout() << "amroperator, phase = " << iphase << " done" << endl;
  }

										 
  //    m_ebops[iphase]->applyOp(*m_alias[0], *m_alias[1], m_alias[2], a_homogeneousBC, false);
  //    m_ebops[iphase]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[iphase]);

#if verb
  pout() << "mfconductivityop::amroperator - done" << endl;
#endif
}

void mfconductivityop::AMRResidualNC(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("mfconductivityop::AMRResidualNC");
#if verb
  pout() << "mfconductivityop::amrresidualnc"<< endl;
#endif

  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);

#if verb
  pout() << "mfconductivityop::amrresidualnc - done"<< endl;
#endif
}

void mfconductivityop::AMROperatorNC(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("mfconductivityop::AMROperatorNC");
#if verb
  pout() << "mfconductivityop::amroperatornc"<< endl;
#endif
#if verb
  MayDay::Warning("mfconductivityop::AMROperatorNC - must update BCfirst");
#endif

  this->update_bc(a_phi, a_homogeneousBC);
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_LofPhi);
    mfalias::aliasMF(*m_alias[1], i, a_phiFine);
    mfalias::aliasMF(*m_alias[2], i, a_phi);

    mfconductivityop* finerOp = (mfconductivityop*) a_finerOp;
    m_ebops[i]->AMROperatorNC(*m_alias[0], *m_alias[1], *m_alias[2], a_homogeneousBC, finerOp->m_ebops[i]);

    // m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneousBC, true);

    // mfconductivityop* finerOp = (mfconductivityop*) a_finerOp;
    // m_ebops[i]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[i]);
  }
#if verb
  pout() << "mfconductivityop::amroperatornc - done"<< endl;
#endif
}


void mfconductivityop::AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC){
  CH_TIME("mfconductivityop::AMRResidualNF");
#if verb
  pout() << "mfconductivityop::amrresidualnf"<< endl;
#endif

  this->AMROperatorNF(a_residual, a_phi, a_phiCoarse, a_homogeneousBC);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void mfconductivityop::AMROperatorNF(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     bool                        a_homogeneousBC){
  CH_TIME("mfconductivityop::AMROperatorNF");
#if verb
  pout() << "mfconductivityop::amroperatornf"<< endl;
#endif

#if verb
  MayDay::Warning("mfconductivityop::AMROperatorNC - must update BCfirst");
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



void mfconductivityop::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
					 const LevelData<MFCellFAB>& a_correction,
					 const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("mfconductivityop::AMRUpdateResidual");
#if verb
  pout() << "mfconductivityop::amrupdateresidual"<< endl;
#endif

#if verb
  MayDay::Warning("mfconductivityop::applyOp - the matching condition should be updated first");
#endif
  this->update_bc(a_correction, true);

  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_residual);
    mfalias::aliasMF(*m_alias[1], i, a_correction);
    mfalias::aliasMF(*m_alias[2], i, a_coarseCorrection);
    
    m_ebops[i]->AMRUpdateResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void mfconductivityop::AMRRestrict(LevelData<MFCellFAB>&       a_resCoarse,
				   const LevelData<MFCellFAB>& a_residual,
				   const LevelData<MFCellFAB>& a_correction,
				   const LevelData<MFCellFAB>& a_coarseCorrection,
				   bool                        a_skip_res){
  CH_TIME("mfconductivityop::AMRRestrict");
#if verb
  pout() << "mfconductivityop::amrrestrict"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_resCoarse);
    mfalias::aliasMF(*m_alias[1], i, a_residual);
    mfalias::aliasMF(*m_alias[2], i, a_correction);
    mfalias::aliasMF(*m_alias[3], i, a_coarseCorrection);
    
    m_ebops[i]->AMRRestrict(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_skip_res);
  }
}

void mfconductivityop::AMRProlong(LevelData<MFCellFAB>&       a_correction,
				  const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("mfconductivityop::AMRProlong");
#if verb
  pout() << "mfconductivityop::amrprolong"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_correction);
    mfalias::aliasMF(*m_alias[1], i, a_coarseCorrection);
    
    m_ebops[i]->AMRProlong(*m_alias[0], *m_alias[1]);
  }
}



Real mfconductivityop::AMRNorm(const LevelData<MFCellFAB>& a_coar_resid,
			       const LevelData<MFCellFAB>& a_fine_resid,
			       const int&                  a_ref_rat,
			       const int&                  a_ord){
  CH_TIME("mf_helmholtzop::AMRNorm");
#if verb
  pout() << "mfconductivityop::amrnorm"<< endl;
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

int mfconductivityop::refToCoarser(){
#if verb
  pout() << "mfconductivityop::reftocoarser"<< endl;
#endif
  return m_ref_to_coarser;
}
