/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MfHelmholtzOp.cpp
  @brief  Implementation of CD_MfHelmholtzOp.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <DirichletConductivityDomainBC.H>
#include <MFLevelDataOps.H>
#include <BaseIVFactory.H>
#include <EBAMRDataOps.H>

// Our includes
#include <CD_MfElectrostaticDirichletEbBc.H>
#include <CD_MfHelmholtzOp.H>
#include <CD_MultifluidAlias.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>
  
#define verb 0

MfHelmholtzOp::MfHelmholtzOp(){
#if verb
  pout() << "MfHelmholtzOp::MfHelmholtzOp" << endl;
#endif
}
  

MfHelmholtzOp::~MfHelmholtzOp(){
  for (int i = 0; i < m_alias.size(); i++){
    delete m_alias[i];
    m_alias[i] = NULL;
  }

  //  delete m_jumpbc;

  for (int iphase= 0; iphase < m_phases; iphase++){
    //    delete m_ebbc[iphase];
  }
}

void MfHelmholtzOp::define(const RefCountedPtr<MultiFluidIndexSpace>&                    a_multiFluidIndexSpace,
			      const RefCountedPtr<BaseDomainBCFactory>&     a_dombc,
			      const RefCountedPtr<LevelData<MFCellFAB> >&   a_aco,
			      const RefCountedPtr<LevelData<MFFluxFAB> >&   a_bco,
			      const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_bco_irreg,
			      const MFQuadCFInterp&                         a_quadcfi,
			      const MFFluxReg&                          a_fluxreg,
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
			      const Real&                                   a_lengthScale,
			      const Real&                                   a_dx,
			      const Real&                                   a_dx_coar,
			      const Real&                                   a_alpha,
			      const Real&                                   a_beta,
			      const RealVect&                               a_origin){


  const int numPhases = a_multiFluidIndexSpace->numPhases();
  
  const int num_alias  = 6;

  m_multifluidIndexSpace = a_multiFluidIndexSpace;
  m_ncomp = 1;
  m_relax = a_relax_type;
  m_domain = a_domain;
  m_phases = numPhases;
  m_ebops.resize(numPhases);
  m_ebbc.resize(numPhases);
  m_aCoefeffs.resize(numPhases);
  m_bcoeffs.resize(numPhases);
  m_bcoeffs_irr.resize(numPhases);
  m_alias.resize(num_alias);
  m_ref_to_coarser = a_ref_to_coar;
  m_ghost_phi = a_ghost_phi;
  m_ghost_rhs = a_ghost_rhs;
  m_origin = a_origin;
  m_dx = a_dx;
  m_dirival.resize(numPhases);
  m_multifluid = m_multifluidIndexSpace->numPhases() > 1;
  m_lengthScale = a_lengthScale;

  if(a_has_mg){
    m_mflg_coar_mg = a_mflg_coar_mg;
  }

  for (int i = 0; i < num_alias; i++){
    m_alias[i] = new LevelData<EBCellFAB>();
  }

  // Object for matching boundary conditions. Native EBBC is data-based Dirichlet
  m_jumpbc = RefCountedPtr<JumpBc> (new JumpBc(a_mflg, *a_bco_irreg, a_dx, a_order_ebbc, (a_mflg.getEBLevelGrid(0)).getCFIVS()));

  Vector<EBISLayout> layouts(numPhases);
  Vector<int> comps(numPhases);;


  
  for (int iphase = 0; iphase < numPhases; iphase++){

    const EBLevelGrid& eblg      = a_mflg.getEBLevelGrid(iphase);
    const EBISLayout&  ebisl     = eblg.getEBISL();

    layouts[iphase] = ebisl;
    comps[iphase]   = m_ncomp;


    EBLevelGrid eblg_fine;
    EBLevelGrid eblg_coar;
    EBLevelGrid eblg_mg;

    if(a_has_fine){
      eblg_fine = a_mflg_fine.getEBLevelGrid(iphase);
    }
    if(a_has_coar){
      eblg_coar = a_mflg_coar.getEBLevelGrid(iphase);
    }
    if(a_has_mg){
      eblg_mg   = a_mflg_coar_mg.getEBLevelGrid(iphase);
    }


    RefCountedPtr<EBFluxRegister> fastFR;
    RefCountedPtr<EBQuadCFInterp> quadcfi;
    if(a_has_coar){
      quadcfi = a_quadcfi.getEBQuadCFInterpPointer(iphase);
      CH_assert(!quadcfi.isNull());
    }
    if(a_has_fine){
      fastFR = a_fluxreg.getFluxRegPointer(iphase);
    }


    // Coefficients
    m_aCoefeffs[iphase]     = RefCountedPtr<LevelData<EBCellFAB> >        (new LevelData<EBCellFAB>());
    m_bcoeffs[iphase]     = RefCountedPtr<LevelData<EBFluxFAB> >        (new LevelData<EBFluxFAB>());
    m_bcoeffs_irr[iphase] = RefCountedPtr<LevelData<BaseIVFAB<Real> > > (new LevelData<BaseIVFAB<Real> >());

    
    // Domain BC
    ConductivityBaseDomainBC* bc = (ConductivityBaseDomainBC*) a_dombc->create(a_domain, ebisl, a_dx*RealVect::Unit);
    RefCountedPtr<ConductivityBaseDomainBC> dbc(bc);




    // Create storage for data-based dirichlet boundary conditions
    LayoutData<IntVectSet> ivs(eblg.getDBL());
    for (DataIterator dit = ivs.dataIterator(); dit.ok(); ++dit){
      ivs[dit()] = ebisl[dit()].getIrregIVS(eblg.getDBL().get(dit()));
    }

    BaseIVFactory<Real> ivfact(ebisl, ivs);
    m_dirival[iphase] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> > (eblg.getDBL(), 1, IntVect::Zero, ivfact));



    // EB BC
    m_ebbc[iphase] = RefCountedPtr<MfElectrostaticDirichletEbBc> (new MfElectrostaticDirichletEbBc(a_domain,
												 ebisl,
												 a_dx*RealVect::Unit,
												 &a_ghost_phi,
												 &a_ghost_rhs,
												 iphase));
    m_ebbc[iphase]->setJumpObject(m_jumpbc);
    m_ebbc[iphase]->setOrder(a_order_ebbc);
    m_ebbc[iphase]->defineIVS(a_mflg);

    EBLevelDataOps::setVal(*m_dirival[iphase], 0.0);
    m_ebbc[iphase]->setData(m_dirival[iphase]);

    MultifluidAlias::aliasMF(*m_aCoefeffs[iphase],     iphase, *a_aco);
    MultifluidAlias::aliasMF(*m_bcoeffs[iphase],     iphase, *a_bco);
    MultifluidAlias::aliasMF(*m_bcoeffs_irr[iphase], iphase, *a_bco_irreg);

    const RefCountedPtr<LevelData<EBCellFAB> >&        aco     = m_aCoefeffs[iphase];
    const RefCountedPtr<LevelData<EBFluxFAB> >&        bco     = m_bcoeffs[iphase];
    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& bco_irr = m_bcoeffs_irr[iphase];
    const RefCountedPtr<MfElectrostaticDirichletEbBc>&  ebbc    = m_ebbc[iphase];

    const Real alpha = a_alpha;
    const Real beta  = a_beta;

#if verb
    pout() << "MfHelmholtzOp::creating oper" << endl;
#endif
    m_ebops[iphase] = RefCountedPtr<EbHelmholtzOp> (new EbHelmholtzOp(eblg_fine,
									    eblg,
									    eblg_coar,
									    eblg_mg,
									    quadcfi,
									    fastFR,
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
#if verb
    pout() << "MfHelmholtzOp::done creating oper" << endl;
#endif
  }

  MFCellFactory* factory = new MFCellFactory(layouts, comps);
  RefCountedPtr<DataFactory<MFCellFAB> > fac(factory);
  m_ops.define(fac);

  EBArith::getMultiColors(m_colors);

#if 0 // Define aggregate stencils in JumpBc object. Only for the new JumpBc class
  LevelData<MFCellFAB> dummy(a_aco->disjointBoxLayout(), 1, a_ghost_phi*IntVect::Unit, *factory);
  for (int iphase = 0; iphase < numPhases; iphase++){
    m_jumpbc->define_agg_stencils(*m_dirival[iphase], dummy, iphase);
  }
#endif
}

void MfHelmholtzOp::setJump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump){
  m_jump = a_jump;
}

void MfHelmholtzOp::setTime(Real* a_time){
  m_time = a_time;
}

void MfHelmholtzOp::setDirichletEbBc(const ElectrostaticEbBc& a_ebbc){
  m_electrostaticEbBc = a_ebbc;

  this->set_bc_from_levelset();
}

void MfHelmholtzOp::update_bc(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
  CH_TIME("MfHelmholtzOp::update_bc");
#if verb
  pout() << "MfHelmholtzOp::update_bc"<< endl;
#endif

  //  const Real t0 = MPI_Wtime();
  LevelData<MFCellFAB>* phi = const_cast<LevelData<MFCellFAB>* > (&a_phi);
  phi->exchange();
  //  const Real t1 = MPI_Wtime();

  //  this->set_bc_from_levelset();
  this->set_bc_from_matching(a_phi, a_homogeneous);
  //  const Real t2 = MPI_Wtime();

  // pout() << "update bc exchange = " << 100.*(t1-t0)/(t2-t0) << endl;
  // pout() << "update bc matching = " << 100.*(t2-t1)/(t2-t0) << endl;

#if verb
  pout() << "MfHelmholtzOp::update_bc - done" << endl;
#endif
}

void MfHelmholtzOp::update_bc(const LevelData<MFCellFAB>& a_phi, DataIterator& a_dit, const bool a_homogeneous){
  LevelData<MFCellFAB>* phi = const_cast<LevelData<MFCellFAB>* > (&a_phi);
  phi->exchange();
  this->set_bc_from_matching(a_phi, a_dit, a_homogeneous);
}

void MfHelmholtzOp::set_bc_from_levelset(){
  CH_TIME("MfHelmholtzOp::set_bc_from_levelset");
#if verb
  pout() << "MfHelmholtzOp::set_bc_from_levelset"<< endl;
#endif


  const int comp = 0;

  const Real physDx = m_dx/m_lengthScale;
  
  for (int iphase = 0; iphase < m_phases; iphase++){
    LevelData<BaseIVFAB<Real> >& val = *m_dirival[iphase];

    for (DataIterator dit = val.dataIterator(); dit.ok(); ++dit){
      const IntVectSet& ivs  = val[dit()].getIVS();
      const EBGraph& ebgraph = val[dit()].getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, physDx, m_origin);

	const std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction> >& electrodeBcs = m_electrostaticEbBc.getBcs();
	
	if(electrodeBcs.size() > 0){
	  int  func = -1;
	  Real dist = std::numeric_limits<Real>::infinity();

	  for (int i = 0; i < electrodeBcs.size(); i++){
	    Real cur_val = (electrodeBcs[i].first.getImplicitFunction())->value(pos);
	    if(std::abs(cur_val) < std::abs(dist)){
	      func  = i;
	      dist  = cur_val;
	    }
	  }


	  ElectrostaticEbBc::BcFunction bcFunction = electrodeBcs[func].second;
	  
	  val[dit()](vof, comp) = bcFunction(pos, 0.0);
	}
      }
    }
  }
}

void MfHelmholtzOp::set_bc_from_matching(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
  DataIterator dit = a_phi.dataIterator();
  this->set_bc_from_matching(a_phi, dit, a_homogeneous);
}

void MfHelmholtzOp::set_bc_from_matching(const LevelData<MFCellFAB>& a_phi, DataIterator& a_dit, const bool a_homogeneous){
  CH_TIME("MfHelmholtzOp::set_bc_from_matching");
#if verb
  pout() << "MfHelmholtzOp::set_bc_from_matching"<< endl;
#endif

  if(m_multifluid){
    for (int iphase = 0; iphase <= 1; iphase++){
      m_jumpbc->matchBc(*m_dirival[iphase], *m_jump, a_phi, a_dit, a_homogeneous);
    }
  }

#if verb
  pout() << "MfHelmholtzOp::set_bc_from_matching - done"<< endl;
#endif
}

void MfHelmholtzOp::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta){
  CH_TIME("MfHelmholtzOp::setAlphaAndBeta");
#if verb
  pout() << "MfHelmholtzOp::setalaphandbeta"<< endl;
#endif
  for (int iphase = 0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setAlphaAndBeta(a_alpha, a_beta);
  }
}

void MfHelmholtzOp::diagonalScale(LevelData<MFCellFAB>& a_rhs){
  CH_TIME("MfHelmholtzOp::diagonalScale");
#if verb
  pout() << "MfHelmholtzOp::setdiagonalscale"<< endl;
#endif
  MayDay::Abort("MfHelmholtzOp::DiagonalScale - where did i get called?");

  // // Operator diagonal scale
  for (int iphase = 0; iphase < m_phases; iphase++){
    MultifluidAlias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->diagonalScale(*m_alias[0], true);
  }

  MayDay::Abort("MfHelmholtzOp::diagonalScale - add scaling for area fractions");
  //  MFLevelDataOps::kappaWeight(a_rhs); // This came from MFPoissonOp
}

void MfHelmholtzOp::divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "MfHelmholtzOp::dividebyidenticyoef"<< endl;
#endif

  for (int iphase = 0; iphase < m_phases; iphase++){
    MultifluidAlias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->divideByIdentityCoef(*m_alias[0]);
  }
}

void MfHelmholtzOp::applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
					 const LevelData<MFCellFAB>& a_phi){
#if verb
  pout() << "MfHelmholtzOp::applyopnoboundary"<< endl;
#endif
  CH_TIME("MfHelmholtzOp::applyOpNoBoundary");

  this->update_bc(a_phi, true);
  for (int iphase=0; iphase < m_phases; iphase++){
    MultifluidAlias::aliasMF(*m_alias[0], iphase, a_opPhi);
    MultifluidAlias::aliasMF(*m_alias[1], iphase, a_phi);
    m_ebops[iphase]->applyOpNoBoundary(*m_alias[0], *m_alias[1]);
  }
}

void MfHelmholtzOp::setTime(Real a_oldTime, Real a_mu, Real a_dt){
  CH_TIME("MfHelmholtzOp::setTime");
#if verb
  pout() << "MfHelmholtzOp::settime"<< endl;
#endif
  for (int iphase=0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setTime(a_oldTime, a_mu, a_dt);
  }
}

void MfHelmholtzOp::residual(LevelData<MFCellFAB>&        a_lhs,
				const LevelData<MFCellFAB>&  a_phi,
				const LevelData<MFCellFAB>&  a_rhs,
				bool                         a_homogeneous){
  CH_TIME("MfHelmholtzOp::residual");
#if verb
  pout() << "MfHelmholtzOp::residual"<< endl;
#endif
  
  this->applyOp(a_lhs, a_phi, a_homogeneous);

  // Make residual rhs - L(phi)
  this->scale(a_lhs,      -1.0);
  this->incr(a_lhs, a_rhs, 1.0);

}

void MfHelmholtzOp::preCond(LevelData<MFCellFAB>&       a_correction,
			       const LevelData<MFCellFAB>& a_residual){
  CH_TIME("MfHelmholtzOp::preCond");
#if verb
  pout() << "MfHelmholtzOp::precond"<< endl;
#endif

  this->relax(a_correction, a_residual, 40);
}

void MfHelmholtzOp::applyOp(LevelData<MFCellFAB>&        a_lhs,
			       const LevelData<MFCellFAB>&  a_phi,
			       bool                         a_homogeneous){
  CH_TIME("MfHelmholtzOp::applyOp");
#if verb
  pout() << "MfHelmholtzOp::applyop"<< endl;
#endif

  DataIterator dit = a_lhs.dataIterator();
  this->applyOp(a_lhs, a_phi, dit, a_homogeneous);

}

void MfHelmholtzOp::applyOp(LevelData<MFCellFAB>&        a_lhs,
			       const LevelData<MFCellFAB>&  a_phi,
			       DataIterator&                a_dit,
			       bool                         a_homogeneous){
  CH_TIME("MfHelmholtzOp::applyOp");
#if verb
  pout() << "MfHelmholtzOp::applyop"<< endl;
#endif

  this->update_bc(a_phi, a_dit, a_homogeneous);

  for (int i=0; i < m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_lhs);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_phi);
    
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneous, true, a_dit);
  }
}

void MfHelmholtzOp::create(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("MfHelmholtzOp::create");
#if verb
  pout() << "MfHelmholtzOp::create"<< endl;
#endif
  m_ops.create(a_lhs, a_rhs);
}

void MfHelmholtzOp::createCoarsened(LevelData<MFCellFAB>&       a_lhs,
				       const LevelData<MFCellFAB>& a_rhs,
				       const int&                  a_refRat) {
  CH_TIME("MfHelmholtzOp::createCoarsened");
#if verb
  pout() << "MfHelmholtzOp::createCoarsened"<< endl;
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
  pout() << "MfHelmholtzOp::createCoarsened - done"<< endl;
#endif
}

void MfHelmholtzOp::assign(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "MfHelmholtzOp::assign"<< endl;
#endif
  m_ops.assign(a_lhs, a_rhs);
#if verb
  pout() << "MfHelmholtzOp::assign - doen"<< endl;
#endif
}


Real MfHelmholtzOp::dotProduct(const LevelData<MFCellFAB>& a_data1,
				  const LevelData<MFCellFAB>& a_data2){
#if verb
  pout() << "MfHelmholtzOp::dotproduct"<< endl;
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

void MfHelmholtzOp::incr(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_rhs,
			    Real                        a_scale){
#if verb
  pout() << "MfHelmholtzOp::incr"<< endl;
#endif
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
#if verb
  pout() << "MfHelmholtzOp::incr - done"<< endl;
#endif
}

void MfHelmholtzOp::axby(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_x,
			    const LevelData<MFCellFAB>& a_y,
			    Real a,
			    Real b){
#if verb
  pout() << "MfHelmholtzOp::axby"<< endl;
#endif

  m_ops.axby(a_lhs, a_x, a_y, a, b);
}

void MfHelmholtzOp::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale){
  CH_TIME("MfHelmholtzOp::scale");
#if verb
  pout() << "MfHelmholtzOp::scale"<< endl;
#endif
  m_ops.scale(a_lhs, a_scale);
#if verb
  pout() << "MfHelmholtzOp::scale - doen"<< endl;
#endif
}

Real MfHelmholtzOp::norm(const LevelData<MFCellFAB>& a_x, int a_ord){
  CH_TIME("MfHelmholtzOp::norm");
#if verb
  pout() << "MfHelmholtzOp::norm"<< endl;
#endif
  Real volume;
  //  Real rtn = this->kappaNorm(volume, a_x, a_ord);
#if verb
  pout() << "MfHelmholtzOp::norm - done"<< endl;
#endif

#if 1
  Real rtn = 0.0;
  for (int iphase=0; iphase<m_phases;iphase++){
    LevelData<EBCellFAB> alias;
    MultifluidAlias::aliasMF(alias, iphase, a_x);

    Real phaseNorm = m_ebops[iphase]->norm(alias, 0);
    rtn = Max(rtn, phaseNorm);
  }
#endif

  return rtn;
}

Real MfHelmholtzOp::kappaNorm(Real& a_volume, const LevelData<MFCellFAB>& a_data, int a_p) const {
#if verb
  pout() << "MfHelmholtzOp::kappaNorm"<< endl;
#endif
  
  Real accum = 0.0;
  
  a_volume = 0.0;
  int ncomp=a_data.nComp();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit){
    DataIndex d = dit();
    Vector<Real> cur(ncomp,0.0);
    Real curVolume = 0.0;
    Vector<Real> phaseCur(ncomp,0.0);
    Real phaseVolume;

    
    for (int i=0; i<m_phases; i++){
      const EBCellFAB& data = a_data[d].getPhase(i);
      const Box& box = a_data.getBoxes().get(d);
      phaseCur = EBLevelDataOps::vectorSumKappaPow(phaseVolume,data,box,
						   EBLEVELDATAOPS_ALLVOFS,m_domain,a_p);


      for (int i=0; i<ncomp; i++){
	if (a_p == 0){
	  if(phaseCur[i] > cur[i]){
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
      if(a_p == 0){
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
  if(a_p == 0){
    Real recv;
    int result;
    
    result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    accum = recv;
  }
  else{
    Real recv;
    int result;
    
    result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    accum = recv;
    
    result = MPI_Allreduce(&a_volume, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    a_volume = recv;
  }
#endif
  
  if(a_p != 0){
    if (a_volume > 0.0){
      accum = accum / a_volume;
    }
    else{
      accum = 0;
    }
    Real invPow = 1.0/a_p;
    accum = pow(accum,invPow);
  }

  return accum;
}


void MfHelmholtzOp::setToZero(LevelData<MFCellFAB>& a_x){
  CH_TIME("MfHelmholtzOp::setToZero");
#if verb
  pout() << "MfHelmholtzOp::setToZero"<< endl;
#endif
  m_ops.setToZero(a_x);
#if verb
  pout() << "MfHelmholtzOp::setToZero - done"<< endl;
#endif
}

void MfHelmholtzOp::relax(LevelData<MFCellFAB>&       a_e,
			     const LevelData<MFCellFAB>& a_residual,
			     int                         a_iterations){
  CH_TIME("MfHelmholtzOp::relax");
#if verb
  pout() << "MfHelmholtzOp::relax"<< endl;
#endif

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
  //       MayDay::Error("MfHelmholtzOp: Invalid relaxation type");
  //     }
  //   }
  // #else
  //   const bool homogeneous = true;
  //   for (int i = 0; i < a_iterations; i++){
  //     this->update_bc(a_e, homogeneous);

  //     for (int iphase = 0; iphase < m_phases; iphase++){
  //       MultifluidAlias::aliasMF(*m_alias[0], iphase, a_e);
  //       MultifluidAlias::aliasMF(*m_alias[1], iphase, a_residual);

  //       m_ebops[iphase]->lazyGauSai(*m_alias[0], *m_alias[1]);
  //     }
  //   }
  // #endif

  const bool homogeneous = true;
  
  if(!m_multifluid){ // Single-fluid code
    for (int iphase = 0; iphase < m_phases; iphase++){
      MultifluidAlias::aliasMF(*m_alias[0], iphase, a_e);
      MultifluidAlias::aliasMF(*m_alias[1], iphase, a_residual);

      for (int i = 0; i < a_iterations; i++){
	this->update_bc(a_e, homogeneous);
	if(m_relax == 0){
	  m_ebops[iphase]->relaxPoiJac(*m_alias[0], *m_alias[1], 1);

	}
	else if(m_relax == 1){
	  m_ebops[iphase]->relaxGauSai(*m_alias[0], *m_alias[1], 1);
	}
	else if(m_relax == 2){
	  m_ebops[iphase]->relaxGSRBFast(*m_alias[0], *m_alias[1], 1);
	}
	else {
	  MayDay::Abort("MfHelmholtzOp::relax - unknown relaxation type requested");
	}
      }
    }
  }
  else { // Multifluid code
#if 0 // Original code
    for (int i = 0; i < a_iterations; i++){
      this->update_bc(a_e, homogeneous);

      for (int iphase = 0; iphase <= 1; iphase++){
	MultifluidAlias::aliasMF(*m_alias[0], iphase, a_e);
	MultifluidAlias::aliasMF(*m_alias[1], iphase, a_residual);

#if 0 // Original code
	//	m_ebops[iphase]->lazyGauSai(*m_alias[0], *m_alias[1]);
	m_ebops[iphase]->relaxGauSai(*m_alias[0], *m_alias[1], 1);
#else // test
	m_ebops[iphase]->relaxGSRBFast(*m_alias[0], *m_alias[1], 1);
#endif
      }
    }
#else // Place to put optimized code
    const DisjointBoxLayout& dbl = a_e.disjointBoxLayout();
    
    LevelData<MFCellFAB> lphi;
    this->create(lphi, a_residual);
    
    const bool homogeneous = true;

    MultifluidAlias::aliasMF(*m_alias[0], 0, a_e);        
    MultifluidAlias::aliasMF(*m_alias[1], 0, a_residual);
    MultifluidAlias::aliasMF(*m_alias[2], 0, lphi);
    MultifluidAlias::aliasMF(*m_alias[3], 1, a_e);
    MultifluidAlias::aliasMF(*m_alias[4], 1, a_residual);
    MultifluidAlias::aliasMF(*m_alias[5], 1, lphi);

    
    for (int i = 0; i < a_iterations; i++){
      //      const Real t0 = MPI_Wtime();
      this->update_bc(a_e, homogeneous);
      //      const Real t1 = MPI_Wtime();


      for (int icolor = 0; icolor < m_colors.size(); icolor++){
#if 0
	// Get coarse-fine boundary conditions
	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  for (int dir = 0; dir < SpaceDim; dir++){
	    for (SideIterator sit; sit.ok(); sit.next()){
	      m_ebops[0]->applyHomogeneousCFBCs((*m_alias[0])[dit()], dit(), dir, sit());
	      m_ebops[1]->applyHomogeneousCFBCs((*m_alias[3])[dit()], dit(), dir, sit());
	    }
	  }
	}
#endif

#if verb
	pout() << "mfconducitivyop::relax - applying operator gas" << endl;
#endif
	m_ebops[0]->applyOp(*m_alias[2], *m_alias[0], NULL, true, true, false);
#if verb
	pout() << "mfconducitivyop::relax - applying operator solid" << endl;
#endif
	m_ebops[1]->applyOp(*m_alias[5], *m_alias[3], NULL, true, true, false);


	// Intersection boxes will be relaxed twice. Can't avoid this. 
	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  m_ebops[0]->gsrbColor((*m_alias[0])[dit()],
				(*m_alias[2])[dit()],
				(*m_alias[1])[dit()],
				dbl.get(dit()),
				dit(),
				m_colors[icolor]);
	  m_ebops[1]->gsrbColor((*m_alias[3])[dit()],
				(*m_alias[5])[dit()],
				(*m_alias[4])[dit()],
				dbl.get(dit()),
				dit(),
				m_colors[icolor]);
	}
#if 1
	if((icolor-1) % 2 == 0 && icolor - 1 < m_colors.size()){
	  m_alias[0]->exchange();
	  m_alias[3]->exchange();
	}
#endif
      }
      //      const Real t2 = MPI_Wtime();
      //      pout() << "update bc time = " << 100.*(t1 - t0)/(t2-t0) << "\t relax time = " << 100.*(t2 - t1)/(t2-t0) << endl;
    }
#endif
  }

#if verb
  pout() << "MfHelmholtzOp::relax - done"<< endl;
#endif
}

void MfHelmholtzOp::levelJacobi(LevelData<MFCellFAB>&       a_phi,
				   const LevelData<MFCellFAB>& a_rhs,
				   const int                   a_iterations){
  CH_TIME("MfHelmholtzOp::levelJacobi");
  bool homogeneous = true;

  for (int iphase = 0; iphase < m_phases; iphase++){
    MultifluidAlias::aliasMF(*m_alias[0], iphase, a_phi);
    MultifluidAlias::aliasMF(*m_alias[1], iphase, a_rhs);

#if 1 // This is the only way we can make it converge for now
      // m_ebops[iphase]->relaxPoiJac(*m_alias[0], *m_alias[1], 1);
    m_ebops[iphase]->relaxGauSai(*m_alias[0], *m_alias[1], 1);
    //m_ebops[iphase]->lazyGauSai(*m_alias[0], *m_alias[1]);
    /// m_ebops[iphase]->relaxGSRBFast(*m_alias[0], *m_alias[1], 1);
#else
    m_ebops[iphase]->relax_mf(*m_alias[0], *m_alias[1], 1);
#endif
  }
}

void MfHelmholtzOp::createCoarser(LevelData<MFCellFAB>&       a_coarse,
				     const LevelData<MFCellFAB>& a_fine,
				     bool                        ghosted){
  CH_TIME("MfHelmholtzOp::createCoarser");
#if verb
  pout() << "MfHelmholtzOp::createCoarser"<< endl;
#endif

  DisjointBoxLayout dbl = m_mflg_coar_mg.getGrids();

  Vector<int> comps(m_phases, m_ncomp);
  Vector<EBISLayout> ebisl(m_phases);
  for (int iphase = 0; iphase < m_phases; iphase++){
    ebisl[iphase] = m_mflg_coar_mg.getEBLevelGrid(iphase).getEBISL();
  }

  MFCellFactory factory(ebisl, comps);
  a_coarse.define(dbl, m_ncomp, a_fine.ghostVect(), factory);
}

void MfHelmholtzOp::restrictResidual(LevelData<MFCellFAB>&       a_resCoarse,
					LevelData<MFCellFAB>&       a_phiFine,
					const LevelData<MFCellFAB>& a_rhsFine){
  CH_TIME("MfHelmholtzOp::restrictResidual");
#if verb
  pout() << "MfHelmholtzOp::restrictResidual"<< endl;
#endif
  for (int i=0; i < m_phases; i++) {
    MultifluidAlias::aliasMF(*m_alias[0], i, a_resCoarse);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_phiFine);
    MultifluidAlias::aliasMF(*m_alias[2], i, a_rhsFine);
    
    m_ebops[i]->restrictResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void MfHelmholtzOp::prolongIncrement(LevelData<MFCellFAB>&       a_phiThisLevel,
					const LevelData<MFCellFAB>& a_correctCoarse){
  CH_TIME("MfHelmholtzOp::prolongIncrement");
#if verb
  pout() << "MfHelmholtzOp::prolongIncrement"<< endl;
#endif
  for (int i=0; i < m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_phiThisLevel);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_correctCoarse);
    
    m_ebops[i]->prolongIncrement(*m_alias[0], *m_alias[1]);
  }
}

void MfHelmholtzOp::AMRResidual(LevelData<MFCellFAB>&       a_residual,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   const LevelData<MFCellFAB>& a_rhs,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("MfHelmholtzOp::AMRResidual");
#if verb
  pout() << "MfHelmholtzOp::amrresidual"<< endl;
#endif

  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void MfHelmholtzOp::AMROperator(LevelData<MFCellFAB>&       a_LofPhi,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("MfHelmholtzOp::AMROperator");
#if verb
  pout() << "MfHelmholtzOp::amroperator" << endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);

  for (int iphase = 0; iphase < m_phases; iphase++){
    MultifluidAlias::aliasMF(*m_alias[0], iphase, a_LofPhi);
    MultifluidAlias::aliasMF(*m_alias[1], iphase, a_phiFine);
    MultifluidAlias::aliasMF(*m_alias[2], iphase, a_phi);
    MultifluidAlias::aliasMF(*m_alias[3], iphase, a_phiCoarse);

    MfHelmholtzOp* finerOp = (MfHelmholtzOp*) a_finerOp;
#if verb
    pout() << "MfHelmholtzOp::AMROperator - apply EbHelmholtzOps" << endl;
#endif
    m_ebops[iphase]->AMROperator(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_homogeneousBC, finerOp->m_ebops[iphase]);
    // m_ebops[iphase]->applyOp(*m_alias[0], *m_alias[1], m_alias[2], a_homogeneousBC, false);
    // m_ebops[iphase]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[iphase]);
#if verb
    pout() << "MfHelmholtzOp::AMROperator - apply EbHelmholtzOps - done" << endl;
#endif
  }

#if verb
  pout() << "MfHelmholtzOp::amroperator - done" << endl;
#endif
}

void MfHelmholtzOp::AMRResidualNC(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("MfHelmholtzOp::AMRResidualNC");
#if verb
  pout() << "MfHelmholtzOp::amrresidualnc"<< endl;
#endif

  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);

#if verb
  pout() << "MfHelmholtzOp::amrresidualnc - done"<< endl;
#endif
}

void MfHelmholtzOp::AMROperatorNC(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("MfHelmholtzOp::AMROperatorNC");
#if verb
  pout() << "MfHelmholtzOp::amroperatornc"<< endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);
  
  for (int i=0; i < m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_LofPhi);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_phiFine);
    MultifluidAlias::aliasMF(*m_alias[2], i, a_phi);

    MfHelmholtzOp* finerOp = (MfHelmholtzOp*) a_finerOp;
    m_ebops[i]->AMROperatorNC(*m_alias[0], *m_alias[1], *m_alias[2], a_homogeneousBC, finerOp->m_ebops[i]);

    // m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneousBC, true);

    // MfHelmholtzOp* finerOp = (MfHelmholtzOp*) a_finerOp;
    // m_ebops[i]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[i]);
  }
#if verb
  pout() << "MfHelmholtzOp::amroperatornc - done"<< endl;
#endif
}


void MfHelmholtzOp::AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC){
  CH_TIME("MfHelmholtzOp::AMRResidualNF");
#if verb
  pout() << "MfHelmholtzOp::amrresidualnf"<< endl;
#endif

  this->AMROperatorNF(a_residual, a_phi, a_phiCoarse, a_homogeneousBC);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void MfHelmholtzOp::AMROperatorNF(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     bool                        a_homogeneousBC){
  CH_TIME("MfHelmholtzOp::AMROperatorNF");
#if verb
  pout() << "MfHelmholtzOp::amroperatornf"<< endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);
  
  for (int i=0; i<m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_LofPhi);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_phi);
    MultifluidAlias::aliasMF(*m_alias[2], i, a_phiCoarse);

    m_ebops[i]->AMROperatorNF(*m_alias[0], *m_alias[1], *m_alias[2], a_homogeneousBC);
    //    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], m_alias[2], a_homogeneousBC, false);
  }
}



void MfHelmholtzOp::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
					 const LevelData<MFCellFAB>& a_correction,
					 const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("MfHelmholtzOp::AMRUpdateResidual");
#if verb
  pout() << "MfHelmholtzOp::amrupdateresidual"<< endl;
#endif

  this->update_bc(a_correction, true);

  for (int i=0; i < m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_residual);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_correction);
    MultifluidAlias::aliasMF(*m_alias[2], i, a_coarseCorrection);
    
    m_ebops[i]->AMRUpdateResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void MfHelmholtzOp::AMRRestrict(LevelData<MFCellFAB>&       a_resCoarse,
				   const LevelData<MFCellFAB>& a_residual,
				   const LevelData<MFCellFAB>& a_correction,
				   const LevelData<MFCellFAB>& a_coarseCorrection,
				   bool                        a_skip_res){
  CH_TIME("MfHelmholtzOp::AMRRestrict");
#if verb
  pout() << "MfHelmholtzOp::amrrestrict"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_resCoarse);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_residual);
    MultifluidAlias::aliasMF(*m_alias[2], i, a_correction);
    MultifluidAlias::aliasMF(*m_alias[3], i, a_coarseCorrection);
    
    m_ebops[i]->AMRRestrict(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_skip_res);
  }
}

void MfHelmholtzOp::AMRProlong(LevelData<MFCellFAB>&       a_correction,
				  const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("MfHelmholtzOp::AMRProlong");
#if verb
  pout() << "MfHelmholtzOp::amrprolong"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_correction);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_coarseCorrection);
    
    m_ebops[i]->AMRProlong(*m_alias[0], *m_alias[1]);
  }
}



Real MfHelmholtzOp::AMRNorm(const LevelData<MFCellFAB>& a_coar_resid,
			       const LevelData<MFCellFAB>& a_fine_resid,
			       const int&                  a_ref_rat,
			       const int&                  a_ord){
  CH_TIME("mf_helmholtzop::AMRNorm");
#if verb
  pout() << "MfHelmholtzOp::amrnorm"<< endl;
#endif
  Real m = 0;

  for (int i = 0; i < m_phases; i++){
    MultifluidAlias::aliasMF(*m_alias[0], i, a_coar_resid);
    MultifluidAlias::aliasMF(*m_alias[1], i, a_fine_resid);

    Real norm = m_ebops[i]->AMRNorm(*m_alias[0], *m_alias[1], a_ref_rat, a_ord);

    m = Max(m, norm);
  }

  return m;
}

int MfHelmholtzOp::refToCoarser(){
#if verb
  pout() << "MfHelmholtzOp::reftocoarser"<< endl;
#endif
  return m_ref_to_coarser;
}

#include <CD_NamespaceFooter.H>
