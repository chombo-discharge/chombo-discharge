/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MfConductivityOp.cpp
  @brief  Implementation of CD_MfConductivityOp.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <DirichletConductivityDomainBC.H>
#include <MFLevelDataOps.H>
#include <BaseIVFactory.H>
#include <EBAMRDataOps.H>

// Our includes
#include <mfdirichletconductivityebbc.H>
#include <CD_MfConductivityOp.H>
#include <mfalias.H>
#include <data_ops.H>
#include <CD_NamespaceHeader.H>
  
#define verb 0

MfConductivityOp::MfConductivityOp(){
#if verb
  pout() << "MfConductivityOp::MfConductivityOp" << endl;
#endif
}
  

MfConductivityOp::~MfConductivityOp(){
  for (int i = 0; i < m_alias.size(); i++){
    delete m_alias[i];
    m_alias[i] = NULL;
  }

  //  delete m_jumpbc;

  for (int iphase= 0; iphase < m_phases; iphase++){
    //    delete m_ebbc[iphase];
  }
}

void MfConductivityOp::define(const RefCountedPtr<mfis>&                    a_mfis,
			      const RefCountedPtr<BaseDomainBCFactory>&     a_dombc,
			      const RefCountedPtr<LevelData<MFCellFAB> >&   a_aco,
			      const RefCountedPtr<LevelData<MFFluxFAB> >&   a_bco,
			      const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_bco_irreg,
			      const MFQuadCFInterp&                         a_quadcfi,
			      const MFFastFluxReg&                          a_fluxreg,
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


  const int num_phases = a_mfis->num_phases();
  
  const int num_alias  = 6;

  m_multifluidIndexSpace = a_mfis;
  m_ncomp = 1;
  m_relax = a_relax_type;
  m_domain = a_domain;
  m_phases = num_phases;
  m_ebops.resize(num_phases);
  m_ebbc.resize(num_phases);
  m_aCoefficienteffs.resize(num_phases);
  m_bcoeffs.resize(num_phases);
  m_bcoeffs_irr.resize(num_phases);
  m_alias.resize(num_alias);
  m_ref_to_coarser = a_ref_to_coar;
  m_ghost_phi = a_ghost_phi;
  m_ghost_rhs = a_ghost_rhs;
  m_origin = a_origin;
  m_dx = a_dx;
  m_dirival.resize(num_phases);
  m_multifluid = m_multifluidIndexSpace->num_phases() > 1;
  m_lengthScale = a_lengthScale;

  if(a_has_mg){
    m_mflg_coar_mg = a_mflg_coar_mg;
  }

  for (int i = 0; i < num_alias; i++){
    m_alias[i] = new LevelData<EBCellFAB>();
  }

  // Object for matching boundary conditions. Native EBBC is data-based Dirichlet
  m_jumpbc = RefCountedPtr<jump_bc> (new jump_bc(a_mflg, *a_bco_irreg, a_dx, a_order_ebbc, (a_mflg.getEBLevelGrid(0)).getCFIVS()));

  Vector<EBISLayout> layouts(num_phases);
  Vector<int> comps(num_phases);;


  
  for (int iphase = 0; iphase < num_phases; iphase++){

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
      quadcfi = a_quadcfi.getNWOEBQuadCFInterp_ptr(iphase);
      CH_assert(!quadcfi.isNull());
    }
    if(a_has_fine){
      fastFR = a_fluxreg.get_fastfr_ptr(iphase);
    }


    // Coefficients
    m_aCoefficienteffs[iphase]     = RefCountedPtr<LevelData<EBCellFAB> >        (new LevelData<EBCellFAB>());
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
    m_ebbc[iphase] = RefCountedPtr<mfdirichletconductivityebbc> (new mfdirichletconductivityebbc(a_domain,
												 ebisl,
												 a_dx*RealVect::Unit,
												 &a_ghost_phi,
												 &a_ghost_rhs,
												 iphase));
    m_ebbc[iphase]->set_jump_object(m_jumpbc);
    m_ebbc[iphase]->setOrder(a_order_ebbc);
    m_ebbc[iphase]->define_ivs(a_mflg);

    EBLevelDataOps::setVal(*m_dirival[iphase], 0.0);
    m_ebbc[iphase]->setData(m_dirival[iphase]);

    mfalias::aliasMF(*m_aCoefficienteffs[iphase],     iphase, *a_aco);
    mfalias::aliasMF(*m_bcoeffs[iphase],     iphase, *a_bco);
    mfalias::aliasMF(*m_bcoeffs_irr[iphase], iphase, *a_bco_irreg);

    const RefCountedPtr<LevelData<EBCellFAB> >&        aco     = m_aCoefficienteffs[iphase];
    const RefCountedPtr<LevelData<EBFluxFAB> >&        bco     = m_bcoeffs[iphase];
    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& bco_irr = m_bcoeffs_irr[iphase];
    const RefCountedPtr<mfdirichletconductivityebbc>&  ebbc    = m_ebbc[iphase];

    const Real alpha = a_alpha;
    const Real beta  = a_beta;

#if verb
    pout() << "MfConductivityOp::creating oper" << endl;
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
    pout() << "MfConductivityOp::done creating oper" << endl;
#endif
  }

  MFCellFactory* factory = new MFCellFactory(layouts, comps);
  RefCountedPtr<DataFactory<MFCellFAB> > fac(factory);
  m_ops.define(fac);

  EBArith::getMultiColors(m_colors);

#if 0 // Define aggregate stencils in jump_bc object. Only for the new jump_bc class
  LevelData<MFCellFAB> dummy(a_aco->disjointBoxLayout(), 1, a_ghost_phi*IntVect::Unit, *factory);
  for (int iphase = 0; iphase < num_phases; iphase++){
    m_jumpbc->define_agg_stencils(*m_dirival[iphase], dummy, iphase);
  }
#endif
}

void MfConductivityOp::set_jump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump){
  m_jump = a_jump;
}

void MfConductivityOp::setTime(Real* a_time){
  m_time = a_time;
}

void MfConductivityOp::setDirichletEbBc(const ElectrostaticEbBc& a_ebbc){
  m_electrostaticEbBc = a_ebbc;

  this->set_bc_from_levelset();
}

void MfConductivityOp::update_bc(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
  CH_TIME("MfConductivityOp::update_bc");
#if verb
  pout() << "MfConductivityOp::update_bc"<< endl;
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
  pout() << "MfConductivityOp::update_bc - done" << endl;
#endif
}

void MfConductivityOp::update_bc(const LevelData<MFCellFAB>& a_phi, DataIterator& a_dit, const bool a_homogeneous){
  LevelData<MFCellFAB>* phi = const_cast<LevelData<MFCellFAB>* > (&a_phi);
  phi->exchange();
  this->set_bc_from_matching(a_phi, a_dit, a_homogeneous);
}

void MfConductivityOp::set_bc_from_levelset(){
  CH_TIME("MfConductivityOp::set_bc_from_levelset");
#if verb
  pout() << "MfConductivityOp::set_bc_from_levelset"<< endl;
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

	const std::vector<std::pair<electrode, ElectrostaticEbBc::BcFunction> >& electrodeBcs = m_electrostaticEbBc.getBcs();
	
	if(electrodeBcs.size() > 0){
	  int  func = -1;
	  Real dist = std::numeric_limits<Real>::infinity();

	  for (int i = 0; i < electrodeBcs.size(); i++){
	    Real cur_val = (electrodeBcs[i].first.get_function())->value(pos);
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

void MfConductivityOp::set_bc_from_matching(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneous){
  DataIterator dit = a_phi.dataIterator();
  this->set_bc_from_matching(a_phi, dit, a_homogeneous);
}

void MfConductivityOp::set_bc_from_matching(const LevelData<MFCellFAB>& a_phi, DataIterator& a_dit, const bool a_homogeneous){
  CH_TIME("MfConductivityOp::set_bc_from_matching");
#if verb
  pout() << "MfConductivityOp::set_bc_from_matching"<< endl;
#endif

  if(m_multifluid){
    for (int iphase = 0; iphase <= 1; iphase++){
      m_jumpbc->match_bc(*m_dirival[iphase], *m_jump, a_phi, a_dit, a_homogeneous);
    }
  }

#if verb
  pout() << "MfConductivityOp::set_bc_from_matching - done"<< endl;
#endif
}

void MfConductivityOp::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta){
  CH_TIME("MfConductivityOp::setAlphaAndBeta");
#if verb
  pout() << "MfConductivityOp::setalaphandbeta"<< endl;
#endif
  for (int iphase = 0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setAlphaAndBeta(a_alpha, a_beta);
  }
}

void MfConductivityOp::diagonalScale(LevelData<MFCellFAB>& a_rhs){
  CH_TIME("MfConductivityOp::diagonalScale");
#if verb
  pout() << "MfConductivityOp::setdiagonalscale"<< endl;
#endif
  MayDay::Abort("MfConductivityOp::DiagonalScale - where did i get called?");

  // // Operator diagonal scale
  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->diagonalScale(*m_alias[0], true);
  }

  MayDay::Abort("MfConductivityOp::diagonalScale - add scaling for area fractions");
  //  MFLevelDataOps::kappaWeight(a_rhs); // This came from MFPoissonOp
}

void MfConductivityOp::divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "MfConductivityOp::dividebyidenticyoef"<< endl;
#endif

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_rhs);

    m_ebops[iphase]->divideByIdentityCoef(*m_alias[0]);
  }
}

void MfConductivityOp::applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
					 const LevelData<MFCellFAB>& a_phi){
#if verb
  pout() << "MfConductivityOp::applyopnoboundary"<< endl;
#endif
  CH_TIME("MfConductivityOp::applyOpNoBoundary");

  this->update_bc(a_phi, true);
  for (int iphase=0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_opPhi);
    mfalias::aliasMF(*m_alias[1], iphase, a_phi);
    m_ebops[iphase]->applyOpNoBoundary(*m_alias[0], *m_alias[1]);
  }
}

void MfConductivityOp::setTime(Real a_oldTime, Real a_mu, Real a_dt){
  CH_TIME("MfConductivityOp::setTime");
#if verb
  pout() << "MfConductivityOp::settime"<< endl;
#endif
  for (int iphase=0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setTime(a_oldTime, a_mu, a_dt);
  }
}

void MfConductivityOp::residual(LevelData<MFCellFAB>&        a_lhs,
				const LevelData<MFCellFAB>&  a_phi,
				const LevelData<MFCellFAB>&  a_rhs,
				bool                         a_homogeneous){
  CH_TIME("MfConductivityOp::residual");
#if verb
  pout() << "MfConductivityOp::residual"<< endl;
#endif
  
  this->applyOp(a_lhs, a_phi, a_homogeneous);

  // Make residual rhs - L(phi)
  this->scale(a_lhs,      -1.0);
  this->incr(a_lhs, a_rhs, 1.0);

}

void MfConductivityOp::preCond(LevelData<MFCellFAB>&       a_correction,
			       const LevelData<MFCellFAB>& a_residual){
  CH_TIME("MfConductivityOp::preCond");
#if verb
  pout() << "MfConductivityOp::precond"<< endl;
#endif

  this->relax(a_correction, a_residual, 40);
}

void MfConductivityOp::applyOp(LevelData<MFCellFAB>&        a_lhs,
			       const LevelData<MFCellFAB>&  a_phi,
			       bool                         a_homogeneous){
  CH_TIME("MfConductivityOp::applyOp");
#if verb
  pout() << "MfConductivityOp::applyop"<< endl;
#endif

  DataIterator dit = a_lhs.dataIterator();
  this->applyOp(a_lhs, a_phi, dit, a_homogeneous);

}

void MfConductivityOp::applyOp(LevelData<MFCellFAB>&        a_lhs,
			       const LevelData<MFCellFAB>&  a_phi,
			       DataIterator&                a_dit,
			       bool                         a_homogeneous){
  CH_TIME("MfConductivityOp::applyOp");
#if verb
  pout() << "MfConductivityOp::applyop"<< endl;
#endif

  this->update_bc(a_phi, a_dit, a_homogeneous);

  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_lhs);
    mfalias::aliasMF(*m_alias[1], i, a_phi);
    
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneous, true, a_dit);
  }
}

void MfConductivityOp::create(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("MfConductivityOp::create");
#if verb
  pout() << "MfConductivityOp::create"<< endl;
#endif
  m_ops.create(a_lhs, a_rhs);
}

void MfConductivityOp::createCoarsened(LevelData<MFCellFAB>&       a_lhs,
				       const LevelData<MFCellFAB>& a_rhs,
				       const int&                  a_refRat) {
  CH_TIME("MfConductivityOp::createCoarsened");
#if verb
  pout() << "MfConductivityOp::createCoarsened"<< endl;
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
  pout() << "MfConductivityOp::createCoarsened - done"<< endl;
#endif
}

void MfConductivityOp::assign(LevelData<MFCellFAB>&       a_lhs,
			      const LevelData<MFCellFAB>& a_rhs){
#if verb
  pout() << "MfConductivityOp::assign"<< endl;
#endif
  m_ops.assign(a_lhs, a_rhs);
#if verb
  pout() << "MfConductivityOp::assign - doen"<< endl;
#endif
}


Real MfConductivityOp::dotProduct(const LevelData<MFCellFAB>& a_data1,
				  const LevelData<MFCellFAB>& a_data2){
#if verb
  pout() << "MfConductivityOp::dotproduct"<< endl;
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

void MfConductivityOp::incr(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_rhs,
			    Real                        a_scale){
#if verb
  pout() << "MfConductivityOp::incr"<< endl;
#endif
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
#if verb
  pout() << "MfConductivityOp::incr - done"<< endl;
#endif
}

void MfConductivityOp::axby(LevelData<MFCellFAB>&       a_lhs,
			    const LevelData<MFCellFAB>& a_x,
			    const LevelData<MFCellFAB>& a_y,
			    Real a,
			    Real b){
#if verb
  pout() << "MfConductivityOp::axby"<< endl;
#endif

  m_ops.axby(a_lhs, a_x, a_y, a, b);
}

void MfConductivityOp::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale){
  CH_TIME("MfConductivityOp::scale");
#if verb
  pout() << "MfConductivityOp::scale"<< endl;
#endif
  m_ops.scale(a_lhs, a_scale);
#if verb
  pout() << "MfConductivityOp::scale - doen"<< endl;
#endif
}

Real MfConductivityOp::norm(const LevelData<MFCellFAB>& a_x, int a_ord){
  CH_TIME("MfConductivityOp::norm");
#if verb
  pout() << "MfConductivityOp::norm"<< endl;
#endif
  Real volume;
  //  Real rtn = this->kappaNorm(volume, a_x, a_ord);
#if verb
  pout() << "MfConductivityOp::norm - done"<< endl;
#endif

#if 1
  Real rtn = 0.0;
  for (int iphase=0; iphase<m_phases;iphase++){
    LevelData<EBCellFAB> alias;
    mfalias::aliasMF(alias, iphase, a_x);

    Real phaseNorm = m_ebops[iphase]->norm(alias, 0);
    rtn = Max(rtn, phaseNorm);
  }
#endif

  return rtn;
}

Real MfConductivityOp::kappaNorm(Real& a_volume, const LevelData<MFCellFAB>& a_data, int a_p) const {
#if verb
  pout() << "MfConductivityOp::kappaNorm"<< endl;
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


void MfConductivityOp::setToZero(LevelData<MFCellFAB>& a_x){
  CH_TIME("MfConductivityOp::setToZero");
#if verb
  pout() << "MfConductivityOp::setToZero"<< endl;
#endif
  m_ops.setToZero(a_x);
#if verb
  pout() << "MfConductivityOp::setToZero - done"<< endl;
#endif
}

void MfConductivityOp::relax(LevelData<MFCellFAB>&       a_e,
			     const LevelData<MFCellFAB>& a_residual,
			     int                         a_iterations){
  CH_TIME("MfConductivityOp::relax");
#if verb
  pout() << "MfConductivityOp::relax"<< endl;
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
  //       MayDay::Error("MfConductivityOp: Invalid relaxation type");
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
	  MayDay::Abort("MfConductivityOp::relax - unknown relaxation type requested");
	}
      }
    }
  }
  else { // Multifluid code
#if 0 // Original code
    for (int i = 0; i < a_iterations; i++){
      this->update_bc(a_e, homogeneous);

      for (int iphase = 0; iphase <= 1; iphase++){
	mfalias::aliasMF(*m_alias[0], iphase, a_e);
	mfalias::aliasMF(*m_alias[1], iphase, a_residual);

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

    mfalias::aliasMF(*m_alias[0], 0, a_e);        
    mfalias::aliasMF(*m_alias[1], 0, a_residual);
    mfalias::aliasMF(*m_alias[2], 0, lphi);
    mfalias::aliasMF(*m_alias[3], 1, a_e);
    mfalias::aliasMF(*m_alias[4], 1, a_residual);
    mfalias::aliasMF(*m_alias[5], 1, lphi);

    
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
  pout() << "MfConductivityOp::relax - done"<< endl;
#endif
}

void MfConductivityOp::levelJacobi(LevelData<MFCellFAB>&       a_phi,
				   const LevelData<MFCellFAB>& a_rhs,
				   const int                   a_iterations){
  CH_TIME("MfConductivityOp::levelJacobi");
  bool homogeneous = true;

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_phi);
    mfalias::aliasMF(*m_alias[1], iphase, a_rhs);

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

void MfConductivityOp::createCoarser(LevelData<MFCellFAB>&       a_coarse,
				     const LevelData<MFCellFAB>& a_fine,
				     bool                        ghosted){
  CH_TIME("MfConductivityOp::createCoarser");
#if verb
  pout() << "MfConductivityOp::createCoarser"<< endl;
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

void MfConductivityOp::restrictResidual(LevelData<MFCellFAB>&       a_resCoarse,
					LevelData<MFCellFAB>&       a_phiFine,
					const LevelData<MFCellFAB>& a_rhsFine){
  CH_TIME("MfConductivityOp::restrictResidual");
#if verb
  pout() << "MfConductivityOp::restrictResidual"<< endl;
#endif
  for (int i=0; i < m_phases; i++) {
    mfalias::aliasMF(*m_alias[0], i, a_resCoarse);
    mfalias::aliasMF(*m_alias[1], i, a_phiFine);
    mfalias::aliasMF(*m_alias[2], i, a_rhsFine);
    
    m_ebops[i]->restrictResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void MfConductivityOp::prolongIncrement(LevelData<MFCellFAB>&       a_phiThisLevel,
					const LevelData<MFCellFAB>& a_correctCoarse){
  CH_TIME("MfConductivityOp::prolongIncrement");
#if verb
  pout() << "MfConductivityOp::prolongIncrement"<< endl;
#endif
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_phiThisLevel);
    mfalias::aliasMF(*m_alias[1], i, a_correctCoarse);
    
    m_ebops[i]->prolongIncrement(*m_alias[0], *m_alias[1]);
  }
}

void MfConductivityOp::AMRResidual(LevelData<MFCellFAB>&       a_residual,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   const LevelData<MFCellFAB>& a_rhs,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("MfConductivityOp::AMRResidual");
#if verb
  pout() << "MfConductivityOp::amrresidual"<< endl;
#endif

  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void MfConductivityOp::AMROperator(LevelData<MFCellFAB>&       a_LofPhi,
				   const LevelData<MFCellFAB>& a_phiFine,
				   const LevelData<MFCellFAB>& a_phi,
				   const LevelData<MFCellFAB>& a_phiCoarse,
				   bool                        a_homogeneousBC,
				   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("MfConductivityOp::AMROperator");
#if verb
  pout() << "MfConductivityOp::amroperator" << endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);

  for (int iphase = 0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_alias[0], iphase, a_LofPhi);
    mfalias::aliasMF(*m_alias[1], iphase, a_phiFine);
    mfalias::aliasMF(*m_alias[2], iphase, a_phi);
    mfalias::aliasMF(*m_alias[3], iphase, a_phiCoarse);

    MfConductivityOp* finerOp = (MfConductivityOp*) a_finerOp;
#if verb
    pout() << "MfConductivityOp::AMROperator - apply EbHelmholtzOps" << endl;
#endif
    m_ebops[iphase]->AMROperator(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_homogeneousBC, finerOp->m_ebops[iphase]);
    // m_ebops[iphase]->applyOp(*m_alias[0], *m_alias[1], m_alias[2], a_homogeneousBC, false);
    // m_ebops[iphase]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[iphase]);
#if verb
    pout() << "MfConductivityOp::AMROperator - apply EbHelmholtzOps - done" << endl;
#endif
  }

#if verb
  pout() << "MfConductivityOp::amroperator - done" << endl;
#endif
}

void MfConductivityOp::AMRResidualNC(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp){
  CH_TIME("MfConductivityOp::AMRResidualNC");
#if verb
  pout() << "MfConductivityOp::amrresidualnc"<< endl;
#endif

  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousBC, a_finerOp);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);

#if verb
  pout() << "MfConductivityOp::amrresidualnc - done"<< endl;
#endif
}

void MfConductivityOp::AMROperatorNC(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phiFine,
				     const LevelData<MFCellFAB>& a_phi,
				     bool                        a_homogeneousBC,
				     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("MfConductivityOp::AMROperatorNC");
#if verb
  pout() << "MfConductivityOp::amroperatornc"<< endl;
#endif

  this->update_bc(a_phi, a_homogeneousBC);
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_LofPhi);
    mfalias::aliasMF(*m_alias[1], i, a_phiFine);
    mfalias::aliasMF(*m_alias[2], i, a_phi);

    MfConductivityOp* finerOp = (MfConductivityOp*) a_finerOp;
    m_ebops[i]->AMROperatorNC(*m_alias[0], *m_alias[1], *m_alias[2], a_homogeneousBC, finerOp->m_ebops[i]);

    // m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, a_homogeneousBC, true);

    // MfConductivityOp* finerOp = (MfConductivityOp*) a_finerOp;
    // m_ebops[i]->reflux(*m_alias[0], *m_alias[3], *m_alias[1], finerOp->m_ebops[i]);
  }
#if verb
  pout() << "MfConductivityOp::amroperatornc - done"<< endl;
#endif
}


void MfConductivityOp::AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     const LevelData<MFCellFAB>& a_rhs,
				     bool                        a_homogeneousBC){
  CH_TIME("MfConductivityOp::AMRResidualNF");
#if verb
  pout() << "MfConductivityOp::amrresidualnf"<< endl;
#endif

  this->AMROperatorNF(a_residual, a_phi, a_phiCoarse, a_homogeneousBC);

  // Make residual = a_rhs - L(phi)
  this->scale(a_residual,      -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void MfConductivityOp::AMROperatorNF(LevelData<MFCellFAB>&       a_LofPhi,
				     const LevelData<MFCellFAB>& a_phi,
				     const LevelData<MFCellFAB>& a_phiCoarse,
				     bool                        a_homogeneousBC){
  CH_TIME("MfConductivityOp::AMROperatorNF");
#if verb
  pout() << "MfConductivityOp::amroperatornf"<< endl;
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



void MfConductivityOp::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
					 const LevelData<MFCellFAB>& a_correction,
					 const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("MfConductivityOp::AMRUpdateResidual");
#if verb
  pout() << "MfConductivityOp::amrupdateresidual"<< endl;
#endif

  this->update_bc(a_correction, true);

  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_residual);
    mfalias::aliasMF(*m_alias[1], i, a_correction);
    mfalias::aliasMF(*m_alias[2], i, a_coarseCorrection);
    
    m_ebops[i]->AMRUpdateResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

void MfConductivityOp::AMRRestrict(LevelData<MFCellFAB>&       a_resCoarse,
				   const LevelData<MFCellFAB>& a_residual,
				   const LevelData<MFCellFAB>& a_correction,
				   const LevelData<MFCellFAB>& a_coarseCorrection,
				   bool                        a_skip_res){
  CH_TIME("MfConductivityOp::AMRRestrict");
#if verb
  pout() << "MfConductivityOp::amrrestrict"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_resCoarse);
    mfalias::aliasMF(*m_alias[1], i, a_residual);
    mfalias::aliasMF(*m_alias[2], i, a_correction);
    mfalias::aliasMF(*m_alias[3], i, a_coarseCorrection);
    
    m_ebops[i]->AMRRestrict(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3], a_skip_res);
  }
}

void MfConductivityOp::AMRProlong(LevelData<MFCellFAB>&       a_correction,
				  const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("MfConductivityOp::AMRProlong");
#if verb
  pout() << "MfConductivityOp::amrprolong"<< endl;
#endif
  
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_correction);
    mfalias::aliasMF(*m_alias[1], i, a_coarseCorrection);
    
    m_ebops[i]->AMRProlong(*m_alias[0], *m_alias[1]);
  }
}



Real MfConductivityOp::AMRNorm(const LevelData<MFCellFAB>& a_coar_resid,
			       const LevelData<MFCellFAB>& a_fine_resid,
			       const int&                  a_ref_rat,
			       const int&                  a_ord){
  CH_TIME("mf_helmholtzop::AMRNorm");
#if verb
  pout() << "MfConductivityOp::amrnorm"<< endl;
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

int MfConductivityOp::refToCoarser(){
#if verb
  pout() << "MfConductivityOp::reftocoarser"<< endl;
#endif
  return m_ref_to_coarser;
}

#include <CD_NamespaceFooter.H>
