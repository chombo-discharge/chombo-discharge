/*!
  @file mf_helmholtz_op.cpp
  @brief Implementation of mf_helmholtz_op.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "mf_helmholtz_op.H"
#include "mfalias.H"

#include <DirichletConductivityDomainBC.H>
#include <MFLevelDataOps.H>

mf_helmholtz_op::mf_helmholtz_op(){

}
  

mf_helmholtz_op::~mf_helmholtz_op(){
  
}

void mf_helmholtz_op::define(const RefCountedPtr<mfis>&                    a_mfis,
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
			     const Real&                                   a_beta){


  const int num_phases = a_mfis->num_phases();

#if 1
  MayDay::Warning("mf_helmholtz_op::mf_helmholtz_op - remember to check how many aliasing holders we need");
#endif
  const int num_alias  = 4;

  m_relax = a_relax_type;
  m_domain = a_domain;
  m_phases = num_phases;
  m_ebops.resize(num_phases);
  m_ebbc.resize(num_phases);
  m_acoeffs.resize(num_phases);
  m_bcoeffs.resize(num_phases);
  m_bcoeffs_irr.resize(num_phases);
  m_cell_alias.resize(num_alias);
  m_flux_alias.resize(num_alias);
  m_eb_alias.resize(num_alias);

  for (int i = 0; i < num_alias; i++){
    m_cell_alias[i] = new LevelData<EBCellFAB>();
    m_flux_alias[i] = new LevelData<EBFluxFAB>();
    m_eb_alias[i]   = new LevelData<BaseIVFAB<Real> >();
  }

  // Object for matching boundary conditions. Native EBBC is data-based Dirichlet
  m_jumpbc = RefCountedPtr<jump_bc> (new jump_bc(a_mflg, *a_bco_irreg, a_dx, a_order_ebbc, (a_mflg.get_eblg(0)).getCFIVS()));

  for (int iphase = 0; iphase < num_phases; iphase++){

    const EBLevelGrid& eblg      = a_mflg.get_eblg(iphase);
    const EBISLayout&  ebisl     = eblg.getEBISL();

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

    const RefCountedPtr<EBQuadCFInterp>& quadcfi = a_quadcfi.get_quadcfi_ptr(iphase);
    if(a_has_coar){
      CH_assert(!quadcfi.isNull());
    }

    m_acoeffs[iphase]     = RefCountedPtr<LevelData<EBCellFAB> >        (new LevelData<EBCellFAB>());
    m_bcoeffs[iphase]     = RefCountedPtr<LevelData<EBFluxFAB> >        (new LevelData<EBFluxFAB>());
    m_bcoeffs_irr[iphase] = RefCountedPtr<LevelData<BaseIVFAB<Real> > > (new LevelData<BaseIVFAB<Real> >());
    m_ebbc[iphase]        = RefCountedPtr<DirichletConductivityEBBC>    (new DirichletConductivityEBBC(a_domain,
												       ebisl,
												       a_dx*RealVect::Unit,
												       &a_ghost_phi,
												       &a_ghost_rhs));

#if 1
    m_ebbc[iphase]->setValue(0.0);
    m_ebbc[iphase]->setOrder(2);
    MayDay::Warning("mf_helmholtz_op::mf_helmholtz_op - remember to fix domain boundary conditions");
    RefCountedPtr<DirichletConductivityDomainBC> dombc = RefCountedPtr<DirichletConductivityDomainBC>
      (new DirichletConductivityDomainBC());
    dombc->setValue(0.0);
#endif


    mfalias::aliasMF(*m_acoeffs[iphase],     iphase, *a_aco);
    mfalias::aliasMF(*m_bcoeffs[iphase],     iphase, *a_bco);
    mfalias::aliasMF(*m_bcoeffs_irr[iphase], iphase, *a_bco_irreg);

    const RefCountedPtr<LevelData<EBCellFAB> >&        aco     = m_acoeffs[iphase];
    const RefCountedPtr<LevelData<EBFluxFAB> >&        bco     = m_bcoeffs[iphase];
    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& bco_irr = m_bcoeffs_irr[iphase];
    const RefCountedPtr<DirichletConductivityEBBC>     ebbc    = m_ebbc[iphase];

    m_ebops[iphase] = RefCountedPtr<EBConductivityOp> (new EBConductivityOp(eblg_fine,
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
									    a_alpha,
									    a_beta,
									    aco,
									    bco,
									    bco_irr,
									    a_ghost_phi,
									    a_ghost_rhs,
									    a_relax_type));
									    
    
  }

}



void mf_helmholtz_op::set_jump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump){
  m_jump = a_jump;
}

void mf_helmholtz_op::set_electrodes(const Vector<electrode>& a_electrodes){
  m_electrodes = a_electrodes;
}

void mf_helmholtz_op::update_bc(){
  
}

void mf_helmholtz_op::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta){
  CH_TIME("mf_helmholtz_op::setAlphaAndBeta");
  for (int iphase = 0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setAlphaAndBeta(a_alpha, a_beta);
#if 0
    MayDay::Warning("mf_helmholtz_op::setAlphaAndBeta - MFPoisson multiplies by acoef, why?");
#endif
  }
}

void mf_helmholtz_op::diagonalScale(LevelData<MFCellFAB>& a_rhs){
  CH_TIME("mf_helmholtz_op::diagonalScale");
  MFLevelDataOps::kappaWeight(a_rhs);
}

void mf_helmholtz_op::divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs){
  
}

void mf_helmholtz_op::applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
					const LevelData<MFCellFAB>& a_phi){
  CH_TIME("mf_helmholtz_op::applyOpNoBoundary");
  for (int iphase=0; iphase < m_phases; iphase++){
    mfalias::aliasMF(*m_cell_alias[0], iphase, a_opPhi);
    mfalias::aliasMF(*m_cell_alias[1], iphase, a_phi);
    m_ebops[iphase]->applyOpNoBoundary(*m_cell_alias[0], *m_cell_alias[1]);
  }
}

void mf_helmholtz_op::setTime(Real a_oldTime, Real a_mu, Real a_dt){
  CH_TIME("mf_helmholtz_op::setTime");
  for (int iphase=0; iphase < m_phases; iphase++){
    m_ebops[iphase]->setTime(a_oldTime, a_mu, a_dt);
  }

#if 1
  MayDay::Warning("mf_helmholtz_op::setTime - jump bc should also be time dependent");
#endif
}

void mf_helmholtz_op::residual(LevelData<MFCellFAB>&        a_lhs,
			       const LevelData<MFCellFAB>&  a_phi,
			       const LevelData<MFCellFAB>&  a_rhs,
			       bool                         a_homogeneous){
  CH_TIME("mf_helmholtz_op::residual");
  
  this->applyOp(a_lhs, a_phi, a_homogeneous);
  this->incr(a_lhs, a_rhs, -1);
  this->scale(a_lhs, -1.0);
}

void mf_helmholtz_op::preCond(LevelData<MFCellFAB>&       a_correction,
			      const LevelData<MFCellFAB>& a_residual){
  CH_TIME("mf_helmholtz_op::preCond");
  this->relax(a_correction, a_residual, 40);
}

void mf_helmholtz_op::applyOp(LevelData<MFCellFAB>&        a_lhs,
			      const LevelData<MFCellFAB>&  a_phi,
			      bool                         a_homogeneous){
  CH_TIME("mf_helmholtz_op::applyOp");
#if 1
  MayDay::Warning("mf_helmholtz_op::applyOp - the matching condition should be updated first");
#endif

  for (int i=0; i<m_phases; i++){
    mfalias::aliasMF(*m_cell_alias[0], i, a_lhs);
    mfalias::aliasMF(*m_cell_alias[1], i, a_phi);
#if 0
    MayDay::Abort("mf_helmholtz_op::applyOp - this only works for mf_helmholtz_op. WE should use the Dirichlet EBBC object");
    
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL,
                        a_homogeneous, true, m_boundaryN+i);
#endif
  }
}

void mf_helmholtz_op::applyOp(LevelData<MFCellFAB>&        a_lhs,
			      const LevelData<MFCellFAB>&  a_phi,
			      DataIterator&                a_dit,
			      bool                         a_homogeneous){
  CH_TIME("mf_helmholtz_op::applyOp(dit)");
#if 1
  MayDay::Warning("mf_helmholtz_op::applyOp - the matching condition should be updated first");
#endif

  for (int i=0; i<m_phases; i++){
    mfalias::aliasMF(*m_cell_alias[0], i, a_lhs);
    mfalias::aliasMF(*m_cell_alias[1], i, a_phi);
#if 0
    MayDay::Abort("mf_helmholtz_op::applyOp - this only works for mf_helmholtz_op. WE should use the Dirichlet EBBC object");
    
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL, // What the fuck is this?
			a_dit, a_homogeneous, true,
                        m_boundaryN+i);
#endif
  }
}

void mf_helmholtz_op::create(LevelData<MFCellFAB>&       a_lhs,
			     const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("mf_helmholtz_op::create");
  MayDay::Warning("mf_helmholtz_op::create - m_ops has not been created");
  m_ops.create(a_lhs, a_rhs);
}

void mf_helmholtz_op::createCoarsened(LevelData<MFCellFAB>&       a_lhs,
				      const LevelData<MFCellFAB>& a_rhs,
				      const int&                  a_refRat) {
  CH_TIME("mf_helmholtz_op::createCoarsened");
  MayDay::Warning("mf_helmholtz_op::createCoarsend - not implemented");
#if 0
  IntVect ghostVect = a_rhs.ghostVect();
  Vector<EBLevelGrid> eblg(m_phases);
  for (int i=0; i<m_phases; i++)
    {
      eblg[i] = m_ebops[i]->getEBLG();
    }

  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, eblg[0].getDBL(), a_refRat);
  ProblemDomain coarDom = coarsen(eblg[0].getDomain(), a_refRat);
  IntVect ghostVec = a_rhs.ghostVect();
  Vector<EBISLayout> ebislv(m_phases);
  for (int i=0; i<m_phases; i++)
    {
      eblg[i].getEBIS()->fillEBISLayout(ebislv[i], dblCoarsenedFine, coarDom , ghostVec[0]);
      if (a_refRat > 2)
        {
          ebislv[i].setMaxRefinementRatio(a_refRat, eblg[i].getEBIS());
        }
    }
  //create coarsened data
  Vector<int> vncomp(m_phases, a_rhs.nComp());
  MFCellFactory fact(ebislv, vncomp);
  a_lhs.define(dblCoarsenedFine, a_rhs.nComp(), ghostVec, fact);
#endif
}

void mf_helmholtz_op::assign(LevelData<MFCellFAB>&       a_lhs,
			     const LevelData<MFCellFAB>& a_rhs){

  MayDay::Warning("mf_helmholtz_op::assign - m_ops not defined (yet");
  m_ops.assign(a_lhs, a_rhs);
}


Real mf_helmholtz_op::dotProduct(const LevelData<MFCellFAB>& a_data1,
				 const LevelData<MFCellFAB>& a_data2){
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

void mf_helmholtz_op::incr(LevelData<MFCellFAB>&       a_lhs,
			   const LevelData<MFCellFAB>& a_rhs,
			   Real                        a_scale){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
}

void mf_helmholtz_op::axby(LevelData<MFCellFAB>&       a_lhs,
			   const LevelData<MFCellFAB>& a_x,
			   const LevelData<MFCellFAB>& a_y,
			   Real a,
			   Real b){

  MayDay::Warning("mf_helmholtz_op::axby - m_ops not defined (yet");
  m_ops.axby(a_lhs, a_x, a_y, a, b);
}

void mf_helmholtz_op::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale){
  CH_TIME("mf_helmholtz_op::scale");
  MayDay::Warning("mf_helmholtz_op::scale - m_ops not defined (yet");
  m_ops.scale(a_lhs, a_scale);
}

Real mf_helmholtz_op::norm(const LevelData<MFCellFAB>& a_x, int a_ord){
  CH_TIME("mf_helmholtz_op::norm");
  Real volume;
  Real rtn = kappaNorm(volume, a_x, a_ord);
  return rtn;
}

void mf_helmholtz_op::setToZero(LevelData<MFCellFAB>& a_x){
  CH_TIME("mf_helmholtz_op::setToZero");
  MayDay::Warning("mf_helmholtz_op::scale - m_ops not defined (yet");
  m_ops.setToZero(a_x);
}

void mf_helmholtz_op::relax(LevelData<MFCellFAB>&       a_e,
			    const LevelData<MFCellFAB>& a_residual,
			    int                         iterations){
  CH_TIME("mf_helmholtz_op::relax");
  
  for (int i=0; i<iterations; i++){
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
      MayDay::Error("mf_helmholtz_op: Invalid relaxation type");
    }
  }
}

void mf_helmholtz_op::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
					const LevelData<MFCellFAB>& a_correction,
					const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("mf_helmholtz_op::AMRUpdateResidual");
  bool homogeneousBC = true;
  //  computeBoundaryN(a_correction, homogeneousBC);
  for (int i=0; i < m_phases; i++){
    // mfalias::aliasMF(*m_alias[0], i, a_residual);
    // mfalias::aliasMF(*m_alias[1], i, a_correction);
    // mfalias::aliasMF(*m_alias[2], i, a_coarseCorrection);
    //    m_ebops[i]->AMRUpdateResidual(*m_alias[0], *m_alias[1], *m_alias[2], m_boundaryN+i);
  }
}

Real mf_helmholtz_op::AMRNorm(const LevelData<MFCellFAB>& a_coar_resid,
			      const LevelData<MFCellFAB>& a_fine_resid,
			      const int&                  a_ref_rat,
			      const int&                  a_ord){
  CH_TIME("mf_helmholtzop::AMRNorm");
  Real m = 0;

  for (int i = 0; i < m_phases; i++){
    // mfalias::aliasMF(*m_alias[0], i, a_coar_resid);
    // mfalias::aliasMF(*m_alias[1], i, a_fine_resid);

    // Real norm = m_ebops[i]->AMRNorm(*m_alias[0], *m_alias[1], a_ref_rat, a_ord);
    Real norm;

    m = Max(m, norm);
  }

  return m;
}
