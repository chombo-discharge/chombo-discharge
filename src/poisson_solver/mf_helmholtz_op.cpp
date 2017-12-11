/*!
  @file mf_helmholtz_op.cpp
  @brief Implementation of mf_helmholtz_op.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "mf_helmholtz_op.H"
#include "mfalias.H"

mf_helmholtz_op::mf_helmholtz_op(){

}
  

mf_helmholtz_op::~mf_helmholtz_op(){
  
}

void mf_helmholtz_op::define(const RefCountedPtr<mfis>&                    a_mfis,
			     const RefCountedPtr<MFQuadCFInterp>&          a_quadcfi,
			     const RefCountedPtr<BaseDomainBC>&            a_dombc,
			     const RefCountedPtr<LevelData<MFCellFAB> >&   a_aco,
			     const RefCountedPtr<LevelData<MFFluxFAB> >&   a_bco,
			     const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_bco_irreg,
			     const MFLevelGrid&                            a_mflg_fine,
			     const MFLevelGrid&                            a_mflg,
			     const MFLevelGrid&                            a_mflg_coar,
			     const MFLevelGrid&                            a_mflg_coar_mg,
			     const DisjointBoxLayout&                      a_dbl,
			     const DisjointBoxLayout&                      a_dbl_finer,
			     const DisjointBoxLayout&                      a_dbl_coarser,
			     const DisjointBoxLayout&                      a_dbl_coar_mg,
			     const ProblemDomain&                          a_domain,
			     const bool&                                   a_layout_changed,
			     const bool&                                   a_has_mg_objects,
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
  MayDay::Warning("Remember to check how many aliasing holders we need");
#endif
  const int num_alias  = 4;
  
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
  
  // Object for matching boundary conditions. Native EBBC is always Dirichlet
  m_jumpbc = RefCountedPtr<jump_bc> (new jump_bc(a_mflg, *a_bco_irreg, a_dx, 2, (a_mflg.get_eblg(0)).getCFIVS()));

  for (int iphase = 0; iphase < num_phases; iphase++){

    const EBLevelGrid& eblg_fine = a_mflg_fine.get_eblg(iphase);
    const EBLevelGrid& eblg      = a_mflg.get_eblg(iphase);
    const EBLevelGrid& eblg_coar = a_mflg_coar.get_eblg(iphase);
    const EBLevelGrid& eblg_mg   = a_mflg_coar_mg.get_eblg(iphase);
    const EBISLayout&  ebisl     = eblg.getEBISL();

    const RefCountedPtr<EBQuadCFInterp>& quadcfi = a_quadcfi->get_quadcfi_ptr(iphase);

    m_acoeffs[iphase]     = RefCountedPtr<LevelData<EBCellFAB> >        (new LevelData<EBCellFAB>());
    m_bcoeffs[iphase]     = RefCountedPtr<LevelData<EBFluxFAB> >        (new LevelData<EBFluxFAB>());
    m_bcoeffs_irr[iphase] = RefCountedPtr<LevelData<BaseIVFAB<Real> > > (new LevelData<BaseIVFAB<Real> >());
    m_ebbc[iphase]        = RefCountedPtr<DirichletConductivityEBBC>    (new DirichletConductivityEBBC(a_domain,
												       ebisl,
												       a_dx*RealVect::Unit,
												       &a_ghost_phi,
												       &a_ghost_rhs));

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
									    a_dombc,
									    ebbc,
									    a_dx,
									    a_dx_coar,
									    a_ref_to_fine,
									    a_ref_to_coar,
									    a_has_fine,
									    a_has_coar,
									    a_has_mg_objects,
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

void mf_helmholtz_op::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta){
  m_alpha = a_alpha;
  m_beta  = a_beta;
}

void mf_helmholtz_op::set_jump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump){
  m_jump = a_jump;
}

void mf_helmholtz_op::set_electrodes(const Vector<electrode>& a_electrodes){
  m_electrodes = a_electrodes;
}

void mf_helmholtz_op::update_bc(){
}

void mf_helmholtz_op::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
					const LevelData<MFCellFAB>& a_correction,
					const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("MFPoissonOp::AMRUpdateResidual");
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
