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

void mf_helmholtz_op::define(const RefCountedPtr<mfis>&         a_mfis,
			     const RefCountedPtr<BaseDomainBC>& a_dombc,
			     const RefCountedPtr<jump_bc>&      a_jumpbc,
			     const MFLevelGrid&                 a_mflg_fine,
			     const MFLevelGrid&                 a_mflg,
			     const MFLevelGrid&                 a_mflg_coar,
			     const MFLevelGrid&                 a_mflg_coar_mg,
			     const DisjointBoxLayout&           a_dbl,
			     const DisjointBoxLayout&           a_dbl_finer,
			     const DisjointBoxLayout&           a_dbl_coarser,
			     const DisjointBoxLayout&           a_dbl_coar_mg,
			     const bool&                        a_layout_changed,
			     const bool&                        a_has_mg_objects,
			     const bool&                        a_has_coarser,
			     const bool&                        a_has_finer,
			     const int&                         a_ref_ratio,
			     const int&                         a_ref_ratio_to_finer,
			     const IntVect&                     a_ghost_phi,
			     const IntVect&                     a_ghost_rhs,
			     const Real&                        a_dx,
			     const Real&                        a_alpha,
			     const Real&                        a_beta){

}

void mf_helmholtz_op::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta){
  m_alpha = a_alpha;
  m_beta  = a_beta;
}

void mf_helmholtz_op::AMRUpdateResidual(LevelData<MFCellFAB>&      a_residual,
					const LevelData<MFCellFAB>& a_correction,
					const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("MFPoissonOp::AMRUpdateResidual");
  bool homogeneousBC = true;
  //  computeBoundaryN(a_correction, homogeneousBC);
  for (int i=0; i < m_phases; i++){
    mfalias::aliasMF(*m_alias[0], i, a_residual);
    mfalias::aliasMF(*m_alias[1], i, a_correction);
    mfalias::aliasMF(*m_alias[2], i, a_coarseCorrection);
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
    mfalias::aliasMF(*m_alias[0], i, a_coar_resid);
    mfalias::aliasMF(*m_alias[1], i, a_fine_resid);

    Real norm = m_ebops[i]->AMRNorm(*m_alias[0], *m_alias[1], a_ref_rat, a_ord);

    m = Max(m, norm);
  }

  return m;
}
