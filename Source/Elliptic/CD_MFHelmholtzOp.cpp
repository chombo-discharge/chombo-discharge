/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOp.cpp
  @brief  Implementation of CD_MFHelmholtzOp.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzOp.H>
#include <CD_MultifluidAlias.H>
#include <CD_NamespaceHeader.H>

constexpr int MFHelmholtzOp::m_numAlias;

MFHelmholtzOp::MFHelmholtzOp(){
  MayDay::Abort("MFHelmholtzOp - weak construction is not allowed");
}

MFHelmholtzOp::~MFHelmholtzOp(){

}

void MFHelmholtzOp::define(){

  // Define aliases. 
  m_alias.resize(m_numAlias);
  for (auto& a : m_alias) {
    a = new LevelData<EBCellFAB>();
  }
}

int MFHelmholtzOp::refToCoarser() {
  return m_refToCoar;
}

unsigned int MFHelmholtzOp::orderOfAccuracy(void) const {
  return 99;
}

void MFHelmholtzOp::enforceCFConsistency(LevelData<MFCellFAB>& a_coarCorr, const LevelData<MFCellFAB>& a_fineCorr) {
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> coarCorr;
    LevelData<EBCellFAB> fineCorr;
    
    MultifluidAlias::aliasMF(coarCorr, op.first, a_coarCorr);
    MultifluidAlias::aliasMF(fineCorr, op.first, a_fineCorr);

    op.second->enforceCFConsistency(coarCorr, fineCorr);
  }
}

void MFHelmholtzOp::incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, Real a_scale){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
}

void MFHelmholtzOp::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale) {
  m_ops.scale(a_lhs, a_scale);
}

void MFHelmholtzOp::setToZero(LevelData<MFCellFAB>& a_lhs) {
  m_ops.setToZero(a_lhs);
}

void MFHelmholtzOp::assign(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) {
  m_ops.assign(a_lhs, a_rhs);
}

Real MFHelmholtzOp::norm(const LevelData<MFCellFAB>& a_lhs, int a_order){
  Real norm = 0.0;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> lhs;

    MultifluidAlias::aliasMF(lhs, op.first, a_lhs);

    const Real opNorm = op.second->norm(lhs, a_order);

    norm = std::max(norm, opNorm);
  }

  return norm;
}

void MFHelmholtzOp::create(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs){
  m_ops.create(a_lhs, a_rhs);
}

void MFHelmholtzOp::preCond(LevelData<MFCellFAB>& a_corr, const LevelData<MFCellFAB>& a_residual) {
  this->relax(a_corr, a_residual, 40);
}

void MFHelmholtzOp::applyOp(LevelData<MFCellFAB>& a_Lphi, const LevelData<MFCellFAB>& a_phi, bool a_homogeneousPhysBc) {
  MayDay::Abort("MFHelmholtzOp::applyOp - not implemented");
}

void MFHelmholtzOp::residual(LevelData<MFCellFAB>& a_residual, const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_rhs, const bool a_homogeneousPhysBC){
  // Compute a_residual = rhs - L(phi) 
  this->applyOp(a_residual, a_phi, a_rhs, a_homogeneousPhysBC);
  this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

void MFHelmholtzOp::axby(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& x, const LevelData<MFCellFAB>& y, const Real a, const Real b) {
  m_ops.axby(a_lhs, x, y, a, b);
}

void MFHelmholtzOp::relax(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, int a_iterations) {
  MayDay::Abort("MFHelmholtzOp::relax - not implemented");
}

void MFHelmholtzOp::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
				      const LevelData<MFCellFAB>& a_correction,
				      const LevelData<MFCellFAB>& a_coarseCorrection){

  // Need to update BC first!
  //this->updateJumpBC(a_correction, true);

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> residual;
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;

    MultifluidAlias::aliasMF(residual,         op.first, a_residual);
    MultifluidAlias::aliasMF(correction,       op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    op.second->AMRUpdateResidual(residual, correction, coarseCorrection);
  }
}

void MFHelmholtzOp::AMRRestrict(LevelData<MFCellFAB>&       a_residualCoarse,
				const LevelData<MFCellFAB>& a_residual,
				const LevelData<MFCellFAB>& a_correction,
				const LevelData<MFCellFAB>& a_coarseCorrection,
				bool                        a_skip_res) {
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> residualCoarse;
    LevelData<EBCellFAB> residual;
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;
    
    MultifluidAlias::aliasMF(residualCoarse,   op.first, a_residualCoarse);
    MultifluidAlias::aliasMF(residual,         op.first, a_residual);
    MultifluidAlias::aliasMF(correction,       op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    op.second->AMRRestrict(residualCoarse, residual, correction, coarseCorrection, a_skip_res);
  }
}

void MFHelmholtzOp::AMRProlong(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_coarseCorrection) {
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;
    
    MultifluidAlias::aliasMF(correction,       op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    op.second->AMRProlong(correction, coarseCorrection);
  }
}

#include <CD_NamespaceFooter.H>
