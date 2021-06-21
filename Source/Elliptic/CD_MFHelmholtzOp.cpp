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

MFHelmholtzOp::MFHelmholtzOp(){
  MayDay::Abort("MFHelmholtzOp - weak construction is not allowed");
}

MFHelmholtzOp::~MFHelmholtzOp(){

}

void MFHelmholtzOp::define(){

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

Real MFHelmholtzOp::dotProduct(const LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) {
  MayDay::Abort("MFHelmholtzOp::dotProduct - not implemented");

  Real accum = 0.0;
  Real volum = 0.0;

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    const MFCellFAB& lhs = a_lhs[dit()];

    Real phaseVolume;
    for (int i = 0; i < lhs.numPhases(); i++){
      const EBCellFAB& data1 = a_lhs[dit()].getPhase(i);
      const EBCellFAB& data2 = a_rhs[dit()].getPhase(i);
      const Box box = a_lhs.disjointBoxLayout()[dit()];

      accum += EBLevelDataOps::sumKappaDotProduct(phaseVolume, data1, data2, box, EBLEVELDATAOPS_ALLVOFS, m_mflg.getDomain());
      volum += phaseVolume;
    }
  }

#ifdef CH_MPI
  Real recv;
  MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm); accum = recv;
  MPI_Allreduce(&volum, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm); volum = recv;
#endif

  Real ret = 0.0;
  if(volum > 0.0){
    ret = accum/volum;
  }

  return ret;
}

void MFHelmholtzOp::create(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs){
  m_ops.create(a_lhs, a_rhs);
}

void MFHelmholtzOp::createCoarser(LevelData<MFCellFAB>& a_coarse, const LevelData<MFCellFAB>& a_fine, bool a_ghosted) {
  MayDay::Warning("MFHelmholtzOp::createCoarser - not implemented");
}

void MFHelmholtzOp::createCoarsened(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const int& a_refRat) {
  MayDay::Warning("MFHelmholtzOp::createCoarsened - not implemented");
}

void MFHelmholtzOp::preCond(LevelData<MFCellFAB>& a_corr, const LevelData<MFCellFAB>& a_residual) {
  this->relax(a_corr, a_residual, 40);
}

void MFHelmholtzOp::applyOp(LevelData<MFCellFAB>& a_Lphi, const LevelData<MFCellFAB>& a_phi, bool a_homogeneousPhysBC) {
  this->applyOp(a_Lphi, a_phi, nullptr, a_homogeneousPhysBC, true);
}

void MFHelmholtzOp::applyOp(LevelData<MFCellFAB>&             a_Lphi,
			    const LevelData<MFCellFAB>&       a_phi,
			    const LevelData<MFCellFAB>* const a_phiCoar,
			    const bool                        a_homogeneousPhysBC,
			    const bool                        a_homogeneousCFBC){
  
  // We MUST have updated ghost cells before applying the matching conditions. 
  this->interpolateCF(a_phi, a_phiCoar, a_homogeneousCFBC);
  this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  // Apply the operator on each level. 
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB>  Lphi;
    LevelData<EBCellFAB>  phi;
    LevelData<EBCellFAB>* phiCoar = nullptr;

    // applyOp usually interpolates ghost cells, but we did that above in interpolateCF(...). 
    op.second->turnOffBCs();
    op.second->applyOp(Lphi, phi, phiCoar, a_homogeneousPhysBC, a_homogeneousCFBC);
    op.second->turnOnBCs();
  }
}

void MFHelmholtzOp::interpolateCF(const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>* a_phiCoar, const bool a_homogeneousCF){
  if(m_hasCoar){
    for (auto& op : m_helmOps){
      LevelData<EBCellFAB> phi;
      LevelData<EBCellFAB> phiCoar;
      
      if(a_homogeneousCF){
	MultifluidAlias::aliasMF(phi, op.first, (LevelData<MFCellFAB>&) a_phi);

	op.second->homogeneousCFInterp(phi);
      }
      else{
	MultifluidAlias::aliasMF(phi,     op.first,  (LevelData<MFCellFAB>&) a_phi);
	MultifluidAlias::aliasMF(phiCoar, op.first, *a_phiCoar);

	op.second->inhomogeneousCFInterp(phi, phiCoar);
      }
    }
  }
}

void MFHelmholtzOp::residual(LevelData<MFCellFAB>& a_residual, const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_rhs, const bool a_homogeneousPhysBC){
  // Compute a_residual = rhs - L(phi) 
  this->applyOp(a_residual, a_phi, a_homogeneousPhysBC);
  this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

void MFHelmholtzOp::axby(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& x, const LevelData<MFCellFAB>& y, const Real a, const Real b) {
  m_ops.axby(a_lhs, x, y, a, b);
}

void MFHelmholtzOp::updateJumpBC(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneousPhysBC){
  MayDay::Warning("MFHelmholtzOp::updateJumpBC -- not implemented (yet)");
}

void MFHelmholtzOp::relax(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, int a_iterations) {
  MayDay::Warning("MFHelmholtzOp::relax - not implemented");
}

void MFHelmholtzOp::restrictResidual(LevelData<MFCellFAB>& a_resCoar, LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_rhs) {
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> resCoar;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> rhs;

    MultifluidAlias::aliasMF(resCoar, op.first, a_resCoar);
    MultifluidAlias::aliasMF(phi,     op.first, a_phi);
    MultifluidAlias::aliasMF(rhs,     op.first, a_rhs);

    op.second->restrictResidual(resCoar, phi, rhs);
  }
}

void MFHelmholtzOp::prolongIncrement(LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_correctCoarse) {
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> correctCoarse;

    MultifluidAlias::aliasMF(phi,           op.first, a_phi);
    MultifluidAlias::aliasMF(correctCoarse, op.first, a_correctCoarse);

    op.second->prolongIncrement(phi, correctCoarse);
  }
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

void MFHelmholtzOp::AMRResidual(LevelData<MFCellFAB>&              a_residual,
				const LevelData<MFCellFAB>&        a_phiFine,
				const LevelData<MFCellFAB>&        a_phi,
				const LevelData<MFCellFAB>&        a_phiCoar,
				const LevelData<MFCellFAB>&        a_rhs,
				bool                               a_homogeneousPhysBC,
				AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {

  // Make residual = a_rhs - L(phi)
  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar, a_homogeneousPhysBC, a_finerOp);
  this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

void MFHelmholtzOp::AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
				  const LevelData<MFCellFAB>& a_phi,
				  const LevelData<MFCellFAB>& a_phiCoar,
				  const LevelData<MFCellFAB>& a_rhs,
				  bool                        a_homogeneousPhysBC) {
  
  // Make residual = a_rhs - L(phi)  
  this->AMROperatorNF(a_residual, a_phi, a_phiCoar, a_homogeneousPhysBC);
  this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

void MFHelmholtzOp::AMRResidualNC(LevelData<MFCellFAB>&              a_residual,
				  const LevelData<MFCellFAB>&        a_phiFine,
				  const LevelData<MFCellFAB>&        a_phi,
				  const LevelData<MFCellFAB>&        a_rhs,
				  bool                               a_homogeneousPhysBC,
				  AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {

  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousPhysBC, a_finerOp);
  this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

void MFHelmholtzOp::AMROperatorNF(LevelData<MFCellFAB>&       a_Lphi,
				  const LevelData<MFCellFAB>& a_phi,
				  const LevelData<MFCellFAB>& a_phiCoar,
				  bool                        a_homogeneousPhysBC) {

  this->interpolateCF(a_phi, &a_phiCoar, false);
  this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;

    MultifluidAlias::aliasMF(Lphi,    op.first, a_Lphi);
    MultifluidAlias::aliasMF(phi,     op.first, a_phi);
    MultifluidAlias::aliasMF(phiCoar, op.first, a_phiCoar);

    op.second->turnOffBCs(); // Don't need to interpolate ghost cells again. 
    op.second->AMROperatorNF(Lphi, phi, phiCoar, a_homogeneousPhysBC);
    op.second->turnOnBCs(); 
  }
}

void MFHelmholtzOp::AMROperatorNC(LevelData<MFCellFAB>&              a_Lphi,
				  const LevelData<MFCellFAB>&        a_phiFine,
				  const LevelData<MFCellFAB>&        a_phi,
				  bool                               a_homogeneousPhysBC,
				  AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {
  LevelData<MFCellFAB> phiCoar;
  this->AMROperator(a_Lphi, a_phiFine, a_phi, phiCoar, a_homogeneousPhysBC, a_finerOp);
}

void MFHelmholtzOp::AMROperator(LevelData<MFCellFAB>&              a_Lphi,
				const LevelData<MFCellFAB>&        a_phiFine,
				const LevelData<MFCellFAB>&        a_phi,
				const LevelData<MFCellFAB>&        a_phiCoar,
				const bool                         a_homogeneousPhysBC,
				AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {

  // Update ghost cells and jump conditions first.
  this->interpolateCF(a_phi, &a_phiCoar, false);
  this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phiFine;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;

    MultifluidAlias::aliasMF(Lphi,    op.first, a_Lphi);
    MultifluidAlias::aliasMF(phiFine, op.first, a_phiFine);
    MultifluidAlias::aliasMF(phi,     op.first, a_phi);
    MultifluidAlias::aliasMF(phiCoar, op.first, a_phiCoar);

    MFHelmholtzOp* finerOp = (MFHelmholtzOp*) (a_finerOp);

    op.second->turnOffBCs(); // Don't need to interpolate ghost cells again. 
    op.second->AMROperator(Lphi, phiFine, phi, phiCoar, a_homogeneousPhysBC, (finerOp->m_helmOps).at(op.first));
    op.second->turnOnBCs();
  }
}

#include <CD_NamespaceFooter.H>
