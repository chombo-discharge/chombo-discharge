/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzOp.cpp
  @brief  Implementation of CD_EBHelmholtzOp.H
  @author Robert Marskar
  @tood   Replace EBAMRPoissonOp::staticMaxNorm and don't use EBAMRPoissonOp dependencies
*/

// Chombo includes
#include <EBLevelDataOps.H>
#include <EBCellFactory.H>
#include <EBAMRPoissonOp.H>

// Our includes
#include <CD_EBHelmholtzOp.H>
#include <CD_NamespaceHeader.H>


EBHelmholtzOp::EBHelmholtzOp(const EBLevelGrid &                                  a_eblgFine,
			     const EBLevelGrid &                                  a_eblg,
			     const EBLevelGrid &                                  a_eblgCoar,
			     const EBLevelGrid &                                  a_eblgCoarMG,
			     const RefCountedPtr<EBMultigridInterpolator>&        a_quadCFI,
			     const RefCountedPtr<HelmholtzDomainBc>&              a_domainBC,
			     const RefCountedPtr<HelmholtzEbBc>&                  a_ebBC,
			     const Real    &                                      a_dx,
			     const Real    &                                      a_dxCoar,
			     const int&                                           a_refToFine,
			     const int&                                           a_refToCoar,
			     const bool&                                          a_hasFine,
			     const bool&                                          a_hasCoar,
			     const bool&                                          a_hasMGObjects,
			     const bool&                                          a_layoutChanged,
			     const Real&                                          a_alpha,
			     const Real&                                          a_beta,
			     const RefCountedPtr<LevelData<EBCellFAB> >&          a_acoef,
			     const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
			     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
			     const IntVect&                                       a_ghostCellsPhi,
			     const IntVect&                                       a_ghostCellsRHS,
			     const RelaxationMethod&                              a_relaxationMethod) :
  LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false), // Time-independent
  m_eblg(a_eblg),
  m_relaxationMethod(a_relaxationMethod)
{
  
  
}

EBHelmholtzOp::~EBHelmholtzOp(){

}

void EBHelmholtzOp::residual(LevelData<EBCellFAB>& a_residual, const LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs, bool a_homogeneous) {
  MayDay::Warning("EBHelmholtzOp::residual - not implemented");
}

void EBHelmholtzOp::preCond(LevelData<EBCellFAB>& a_corr, const LevelData<EBCellFAB>& a_residual) {
  MayDay::Warning("EBHelmholtzOp::preCond - not implemented");
}

void EBHelmholtzOp::applyOp(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi, bool a_homogeneousPhysBc) {
  MayDay::Warning("EBHelmholtzOp::applyOp - not implemented");
}

void EBHelmholtzOp::create(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) {
  EBCellFactory fact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), a_rhs.nComp(), a_rhs.ghostVect(), fact);
}

void EBHelmholtzOp::assign(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) {
  a_rhs.localCopyTo(a_lhs);
}

Real EBHelmholtzOp::dotProduct(const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) {
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume, a_lhs, a_rhs, EBLEVELDATAOPS_ALLVOFS, domain);
}

void EBHelmholtzOp::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale) {
  EBLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void EBHelmholtzOp::axby(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_x, const LevelData<EBCellFAB>& a_y, const Real a_a, const Real a_b) {
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}

void EBHelmholtzOp::scale(LevelData<EBCellFAB>& a_lhs, const Real& a_scale) {
  EBLevelDataOps::scale(a_lhs, a_scale);
}

Real EBHelmholtzOp::norm(const LevelData<EBCellFAB>& a_rhs, const int a_order) {
  Real maxNorm = EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);

#ifdef CH_MPI
  Real tmp = 1;
  MPI_Allreduce(&maxNorm, &tmp, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  maxNorm = tmp;
#endif

  return maxNorm;
}

void EBHelmholtzOp::setToZero(LevelData<EBCellFAB>& a_lhs) {
  EBLevelDataOps::setToZero(a_lhs);
}

Real EBHelmholtzOp::getSafety() const {
  Real safety;

  switch(m_relaxationMethod){
  case RelaxationMethod::Jacobi:
    safety = 0.5;
    break;
  case RelaxationMethod::GSRB:
    safety = 1.0;
    break;
  case RelaxationMethod::GSRBFast:
    safety = 1.0;
  default:
    MayDay::Abort("EBHelmholtzOp::getSafety - bad relaxation method requested");
  };

  return safety;
}

#include <CD_NamespaceFooter.H>
