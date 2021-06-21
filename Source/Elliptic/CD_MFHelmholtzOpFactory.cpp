/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOpFactory.cpp
  @brief  Implementation of CD_MFHelmholtzOpFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <BRMeshRefine.H>

// Our includes
#include <CD_MFHelmholtzOpFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzOpFactory::MFHelmholtzOpFactory(const Real&             a_alpha,
					   const Real&             a_beta,
					   const RealVect&         a_probLo,
					   const AmrLevelGrids&    a_amrLevelGrids,
					   const AmrInterpolators& a_amrInterpolators,
					   const AmrFluxRegisters& a_amrFluxRegisters,
					   const AmrCoarseners&    a_amrCoarseners,
					   const AmrRefRatios&     a_amrRefRatios,
					   const AmrResolutions&   a_amrResolutions,
					   const AmrCellData&      a_amrAcoef,
					   const AmrFluxData&      a_amrBcoef,
					   const AmrIrreData&      a_amrBcoefIrreg,
					   const DomainBCFactory&  a_domainBcFactory,
					   const EBBCFactory&      a_ebbcFactory,
					   const IntVect&          a_ghostPhi,
					   const IntVect&          a_ghostRhs,
					   const RelaxType&        a_relaxtionMethod,
					   const ProblemDomain&    a_bottomDomain,
					   const int&              a_jumpOrder,
					   const int&              a_blockingFactor,
					   const AmrLevelGrids&    a_deeperLevelGrids){


}

MFHelmholtzOpFactory::~MFHelmholtzOpFactory(){

}

void MFHelmholtzOpFactory::defineMultigridLevels(){
  MayDay::Abort("MFHelmholtzOpFactory::defineMultigridLevels -- not implemented");
}

bool MFHelmholtzOpFactory::getCoarserLayout(MFLevelGrid& a_coarEblg, const MFLevelGrid& a_fineEblg, const int a_refRat, const int a_blockingFactor) const {
  MayDay::Abort("MFHelmholtzOpFactory::getCoarserLayout -- not implemented");
}

MFHelmholtzOp* MFHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly) {
  return new MFHelmholtzOp();
}

MFHelmholtzOp* MFHelmholtzOpFactory::AMRnewOp(const ProblemDomain& a_domain) {
  return new MFHelmholtzOp();
}

bool MFHelmholtzOpFactory::isCoarser(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() < B.domainBox().numPts();
}

bool MFHelmholtzOpFactory::isFiner(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() > B.domainBox().numPts();
}

int MFHelmholtzOpFactory::refToFiner(const ProblemDomain& a_indexspace) const {
  MayDay::Abort("MFHelmholtzOpFactory::refToFiner -- not implemented");
}

int MFHelmholtzOpFactory::findAmrLevel(const ProblemDomain& a_domain) const{
  MayDay::Abort("MFHelmholtzOpFactory::findAMRLevel -- not implemented");
}

#include <CD_NamespaceFooter.H>

