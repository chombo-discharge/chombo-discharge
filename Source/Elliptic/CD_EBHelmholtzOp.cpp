/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzOp.cpp
  @brief  Implementation of CD_EBHelmholtzOp.H
  @author Robert Marskar
*/

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
			     const RelaxationMethod&                              a_relaxationMethod){

}

EBHelmholtzOp::~EBHelmholtzOp(){

}

#include <CD_NamespaceFooter.H>
