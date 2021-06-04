/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */
/*!
  @file   CD_EBHelmholtzOp.H
  @brief  Declaration of Helmholtz multigrid operators. 
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzOp_H
#define CD_EBHelmholtzOp_H

// Chombo includes
#include <LevelTGA.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making variable-coefficient Helmholtz operators. This is meant to be used with LevelTGA so the operator can be used for TGA time stepping. 
*/
class EBHelmholtzOp : public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB> {
public:

  /*!
    @brief Relaxation method for the operators
  */
  enum class RelaxationMethod {
    Jacobi,
    GSRB,
    GSRBFast,
  };

  /*!
    @brief Disallowed default constructor
  */
  EBHelmholtzOp() = delete;

  /*!
    @brief Disallowed copy constructor
  */
  EBHelmholtzOp(const EBHelmholtzOp& a_other) = delete;
#if 0
  /*!
    @brief Full constructor
  */
  EBConductivityOp(const EBLevelGrid &                                  a_eblgFine,
                   const EBLevelGrid &                                  a_eblg,
                   const EBLevelGrid &                                  a_eblgCoar,
                   const EBLevelGrid &                                  a_eblgCoarMG,
                   const RefCountedPtr<EBQuadCFInterp>&                 a_quadCFI,
                   const RefCountedPtr<ConductivityBaseDomainBC>&       a_domainBC,
                   const RefCountedPtr<ConductivityBaseEBBC>&           a_ebBC,
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
                   const EBHelmholtzOpFactory::RelaxationMethod&                                           a_relaxType);  
  EBHelmholtzOp(

  }
#endif

  
protected:

};

#include <CD_NamespaceFooter.H>

#endif
