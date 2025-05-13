/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzJumpBCFactory.H
  @brief  Declaration of a factory class for making MFHelmholtzJumpBC objects for use in a multifluid MFHelmholtzOp context. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzJumpBCFactory_H
#define CD_MFHelmholtzJumpBCFactory_H

// Our includes
#include <CD_MFHelmholtzJumpBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making jump boundary condition objects in a multifluid MFHelmholtzOp context. 
*/
class MFHelmholtzJumpBCFactory
{
public:
  /*!
    @brief Alias to cut down on typing
  */
  using BcoefPtr = RefCountedPtr<LevelData<MFBaseIVFAB>>;

  /*!
    @brief Alias to cut down on typing
  */
  using AmrMask = RefCountedPtr<LevelData<BaseFab<bool>>>;

  /*!
    @brief Weak constructor. 
  */
  MFHelmholtzJumpBCFactory();

  /*!
    @brief Default constructor
  */
  virtual ~MFHelmholtzJumpBCFactory();

  /*!
    @brief Factory method for creating MFHelmholtzJumpBC objects. 
    @param[in] a_dataLocation Data centering
    @param[in] a_mflg         Grids
    @param[in] a_Bcoef        Helmholtz B-coefficient
    @param[in] a_validCells   Valid grid cells (might be null pointer)
    @param[in] a_dx           Grid resolution
    @param[in] a_order        Stencil approximation order
    @param[in] a_weight       Weighting factor for least squares
    @param[in] a_radius       Stencil radius
    @param[in] a_ghostCF      Number of grid cells that were filled across the CF. 
    @param[in] a_ghostPhi     Number of ghost cells in phi
  */
  virtual RefCountedPtr<MFHelmholtzJumpBC>
  create(const Location::Cell a_dataLocation,
         const MFLevelGrid&   a_mflg,
         const BcoefPtr&      a_Bcoef,
         const AmrMask&       a_validCells,
         const Real           a_dx,
         const int            a_order,
         const int            a_weight,
         const int            a_radius,
         const int            a_ghostCF,
         const IntVect        a_ghostPhi);

  /*!
    @brief Drop BC order if domain size is equal or below this.
  */
  virtual void
  setDomainDropOrder(const int a_domainSize);

protected:
  /*!
    @brief Special flag for dropping stencil order when domains become coarser than this.
  */
  int m_domainDropOrder;
};

#include <CD_NamespaceFooter.H>

#endif
