/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzSaturationChargeJumpBCFactory.H
  @brief  Declaration of a factory class for making MFHelmholtzSaturationChargeJumpBC objects for multifluid Helmholtz. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzSaturationChargeJumpBCFactory_H
#define CD_MFHelmholtzSaturationChargeJumpBCFactory_H

// Our includes
#include <CD_MFHelmholtzJumpBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making saturation charge jump boundary condition objects in multifluid Helmholtz context. 
*/
class MFHelmholtzSaturationChargeJumpBCFactory : public MFHelmholtzJumpBCFactory
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
    @brief Full constructor
    @param[in] a_phase Which phase the saturation charge BC is applied on.
  */
  MFHelmholtzSaturationChargeJumpBCFactory(const phase::which_phase a_phase);

  /*!
    @brief Default constructor
  */
  virtual ~MFHelmholtzSaturationChargeJumpBCFactory();

  /*!
    @brief Factory method for creating MFHelmholtzSaturationChargeJumpBC objects. 
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
    @return Returns an MFHelmholtzJumpBC for encapsulating boundary conditions on multiphase cells. 
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
         const IntVect        a_ghostPhi) override;

protected:
  /*!
    @brief Which phase the saturation BC is applied on
  */
  phase::which_phase m_phase;
};

#include <CD_NamespaceFooter.H>

#endif
