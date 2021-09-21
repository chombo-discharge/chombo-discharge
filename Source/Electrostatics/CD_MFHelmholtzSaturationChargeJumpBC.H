/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzSaturationChargeJumpBC.H
  @brief  Declaration of a jump BC which implements the saturation charge boundary condition. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzSaturationChargeJumpBC_H
#define CD_MFHelmholtzSaturationChargeJumpBC_H

// Our includes
#include <CD_MFHelmholtzJumpBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for computing matching the saturation charge boundary condition in an MFHelmholtzOp context
*/
class MFHelmholtzSaturationChargeJumpBC : public MFHelmholtzJumpBC {
public:

  /*!
    @brief Default constructor
  */
  MFHelmholtzSaturationChargeJumpBC() = delete;

  /*!
    @brief Full constructor. 
    @param[in] a_dataLocation Data centering
    @param[in] a_mflg         Grids
    @param[in] a_Bcoef        Helmholtz B-coefficient
    @param[in] a_dx           Grid resolution
    @param[in] a_order        Stencil approximation order
    @param[in] a_weight       Weighting factor for least squares
    @param[in] a_radius       Stencil radius
    @param[in] a_ghostCF      Number of grid cells that were filled across the CF. 
  */
  MFHelmholtzSaturationChargeJumpBC(const phase::which_phase a_phase,
				    const Location::Cell     a_dataLocation,
				    const MFLevelGrid&       a_mflg,
				    const BcoefPtr&          a_Bcoef,
				    const Real               a_dx,
				    const int                a_order,
				    const int                a_weight,
				    const int                a_radius,
				    const int                a_ghostCF);
  
  /*!
    @brief Destructor (does nothing)
  */
  ~MFHelmholtzSaturationChargeJumpBC();

  /*!
    @brief Match the BC. 
    @details Unlike JumpBC which uses a_jump as the natural parameter, a_jump 
    @param[in] a_phi               Cell-centered data
    @param[in] a_jump              Jump data, i.e. right-hand side of dphi/dn + dphi/m = jump
    @param[in] a_homogeneousPhysBC Homogeneous BC or not
  */
  virtual void matchBC(BaseIVFAB<Real>& a_jump,
		       const MFCellFAB& a_phi,
		       const bool       a_homogeneousPhysBC,
		       const DataIndex& a_dit) const override;

protected:

  /*!
    @brief Phase where the saturation charge BC is applied
  */
  phase::which_phase m_phase;
};

#include <CD_NamespaceFooter.H>

#endif
