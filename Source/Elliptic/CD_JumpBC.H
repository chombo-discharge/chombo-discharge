/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_JumpBC.H
  @brief  Class for computing "jump interface" boundary conditions for multifluid code. 
  @author Robert Marskar
*/

#ifndef CD_JumpBC_H
#define CD_JumpBC_H

// Chombo includes
#include <MFCellFAB.H>
#include <EBISBox.H>
#include <LayoutData.H>
#include <Stencils.H>

// Our includes
#include <CD_MFLevelGrid.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief  Class for computing "jump interface" boundary conditions for multifluid code. 
*/

class JumpBC {
public:

  /*!
    @brief Disallowed weak construction
  */
  JumpBC() = delete;

  /*!
    @brief Defining constructor
    @param[in] a_mflg   Grids
    @param[in] a_dx     Grid resolution
    @param[in] a_order  Stencil approximation order
    @param[in] a_weight Weighting factor for least squares
    @param[in] a_radius Stencil radius
  */
  JumpBC(const MFLevelGrid& a_mflg, const LevelData<MFBaseIVFAB>& a_Bcoef, const Real a_dx, const int a_order, const int a_weight, const int a_radius);

  /*!
    @breif Return stencils
  */
  const LayoutData<MFInterfaceFAB<VoFStencil> >& getStencils() const;


protected:
  
};

#include <CD_NamespaceFooter.H>

#endif
