/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_JumpBC.H
  @brief  Declaration of class for computing "jump interface" boundary conditions for multifluid Helmholtz code. 
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
    @brief Disallowed copy constructor
  */
  JumpBC(const JumpBC& a_other) = delete;

  /*!
    @brief Destructor. Does nothing.
  */
  ~JumpBC();

  /*!
    @brief Disallowed assignment constructor
  */
  JumpBC& operator=(const JumpBC& a_other) = delete;

  /*!
    @brief Define function. Builds stencils. 
  */
  void define(const MFLevelGrid& a_mflg, const LevelData<MFBaseIVFAB>& a_Bcoef, const Real a_dx, const int a_order, const int a_weight, const int a_radius);
  
  /*!
    @brief Return average stencils over EB faces. 
  */
  virtual const BaseIVFAB<VoFStencil>& getAvgStencils(const int a_phase, const DataIndex& a_dit) const;

  /*!
    @brief Get inhomogeneous contribution (part that is proportional to sigma). 
    @param[in] a_phase Phase
    @param[in] a_dit   Grid index
    @return Returns inhomgoeneous contribution to phase a_phase
  */
  virtual const BaseIVFAB<Real>& getInhomogeneousContribution(const int a_phase, const DataIndex& a_dit) const;

  /*!
    @brief Get homogeneous contribution to phase a_phase (part that only affects interior cells, including cells on the other phase). 
    @param[in] a_phase Phase
    @param[in] a_dit   Grid index
    @return Returns homogeneous contribution to phase a_phase
  */
  virtual const BaseIVFAB<Real>& getHomogeneousContribution(const int a_phase, const DataIndex& a_dit) const;
  
protected:

  MFLevelGrid m_mflg;

  /*!
    @brief Computed inhomogeneous contribution
  */
  LayoutData<MFInterfaceFAB<Real> > m_inhomogeneousContrib;

  /*!
    @brief Computed homogeneous contribution
  */
  LayoutData<MFInterfaceFAB<Real> > m_homogeneousContrib;

  /*!
    @brief Average stencil
  */
  LayoutData<MFInterfaceFAB<VoFStencil> > m_avgStencils;
};

#include <CD_NamespaceFooter.H>

#endif
