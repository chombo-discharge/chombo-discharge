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
#include <CD_VofUtils.H>
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
  JumpBC(const MFLevelGrid& a_mflg, const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_Bcoef, const Real a_dx, const int a_order, const int a_weight, const int a_radius);

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
  void define(const MFLevelGrid& a_mflg, const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_Bcoef, const Real a_dx, const int a_order, const int a_weight, const int a_radius);

  /*!
    @brief Match the BC
    @param[in] a_phi               Cell-centered data
    @param[in] a_jump              Jump data, i.e. right-hand side of dphi/dn + dphi/m = jump
    @param[in] a_homogeneousPhysBC Homogeneous BC or not
  */
  void matchBC(const LevelData<MFCellFAB>& a_phi,
	       const LevelData<BaseIVFAB<Real> >& a_jump,
	       const bool a_homogeneousPhysBC) const;

  /*!
    @brief Get phi on the boundary. 
    @param[in] a_phase Phase
    @param[in] a_dit   Grid index
    @return Returns phi on the boundary. On multivalued cells, phi always lies on VolIndex(iv, 0). 
  */
  virtual const BaseIVFAB<Real>& getBndryPhi(const int a_phase, const DataIndex& a_dit) const;
  
protected:

  /*!
    @brief Grids
  */
  MFLevelGrid m_mflg;

  /*!
    @brief B-coefficient
  */
  RefCountedPtr<LevelData<MFBaseIVFAB> > m_Bcoef;

  /*!
    @brief Grid resolution
  */
  Real m_dx;

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Stencil order
  */
  int m_weight;

  /*!
    @brief Stencil order
  */
  int m_radius;

  /*!
    @brief Value of phi on the boundary. 
  */
  LayoutData<MFInterfaceFAB<Real> > m_boundaryPhi;

  /*!
    @brief Average stencil
  */
  LayoutData<MFInterfaceFAB<VoFStencil> > m_avgStencils;

  /*!
    @brief Average weights
  */
  LayoutData<MFInterfaceFAB<Real> > m_avgWeights;

  /*!
    @brief Average b-coefficient
  */
  LayoutData<MFInterfaceFAB<Real> > m_avgBco;

  /*!
    @brief Get stencil with specified order. 
    @param[out] a_stencil      Weight and stencil
    @param[out] a_vof          Grid cell
    @param[out] a_neighborhood Neighboorhood specifier
    @param[in]  a_dit          Data index
    @param[in]  a_order        Attempted order for stencil. This is also the stencil radius. 
    @param[in]  a_phase        Phase. 
    @return If this returns true, a_stencil contains a valid stencil and weight. 
    @note Stencils are computed using weighted least squares with squared inverse distance weighting. 
  */
  virtual bool getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
						  const VolIndex&              a_vof,
						  const VofUtils::Neighborhood a_neighborhood,
						  const DataIndex&             a_dit,
						  const int                    a_order,
						  const int                    a_phase) const;
};

#include <CD_NamespaceFooter.H>

#endif
