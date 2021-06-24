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

  using BcoefPtr = RefCountedPtr<LevelData<MFBaseIVFAB> >;

  /*!
    @brief Disallowed weak construction
  */
  JumpBC() = delete;

  /*!
    @brief Defining constructor
    @param[in] a_mflg    Grids
    @param[in] a_Bcoef   Helmholtz B-coefficient
    @param[in] a_dx      Grid resolution
    @param[in] a_order   Stencil approximation order
    @param[in] a_weight  Weighting factor for least squares
    @param[in] a_radius  Stencil radius
    @param[in] a_ghostCF Number of grid cells that were filled across the CF. 
  */
  JumpBC(const MFLevelGrid& a_mflg,
	 const BcoefPtr&    a_Bcoef,
	 const Real         a_dx,
	 const int          a_order,
	 const int          a_weight,
	 const int          a_radius,
	 const int          a_ghostCF);

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
    @brief Set special flag in case we need to drop stencil radius when the domain is small. 
  */
  void setMG(const bool a_isMGLevel);

  /*!
    @brief Match the BC
    @param[in] a_phi               Cell-centered data
    @param[in] a_jump              Jump data, i.e. right-hand side of dphi/dn + dphi/m = jump
    @param[in] a_homogeneousPhysBC Homogeneous BC or not
  */
  void matchBC(const LevelData<MFCellFAB>&        a_phi,
	       const LevelData<BaseIVFAB<Real> >& a_jump,
	       const bool                         a_homogeneousPhysBC) const;

  /*!
    @brief Get phi on the boundary. 
    @param[in] a_phase Phase
    @param[in] a_dit   Grid index
    @return Returns phi on the boundary. On multivalued cells, phi always lies on VolIndex(iv, 0). 
  */
  virtual const BaseIVFAB<Real>& getBndryPhi(const int a_phase, const DataIndex& a_dit) const;

  /*!
    @brief Return vof iterator which iterates over single-phase cells on a specific phase and box. 
  */
  virtual VoFIterator& getSinglePhaseVofs(const int a_phase, const DataIndex& a_dit) const;

  /*!
    @brief Return vof iterator which iterates over multi-phase cells on a specified phase and box. 
  */
  virtual VoFIterator& getMultiPhaseVofs(const int a_phase, const DataIndex& a_dit) const;

  
  /*!
    @brief Set everything to zero. This is a debugging function.
  */
  void resetBC() const;
  
protected:

  /*!
    @brief Component where stencils are stored. 
  */
  static constexpr int m_comp = 0;

    /*!
    @brief Number of components in stencil
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Grids
  */
  MFLevelGrid m_mflg;

  /*!
    @brief B-coefficient
  */
  RefCountedPtr<LevelData<MFBaseIVFAB> > m_Bcoef;

  /*!
    @brief Is MG level or not
  */
  bool m_isMGLevel;

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
    @brief Number of grid cells that were filled over the CF 
  */
  int m_ghostCF;

  /*!
    @brief Number of phases
  */
  int m_numPhases;

  /*!
    @brief Multiphase or not
  */
  bool m_multiPhase;

  /*!
    @brief Value of phi on the boundary. 
  */
  mutable LayoutData<MFInterfaceFAB<Real> > m_boundaryPhi;
  
  /*!
    @brief Iterators for single-phase cells. The integer index is the phase. 
  */
  mutable std::map<int, std::shared_ptr<LayoutData<VoFIterator> > > m_singlePhaseVofs;

  /*!
    @brief Iterators for multi-phase cells. The integer index is the phase. 
  */
  mutable std::map<int, std::shared_ptr<LayoutData<VoFIterator> > > m_multiPhaseVofs;

  /*!
    @brief IntVectSet over which we run the "averaging"
  */
  LayoutData<IntVectSet> m_ivs;

  /*!
    @brief Average stencil
  */
  LayoutData<MFInterfaceFAB<VoFStencil> > m_avgStencils;

  /*!
    @brief Denom factor
  */
  LayoutData<MFInterfaceFAB<Real> > m_denom;

  /*!
    @brief Average weights
  */
  LayoutData<MFInterfaceFAB<Real> > m_avgWeights;

  /*!
    @brief Define function. Builds stencils. 
  */
  void defineStencils();

  /*!
    @brief Define function. Builds stencils. 
  */
  void defineIterators();

  /*!
    @brief Get stencil with specified order. 
    @param[out] a_stencil      Weight and stencil
    @param[out] a_vof          Grid cell
    @param[out] a_ebisbox      EBISBox
    @param[out] a_neighborhood Neighboorhood specifier
    @param[in]  a_order        Attempted order for stencil. This is also the stencil radius. 
    @return If this returns true, a_stencil contains a valid stencil and weight. 
    @note Stencils are computed using weighted least squares with squared inverse distance weighting. 
  */
  virtual bool getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
						  const VolIndex&              a_vof,
						  const EBISBox&               a_ebisbox,
						  const VofUtils::Neighborhood a_neighborhood,
						  const int                    a_order) const;

  /*!
    @brief Check if stencil is valid. 
    @details This routine makes sure that the stencil a_stencil does not reach over the CF interface and into ghost cells that do not have valid data.
    @param[in] a_stencil Stencil
    @param[in] a_dit     Data index
  */
  inline
  bool isStencilValidCF(const VoFStencil& a_stencil, const DataIndex& a_dit) const;

  /*!
    @brief Match the BC
    @param[in] a_phi               Cell-centered data
    @param[in] a_jump              Jump data, i.e. right-hand side of dphi/dn + dphi/m = jump
    @param[in] a_homogeneousPhysBC Homogeneous BC or not
  */
  inline
  void matchBC(const MFCellFAB&       a_phi,
	       const BaseIVFAB<Real>& a_jump,
	       const bool             a_homogeneousPhysBC,
	       const DataIndex&       a_dit) const;

  /*!
    @brief Apply a stencil and return the result
    @param[in] a_stencil Stencil to apply
    @param[in] a_phi     Data
    @returns Returns S(phi), S being the stencil. 
  */
  inline
  Real applyStencil(const VoFStencil& a_stencil, const EBCellFAB& a_phi) const;

};

#include <CD_NamespaceFooter.H>

#include <CD_JumpBCImplem.H>

#endif
