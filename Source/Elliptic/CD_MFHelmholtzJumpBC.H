/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzJumpBC.H
  @brief  Declaration of class for computing "jump interface" boundary conditions for multifluid Helmholtz code. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzJumpBC_H
#define CD_MFHelmholtzJumpBC_H

// Chombo includes
#include <MFCellFAB.H>
#include <EBISBox.H>
#include <LayoutData.H>
#include <Stencils.H>
#include <AggStencil.H>

// Our includes
#include <CD_Location.H>
#include <CD_MFLevelGrid.H>
#include <CD_VofUtils.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief  Class for computing "jump interface" boundary conditions for multifluid code. 
  @details This class creates stencils for approximating the normal derivative on each side. It uses those stencil to compute a value on the boundary which fulfills
  the discrete version of dphi/dn1 + dphi/dn2 = sigma. This value is used by MFHelmholtzEBBC to compute the fluxes into the cut-cells on each side. 

  When multi-valued cells are involved we use an "average" formulation such that we obtain the surface potential which is equal on all multi-cut surfaces. For this
  reason we store the bare stencils AND the averaged stencils (which are used in the matching). When stencil coefficients are updated then we use the raw stencils
  for recomputing the average stencils. This is done by calling setBco.
*/
class MFHelmholtzJumpBC
{
public:
  /*!
    @brief Alias to cut down on typing
  */
  using BcoefPtr = RefCountedPtr<LevelData<MFBaseIVFAB>>;

  /*!
    @brief Disallowed weak construction
  */
  MFHelmholtzJumpBC() = delete;

  /*!
    @brief Defining constructor.
    @param[in] a_dataLocation Data centering
    @param[in] a_mflg         Grids
    @param[in] a_Bcoef        Helmholtz B-coefficient
    @param[in] a_dx           Grid resolution
    @param[in] a_order        Stencil approximation order
    @param[in] a_weight       Weighting factor for least squares
    @param[in] a_radius       Stencil radius
    @param[in] a_ghostCF      Number of grid cells that were filled across the CF. 
  */
  MFHelmholtzJumpBC(const Location::Cell a_dataLocation,
                    const MFLevelGrid&   a_mflg,
                    const BcoefPtr&      a_Bcoef,
                    const Real           a_dx,
                    const int            a_order,
                    const int            a_weight,
                    const int            a_radius,
                    const int            a_ghostCF);

  /*!
    @brief Disallowed copy constructor
  */
  MFHelmholtzJumpBC(const MFHelmholtzJumpBC& a_other) = delete;

  /*!
    @brief Destructor. Does nothing.
  */
  virtual ~MFHelmholtzJumpBC();

  /*!
    @brief Disallowed assignment constructor
  */
  MFHelmholtzJumpBC&
  operator=(const MFHelmholtzJumpBC& a_other) = delete;

  /*!
    @brief Set B-coefficient. This updates the stencils
    @param[in] a_Bcoef Helmholtz B-coefficient
  */
  void
  setBco(const RefCountedPtr<LevelData<MFBaseIVFAB>>& a_Bcoef);

  /*!
    @brief Return stencil order
  */
  int
  getOrder() const;

  /*!
    @brief Return stencil weight (for least squares)
  */
  int
  getWeight() const;

  /*!
    @brief Return stencil radius (for least squares)
  */
  int
  getRadius() const;

  /*!
    @brief Match the BC
    @param[inout] a_jump              Jump data, i.e. right-hand side of dphi/dn + dphi/m = jump
    @param[in]    a_phi               Cell-centered data
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not
  */
  virtual void
  matchBC(LevelData<BaseIVFAB<Real>>& a_jump, const LevelData<MFCellFAB>& a_phi, const bool a_homogeneousPhysBC) const;

  /*!
    @brief Match the BC
    @param[inout] a_jump              Jump data, i.e. right-hand side of dphi/dn + dphi/m = jump
    @param[in]    a_phi               Cell-centered data
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not
    @param[in]    a_dit               Grid index
  */
  virtual void
  matchBC(BaseIVFAB<Real>& a_jump,
          const MFCellFAB& a_phi,
          const bool       a_homogeneousPhysBC,
          const DataIndex& a_dit) const;

  /*!
    @brief Get phi on the boundary. 
    @param[in] a_phase Phase
    @param[in] a_dit   Grid index
    @return Returns phi on the boundary. On multivalued cells, phi always lies on VolIndex(iv, 0). 
  */
  virtual const BaseIVFAB<Real>&
  getBndryPhi(const int a_phase, const DataIndex& a_dit) const;

  /*!
    @brief Return vof iterator which iterates over single-phase cells on a specific phase and box. 
  */
  virtual VoFIterator&
  getSinglePhaseVofs(const int a_phase, const DataIndex& a_dit) const;

  /*!
    @brief Return vof iterator which iterates over multi-phase cells on a specified phase and box. 
  */
  virtual VoFIterator&
  getMultiPhaseVofs(const int a_phase, const DataIndex& a_dit) const;

  /*!
    @brief Get stencils for computing dphi/dn at the boundary
  */
  virtual const LayoutData<MFInterfaceFAB<VoFStencil>>&
  getGradPhiStencils() const noexcept;

  /*!
    @brief Get constant weights involved when computing dphi/dn at the boundary
  */
  virtual const LayoutData<MFInterfaceFAB<Real>>&
  getGradPhiWeights() const noexcept;

  /*!
    @brief Set everything to zero. This is a debugging function.
  */
  virtual void
  resetBC() const;

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
    @brief Data centering
  */
  Location::Cell m_dataLocation;

  /*!
    @brief Grids
  */
  MFLevelGrid m_mflg;

  /*!
    @brief B-coefficient
  */
  RefCountedPtr<LevelData<MFBaseIVFAB>> m_Bcoef;

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
  mutable LayoutData<MFInterfaceFAB<Real>> m_boundaryPhi;

  /*!
    @brief Iterators for single-phase cells. The integer index is the phase. 
  */
  mutable std::map<int, std::shared_ptr<LayoutData<VoFIterator>>> m_singlePhaseVofs;

  /*!
    @brief Iterators for multi-phase cells. The integer index is the phase. 
  */
  mutable std::map<int, std::shared_ptr<LayoutData<VoFIterator>>> m_multiPhaseVofs;

  /*!
    @brief IntVectSet over which we run the "averaging".
  */
  LayoutData<IntVectSet> m_ivs;

  /*!
    @brief Raw stencils for dphi/dn
  */
  LayoutData<MFInterfaceFAB<VoFStencil>> m_gradPhiStencils;

  /*!
    @brief Raw weights for dphi/dn
  */
  LayoutData<MFInterfaceFAB<Real>> m_gradPhiWeights;

  /*!
    @brief Average stencil
  */
  LayoutData<MFInterfaceFAB<VoFStencil>> m_avgStencils;

  /*!
    @brief Denominator 1/(bp*wp + bq*wq) for all interface cells. 
  */
  LayoutData<MFInterfaceFAB<Real>> m_denom;

  /*!
    @brief Average weights
  */
  LayoutData<MFInterfaceFAB<Real>> m_avgWeights;

  /*!
    @brief Agg stencils for making matching go faster
  */
  LayoutData<RefCountedPtr<AggStencil<EBCellFAB, BaseIVFAB<Real>>>> m_aggStencils[2];

  /*!
    @brief Define function. Builds stencils. 
  */
  void
  defineStencils();

  /*!
    @brief Define function. Builds stencils. 
  */
  void
  defineIterators();

  /*!
    @brief Brief the average stencils
  */
  void
  buildAverageStencils();

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
  virtual bool
  getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
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
  inline bool
  isStencilValidCF(const VoFStencil& a_stencil, const DataIndex& a_dit) const;

  /*!
    @brief Apply a stencil and return the result
    @param[in] a_stencil Stencil to apply
    @param[in] a_phi     Data
    @returns Returns S(phi), S being the stencil. 
  */
  inline Real
  applyStencil(const VoFStencil& a_stencil, const EBCellFAB& a_phi) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_MFHelmholtzJumpBCImplem.H>

#endif
