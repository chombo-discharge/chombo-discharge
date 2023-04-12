/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBC.H
  @brief  Declaration of a base boundary condition class for EB boundary conditions in EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzEBBC_H
#define CD_EBHelmholtzEBBC_H

// Chombo includes
#include <EBLevelGrid.H>
#include <RefCountedPtr.H>
#include <LevelData.H>
#include <BaseIVFAB.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing EB boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzEBBC
{
public:
  /*!
    @brief Default constructor. 
  */
  EBHelmholtzEBBC();

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzEBBC();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzEBBC(const EBHelmholtzEBBC& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzEBBC&
  operator=(const EBHelmholtzEBBC& a_other) = delete;

  /*!
    @brief Define function that is called by EBHelmholtzOp
    @param[in] a_dataLocation    Specification of data centering
    @param[in] a_eblg            Grid where this object is defined. 
    @param[in] a_eblgFiCo        Refinement of this grid.
    @param[in] a_probLo          Lower-left corner of computational domain. Required if BC objects needs physical position. 
    @param[in] a_hasFineAMRLevel Has a finer AMR level or not. False if it doesn't, or if this level is an MG level
    @param[in] a_isMGLevel       Is a multigrid level or not. 
    @param[in] a_dx              Grid resolution. 
    @param[in] a_refToFine       Refinement factor to fine level.
    @param[in] a_ghostCF         Number of ghost cells that are filled across CF interface. 
  */
  void
  define(const Location::Cell                           a_dataLocation,
         const EBLevelGrid&                             a_eblg,
         const EBLevelGrid&                             a_eblgFiCo,
         const RefCountedPtr<LevelData<BaseFab<bool>>>& a_amrValidCells,
         const RealVect&                                a_probLo,
         const Real&                                    a_dx,
         const bool                                     a_hasFineAMRLevel,
         const bool                                     a_isMGLevel,
         const int                                      a_refToFine,
         const int                                      a_ghostCF);

  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @note The EB flux can be represented by dphi/dn = wb * phi_b + sum (wi * phi_i). 
    should input the part of dphi/dn that doesn't use the interior points. 
    @param[inout] a_vofit             Iterator for cut-cells
    @param[in]    a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_Bcoef             EB-centered B-coefficient
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void
  applyEBFluxRelax(VoFIterator&           a_vofit,
                   EBCellFAB&             a_Lphi,
                   const EBCellFAB&       a_phi,
                   const BaseIVFAB<Real>& a_Bcoef,
                   const DataIndex&       a_dit,
                   const Real&            a_beta,
                   const bool&            a_homogeneousPhysBC) const = 0;

  /*!
    @brief Apply full flux through the EB
  */
  virtual void
  applyEBFluxResid(VoFIterator&           a_vofit,
                   EBCellFAB&             a_Lphi,
                   const EBCellFAB&       a_phi,
                   const EBCellFAB* const a_phiFine,
                   const BaseIVFAB<Real>& a_Bcoef,
                   const DataIndex&       a_dit,
                   const Real&            a_beta,
                   const bool&            a_homogeneousPhysBC) const;

  /*!
    @brief Get the stencil for computing the finite-volume approximation to kappa*Div(F) on this grid level. 
    @details The finite volume approximation to kappa*Div(F) is sum(Fluxes)/dx. This routine should return the flux on the EB, i.e. (B*dphi/dn)/dx. The derivative
             can be written dphi/dn = w0*phi0 + sum(wi*phi(i)) where phi0 is a boundary weight and phi(i) are grid points. This routine should return the stencil
	     B*wi*phi(i).
    @note    Multiplication by beta is done in the operator.
  */
  const LayoutData<BaseIVFAB<VoFStencil>>&
  getGradPhiRelaxStencils() const noexcept;

protected:
  /*!
    @brief Component that everything is defined for. Always have m_comp = 0
  */
  constexpr static int m_comp = 0;

  /*!
    @brief Number of components. Always have m_nComp = 1
  */
  constexpr static int m_nComp = 1;

  /*!
    @brief Data centering
  */
  Location::Cell m_dataLocation;

  /*!
    @brief Has finer AMR level or not. This is false for all deeper multigrid levels. 
  */
  bool m_hasFineAMRLevel;

  /*!
    @brief If false, this level is an AMR level and fluxes might need stencils on the fine level also. 
  */
  bool m_isMGLevel;

  /*!
    @brief Refinement factor to fine level (if there is one)
    @details This is used when constructing stencils for the EB fluxes, which should not reach into invalid cells.
  */
  int m_refToFine;

  /*!
    @brief Number of ghost cells that were filled across CF interface
    @details This is passed in from the operator, who knows how many ghost cells were filled by the interpolator. This
    is used so that the stencil-defining objects don't create stencils that reach into invalid ghost cells. Very useful.
  */
  int m_ghostCF;

  /*!
    @brief Grid resolution
  */
  Real m_dx;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;

  /*!
    @brief Level grid
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Refinement of m_eblg 
    @details This should be define donly if m_hasFineAMRLevel and !m_isMGLevel. 
  */
  EBLevelGrid m_eblgFiCo;

  /*!
    @brief If this is an actual AMR level, this contains a list of the cells on this level that are not 
    covered by a finer level (true if not covered, false if covered).
  */
  RefCountedPtr<LevelData<BaseFab<bool>>> m_amrValidCells;

  /*!
    @brief Stencils for computing the flux on a single level. This is a single-level object used for 
    the multigrid relaxation stage. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_gradPhiRelaxStencils;

  /*!
    @brief User define function. 
  */
  virtual void
  define() = 0;

  /*!
    @brief Check if stencil is valid. 
    @details This routine makes sure that the stencil a_stencil does not reach over the CF interface and into ghost cells that do not have valid data.
    @param[in] a_stencil Stencil
    @param[in] a_dit     Data index
  */
  inline bool
  isStencilValidCF(const VoFStencil& a_stencil, const DataIndex& a_dit) const;

  /*!
    @brief Returns physical position at the boundary
  */
  inline RealVect
  getBoundaryPosition(const VolIndex& a_vof, const DataIndex& a_dit) const;

  /*!
    @brief Apply stencil to data holder and return result
  */
  inline Real
  applyStencil(const VoFStencil& a_stencil, const EBCellFAB& a_phi) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBHelmholtzEBBCImplem.H>

#endif
