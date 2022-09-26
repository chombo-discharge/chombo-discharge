/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBC.H
  @brief  Declaration of a domain boundary condition class for EBHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDomainBC_H
#define CD_EBHelmholtzDomainBC_H

// Chombo includes
#include <EBLevelGrid.H>
#include <RefCountedPtr.H>
#include <LevelData.H>
#include <EBFluxFAB.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing domain boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzDomainBC
{
public:
  /*!
    @brief Constructor. 
  */
  EBHelmholtzDomainBC();

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzDomainBC();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzDomainBC(const EBHelmholtzDomainBC& a_other) = delete;

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzDomainBC(const EBHelmholtzDomainBC&& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzDomainBC&
  operator=(const EBHelmholtzDomainBC& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzDomainBC&
  operator=(const EBHelmholtzDomainBC&& a_other) = delete;

  /*!
    @brief Define function
    @param[in] a_eblg   Level grids
    @param[in] a_Bcoef  B-coefficient in the Helmholtz equation
    @param[in] a_probLo Lower-left corner of computational domain
    @param[in] a_dx     Grid resolution
    @note This sets the internals for the BC object. This is called by EBHelmholtzOp when it defines stencils. 
  */
  virtual void
  define(const Location::Cell                       a_dataLocation,
         const EBLevelGrid&                         a_eblg,
         const RefCountedPtr<LevelData<EBFluxFAB>>& a_Bcoef,
         const RealVect&                            a_probLo,
         const Real                                 a_dx);

  /*!
    @brief Get face flux. This is for regular cells. 
    @param[inout] a_faceFlux       Face-centered flux. 
    @param[in]    a_phi            Cell-centered data
    @param[in]    a_Bcoef          Face-centered B-coefficient
    @param[in]    a_dir            Coordinate direction
    @param[in]    a_side           High or low side
    @param[in]    a_dit            Data index
    @param[in]    a_useHomogeneous Homogeneous physical BCs or not. 
    @note a_faceFlux is defined in the first strip of interior cells on the (dir,side) edge of the domain. 
  */
  virtual void
  getFaceFlux(BaseFab<Real>&        a_faceFlux,
              const BaseFab<Real>&  a_phi,
              const BaseFab<Real>&  a_Bcoef,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const DataIndex&      a_dit,
              const bool            a_useHomogeneous) const = 0;

  /*!
    @brief Get face flux. This is for irregular cells. 
    @param[in] a_vof            Grid cell
    @param[in] a_phi            Cell-centered data
    @param[in] a_Bcoef          Face-centered B-coefficient
    @param[in] a_dir            Coordinate direction
    @param[in] a_side           High or low side
    @param[in] a_dit            Data index
    @param[in] a_useHomogeneous Homogeneous physical BCs or not. 
  */
  virtual Real
  getFaceFlux(const VolIndex&       a_vof,
              const EBCellFAB&      a_phi,
              const EBFaceFAB&      a_Bcoef,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const DataIndex&      a_dit,
              const bool            a_useHomogeneous) const = 0;

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
    @brief Data centering. 
  */
  Location::Cell m_dataLocation;

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
    @brief B-coefficent. 
  */
  RefCountedPtr<LevelData<EBFluxFAB>> m_Bcoef;

  /*!
    @brief Returns the cell-centered position at the boundary. 
    @param[in] a_iv   Grid cell
    @param[in] a_dir  Face coordinate direction
    @param[in] a_side Lo/Hi side
  */
  inline RealVect
  getBoundaryPosition(const IntVect& a_iv, const int& a_dir, const Side::LoHiSide& a_side) const;

  /*!
    @brief Utility function which multiplies the flux (or actually, dphi/dn) by the B-coefficient. This is used 
    by the implementation classes when they compute the regular flux at the domain edges/faces. 
    @param[inout] a_flux    On input, contains dphi/dn. On output, contains B*dphi/dn
    @param[in]    a_bco     Helmholtz B-coefficient
    @param[in]    a_dir     Coordinate direction
    @param[in]    a_side    Side
    @note This routine is weird because the flux data holder in the BC classes lives on the cell centers next to the domain
    face/edge. So we need to shift the box over when we multiply. Because both Dirichlet, Neumann, and the Robin Bcs need to do this,
    I'm putting this routine in the parent class. 
  */
  void
  multiplyByBcoef(BaseFab<Real>&       a_flux,
                  const BaseFab<Real>& a_bco,
                  const int            a_dir,
                  const Side::LoHiSide a_side) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBHelmholtzDomainBCImplem.H>

#endif
