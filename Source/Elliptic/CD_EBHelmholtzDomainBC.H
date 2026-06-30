/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  @file   CD_EBHelmholtzDomainBC.H
  @brief  Declaration of a domain boundary condition class for EBHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_EBHELMHOLTZDOMAINBC_H
#define CD_EBHELMHOLTZDOMAINBC_H

// Chombo includes
#include <EBLevelGrid.H>
#include <RefCountedPtr.H>
#include <LevelData.H>
#include <EBFluxFAB.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Base class for passing domain boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzDomainBC
{
public:
  /**
    @brief Constructor. 
  */
  EBHelmholtzDomainBC();

  /**
    @brief Destructor
  */
  virtual ~EBHelmholtzDomainBC();

  /**
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzDomainBC(const EBHelmholtzDomainBC& a_other) = delete;

  /**
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzDomainBC(const EBHelmholtzDomainBC&& a_other) = delete;

  /**
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzDomainBC&
  operator=(const EBHelmholtzDomainBC& a_other) = delete;

  /**
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzDomainBC&
  operator=(const EBHelmholtzDomainBC&& a_other) = delete;

  /**
    @brief Define function
    @param[in] a_eblg   Level grids
    @param[in] a_probLo Lower-left corner of computational domain
    @param[in] a_dx     Grid resolution
    @note This sets the internals for the BC object. This is called by EBHelmholtzOp when it defines stencils. 
    @param[in] a_dataLocation Description
  */
  virtual void
  define(const Location::Cell a_dataLocation, const EBLevelGrid& a_eblg, const RealVect& a_probLo, const Real a_dx);

  /**
    @brief Get face flux. This is for regular cells. 
    @param[in,out] a_faceFlux       Face-centered flux. 
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

  /**
    @brief Get face flux. This is for irregular cells. 
    @param[in] a_vof            Grid cell
    @param[in] a_phi            Cell-centered data
    @param[in] a_Bcoef          Face-centered B-coefficient
    @param[in] a_dir            Coordinate direction
    @param[in] a_side           High or low side
    @param[in] a_dit            Data index
    @param[in] a_useHomogeneous Homogeneous physical BCs or not. 
    @return Face flux
  */
  virtual Real
  getFaceFlux(const VolIndex&       a_vof,
              const EBCellFAB&      a_phi,
              const EBFaceFAB&      a_Bcoef,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const DataIndex&      a_dit,
              const bool            a_useHomogeneous) const = 0;

  /**
    @brief Diagonal-weight multiplier for this boundary face's contribution to the relaxation diagonal.
    @details EBHelmholtzOp::computeDiagWeight builds the relaxation diagonal (m_betaDiagWeight) by taking the
             interior stencil diagonal and, for each domain-boundary face, subtracting this weight times the
             area-fraction-weighted face B-coefficient, i.e. betaWeight -= getDiagWeight() * sum(areaFrac * B/dx^2).
             This weight is therefore the coefficient with which the boundary cell's own value phi_i enters the
             discretized boundary flux returned by getFaceFlux() -- it is how strongly the boundary flux feeds back
             onto the diagonal. It is dimensionless and MUST be kept consistent with this subclass's getFaceFlux()
             discretization: the two encode the same boundary stencil, so if that discretization changes (e.g. a
             higher-order or centroid-based flux) this value must change with it.

             How to choose it: read off the coefficient of phi_i in getFaceFlux(). For the shipped discretizations:
               - Dirichlet -> 2. The boundary value sits half a cell from the cell centre, so the one-sided difference
                 (phi_b - phi_i)/(dx/2) carries a factor of 2 on phi_i.
               - Neumann   -> 0. The flux is prescribed independently of phi_i, so it does not feed the diagonal.

             This is pure virtual on purpose: every concrete domain BC must state its diagonal contribution
             explicitly, so a new BC cannot silently inherit a wrong value. Note this only affects the relaxation
             (smoother) diagonal, not the operator apply, so an inconsistent value degrades convergence rather than
             changing the converged solution.
    @param[in] a_dir  Coordinate direction of the face.
    @param[in] a_side High or low side.
    @return Multiplier applied to B·area/dx² when building the diagonal relaxation weight.
  */
  virtual Real
  getDiagWeight(const int a_dir, const Side::LoHiSide a_side) const = 0;

protected:
  /**
    @brief Component that everything is defined for. Always have m_comp = 0
  */
  constexpr static int m_comp = 0;

  /**
    @brief Number of components. Always have m_nComp = 1
  */
  constexpr static int m_nComp = 1;

  /**
    @brief Data centering. 
  */
  Location::Cell m_dataLocation;

  /**
    @brief Grid resolution
  */
  Real m_dx;

  /**
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;

  /**
    @brief Level grid
  */
  EBLevelGrid m_eblg;

  /**
    @brief Returns the cell-centered position at the boundary. 
    @param[in] a_iv   Grid cell
    @param[in] a_dir  Face coordinate direction
    @param[in] a_side Lo/Hi side
    @return Boundary position
  */
  inline RealVect
  getBoundaryPosition(const IntVect& a_iv, const int& a_dir, const Side::LoHiSide& a_side) const;

  /**
    @brief Utility function which multiplies the flux (or actually, dphi/dn) by the B-coefficient. This is used 
    by the implementation classes when they compute the regular flux at the domain edges/faces. 
    @param[in,out] a_flux    On input, contains dphi/dn. On output, contains B*dphi/dn
    @param[in]    a_bco     Helmholtz B-coefficient
    @param[in]    a_dir     Coordinate direction
    @param[in]    a_side    Side
    @note This routine is weird because the flux data holder in the BC classes lives on the cell centers next to the domain
    face/edge. So we need to shift the box over when we multiply. Because both Dirichlet, Neumann, and the Robin Bcs need to do this,
    I'm putting this routine in the parent class. 
  */
  static void
  multiplyByBcoef(BaseFab<Real>& a_flux, const BaseFab<Real>& a_bco, const int a_dir, const Side::LoHiSide a_side);
};

#include <CD_NamespaceFooter.H>

#include <CD_EBHelmholtzDomainBCImplem.H>

#endif
