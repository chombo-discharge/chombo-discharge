/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBC.H
  @brief  Declaration of a boundary condition class for EBHelmholtzOp
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
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing EB boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzDomainBC {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzDomainBC();

  /*!
    @brief Destructor
  */
  ~EBHelmholtzDomainBC();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzDomainBC(const EBHelmholtzDomainBC& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzDomainBC& operator=(const EBHelmholtzDomainBC& a_other) = delete;

  /*!
    @brief Define function
    @param[in] a_eblg   Level grids
    @param[in] a_Bcoef  B-coefficient for Helmholtz equation
    @param[in] a_probLo Lower-left corner of computational domain
    @param[in] a_dx     Grid resolution
  */
  void define(const EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<EBFluxFAB> >& a_Bcoef, const RealVect& m_probLo, const Real a_dx);

  /*!
    @brief Get face flux. This is for regular cells
    @param[inout] a_faceFlux       Face-centered flux. 
    @param[in]    a_phi            Cell-centered data
    @param[in]    a_dir            Coordinate direction
    @param[in]    a_side           High or low side
    @param[in]    a_dit            Data index
    @param[in]    a_useHomogeneous Homogeneous physical BCs or not. 
    @note a_faceFlux is defined in the first strip of interior cells on the (dir,side) edge of the domain. 
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
			   const BaseFab<Real>&  a_phi,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const DataIndex&      a_dit,
			   const bool            a_useHomogeneous) const = 0;

  /*!
    @brief Get face flux. This is for irregular cells. 
    @param[in] a_vof            Grid cell
    @param[in] a_phi            Cell-centered data
    @param[in] a_dir            Coordinate direction
    @param[in] a_side           High or low side
    @param[in] a_dit            Data index
    @param[in] a_useHomogeneous Homogeneous physical BCs or not. 
  */
  virtual Real getFaceFlux(const VolIndex&       a_vof,
			   const EBCellFAB&      a_phi,
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
  RefCountedPtr<LevelData<EBFluxFAB> > m_Bcoef;

  /*!
    @brief Returns the cell-centered position at the boundary. 
  */
  RealVect getBoundaryPosition(const IntVect& a_iv, const int& a_dir, const Side::LoHiSide& a_side) const;
};

#include <CD_NamespaceFooter.H>

#endif
