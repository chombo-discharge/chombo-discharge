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
  */
  virtual void define(const EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<EBFluxFAB> >& a_Bcoef, const RealVect& m_probLo, const Real a_dx);

  /*!
    @brief Get face flux. This is for regular cells
  */
  virtual void getFaceFlux(FArrayBox& a_faceFlux, const FArrayBox& a_phi, const int& dir, const Side::LoHiSide& a_side) = 0;

  /*!
    @brief Get face flux. This is for irregular cells. 
  */
  virtual void getFaceFlux(Real&                 a_faceFlux,
			   const VolIndex&       a_vof,
			   const EBCellFAB&      a_phi,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const DataIndex&      a_dit) = 0;			   

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
  RealVect getBoundaryPosition(const VolIndex& a_vof, const int& a_dir, const Side::LoHiSide& a_side) const;
};

#include <CD_NamespaceFooter.H>

#endif
