/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzElectrostaticDomainBC.H
  @brief  Declaration of a domain boundary condition class for EBHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzElectrostaticDomainBC_H
#define CD_EBHelmholtzElectrostaticDomainBC_H

// Our includes
#include <CD_ElectrostaticDomainBc.H>
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing domain boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzElectrostaticDomainBC : public EBHelmholtzDomainBC {
public:

  /*!
    @brief Not allowed, use the full constructor. 
  */
  EBHelmholtzElectrostaticDomainBC() = delete;

  /*!
    @brief Not allowed, use the full constructor. 
  */
  EBHelmholtzElectrostaticDomainBC(const ElectrostaticDomainBc& a_electrostaticBCs);

  /*!
    @brief Destructor
  */
  ~EBHelmholtzElectrostaticDomainBC();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzElectrostaticDomainBC(const EBHelmholtzElectrostaticDomainBC& a_other) = delete;

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzElectrostaticDomainBC(const EBHelmholtzElectrostaticDomainBC&& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzElectrostaticDomainBC& operator=(const EBHelmholtzElectrostaticDomainBC& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzElectrostaticDomainBC& operator=(const EBHelmholtzElectrostaticDomainBC&& a_other) = delete;

  /*!
    @brief Get face flux. This is for regular cells. 
    @param[inout] a_faceFlux       Face-centered flux. 
    @param[in]    a_phi            Cell-centered data
    @param[in]    a_dir            Coordinate direction
    @param[in]    a_side           High or low side
    @param[in]    a_dit            Data index
    @param[in]    a_useHomogeneous Homogeneous physical BCs or not. 
    @note a_faceFlux is defined in the first strip of interior cells on the (dir,side) edge of the domain. 
  */
  void getFaceFlux(BaseFab<Real>&        a_faceFlux,
		   const BaseFab<Real>&  a_phi,
		   const int&            a_dir,
		   const Side::LoHiSide& a_side,
		   const DataIndex&      a_dit,
			   const bool            a_useHomogeneous) const override final;

  /*!
    @brief Get face flux. This is for irregular cells. 
    @param[in] a_vof            Grid cell
    @param[in] a_phi            Cell-centered data
    @param[in] a_dir            Coordinate direction
    @param[in] a_side           High or low side
    @param[in] a_dit            Data index
    @param[in] a_useHomogeneous Homogeneous physical BCs or not. 
  */
  Real getFaceFlux(const VolIndex&       a_vof,
		   const EBCellFAB&      a_phi,
		   const int&            a_dir,
		   const Side::LoHiSide& a_side,
		   const DataIndex&      a_dit,
		   const bool            a_useHomogeneous) const override final;

protected:
  
  /*!
    @brief Functions for electrostatic domain boundary conditions
  */
  ElectrostaticDomainBc m_electrostaticBCs;
  
};

#include <CD_NamespaceFooter.H>

#endif
