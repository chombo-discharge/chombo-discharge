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
  @brief Class for passing electrostatic domain boundary conditions into EBHelmholtzOp. 
  @details This calss is a wrapper for Dirichlet/Neumann boundary conditions where the values/functions are passed in through ElectrostaticDomainBc. This
  class will simply instantiate one such BC object for each edge/face in the computational domain, and call the appropriate BC function when EBHelmholtzOp
  asks for it. 
  @note For this class to operate correctly with time-dependent boundary conditions, the user must ensure that
  the time dependence is captured by reference. FieldSolverGMG already does this, by taking the user's boundary condition function
  and replacing the time-argument by a reference to FieldSolver::m_time. 
*/
class EBHelmholtzElectrostaticDomainBC : public EBHelmholtzDomainBC
{
public:
  /*!
    @brief Not allowed, use the full constructor. 
  */
  EBHelmholtzElectrostaticDomainBC() = delete;

  /*!
    @brief Full constructor.
    @param[in] a_electrostaticBCs Electrostatic domain boundary condition functions (in time and space)
  */
  EBHelmholtzElectrostaticDomainBC(const ElectrostaticDomainBc& a_electrostaticBCs);

  /*!
    @brief Disallowed -- don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzElectrostaticDomainBC(const EBHelmholtzElectrostaticDomainBC& a_other) = delete;

  /*!
    @brief Disallowed -- don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzElectrostaticDomainBC(const EBHelmholtzElectrostaticDomainBC&& a_other) = delete;

  /*!
    @brief Destructor
  */
  ~EBHelmholtzElectrostaticDomainBC();

  /*!
    @brief Disallowed - don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzElectrostaticDomainBC&
  operator=(const EBHelmholtzElectrostaticDomainBC& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzElectrostaticDomainBC&
  operator=(const EBHelmholtzElectrostaticDomainBC&& a_other) = delete;

  /*!
    @brief Define function
    @param[in] a_dataLocation Data location (either Location::Cell::Center or Location::Cell::Centroid). 
    @param[in] a_eblg         Level grids
    @param[in] a_probLo       Lower-left corner of computational domain
    @param[in] a_dx           Grid resolution
    @note This overrides EBHelmholtzDomainBC because EBHelmholtzElectrostaticDomainBC is just a wrapper class for Dirichlet/Neumann, and those objects need defining. 
  */
  void
  define(const Location::Cell a_dataLocation,
         const EBLevelGrid&   a_eblg,
         const RealVect&      a_probLo,
         const Real           a_dx) override final;

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
  void
  getFaceFlux(BaseFab<Real>&        a_faceFlux,
              const BaseFab<Real>&  a_phi,
              const BaseFab<Real>&  a_Bcoef,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const DataIndex&      a_dit,
              const bool            a_useHomogeneous) const override final;

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
  Real
  getFaceFlux(const VolIndex&       a_vof,
              const EBCellFAB&      a_phi,
              const EBFaceFAB&      a_Bcoef,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const DataIndex&      a_dit,
              const bool            a_useHomogeneous) const override final;

protected:
  /*!
    @brief Electrostatic domain boundary conditions
  */
  ElectrostaticDomainBc m_electrostaticBCs;

  /*!
    @brief BC objects on each domain side. 
    @details These are instantiated in the constructor where we construct Neumann/Dirichlet on each domain side. 
  */
  std::map<ElectrostaticDomainBc::DomainSide, std::shared_ptr<EBHelmholtzDomainBC>> m_bcObjects;
};

#include <CD_NamespaceFooter.H>

#endif
