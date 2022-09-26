/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEddingtonDomainBC.H
  @brief  Declaration of a domain boundary condition class for EBHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzEddingtonDomainBC_H
#define CD_EBHelmholtzEddingtonDomainBC_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_EddingtonSP1DomainBc.H>
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for passing EddingtonSP1 domain boundary conditions into EBHelmholtzOp. 
  @details This calss is a wrapper for Dirichlet/Neumann boundary conditions where the values/functions are passed in through EddingtonDomainBc
  @note For this class to operate correctly with time-dependent boundary conditions, the user must ensure that
  the time dependence is captured by reference. The time argument in ElectrostaticDomainBc::BcFunction by this class,
  so whatever function is put in there, it must capture time by reference. 
*/
class EBHelmholtzEddingtonSP1DomainBC : public EBHelmholtzDomainBC
{
public:
  /*!
    @brief Not allowed, use the full constructor. 
  */
  EBHelmholtzEddingtonSP1DomainBC() = delete;

  /*!
    @brief Full constructor
    @param[in] a_eddingtonBCs EddingtonSP1 boundary conditions on domain sides
    @param[in] a_species      Radiative transfer species
    @param[in] a_r1           "First reflection coefficient"
    @param[in] a_r2           "Second reflection coefficient"
  */
  EBHelmholtzEddingtonSP1DomainBC(const EddingtonSP1DomainBc&     a_eddingtonBCs,
                                  const RefCountedPtr<RtSpecies>& a_species,
                                  const Real                      a_r1,
                                  const Real                      a_r2);

  /*!
    @brief Disallowed -- don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzEddingtonSP1DomainBC(const EBHelmholtzEddingtonSP1DomainBC& a_other) = delete;

  /*!
    @brief Disallowed -- don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzEddingtonSP1DomainBC(const EBHelmholtzEddingtonSP1DomainBC&& a_other) = delete;

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzEddingtonSP1DomainBC();

  /*!
    @brief Disallowed - don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzEddingtonSP1DomainBC&
  operator=(const EBHelmholtzEddingtonSP1DomainBC& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
    @param[in] a_other Other operator
  */
  EBHelmholtzEddingtonSP1DomainBC&
  operator=(const EBHelmholtzEddingtonSP1DomainBC&& a_other) = delete;

  /*!
    @brief Define function
    @param[in] a_dataLocation Data location
    @param[in] a_eblg         Level grids
    @param[in] a_probLo       Lower-left corner of computational domain
    @param[in] a_dx           Grid resolution
    @note Overriding because this is a wrapper class for Dirichlet/Neumann, and those objects need defining. 
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
    @brief Eddington domain boundary conditions
  */
  EddingtonSP1DomainBc m_eddingtonBCs;

  /*!
    @brief Radiative transfer species
  */
  RefCountedPtr<RtSpecies> m_species;

  /*!
    @brief "First reflection coefficient"
  */
  Real m_r1;

  /*!
    @brief "Second reflection coefficient"
  */
  Real m_r2;

  /*!
    @brief BC objects on each domain side. 
    @details These are instantiated in the constructor where we construct Neumann/Dirichlet on each domain side. 
  */
  std::map<EddingtonSP1DomainBc::DomainSide, std::shared_ptr<EBHelmholtzDomainBC>> m_bcObjects;
};

#include <CD_NamespaceFooter.H>

#endif
