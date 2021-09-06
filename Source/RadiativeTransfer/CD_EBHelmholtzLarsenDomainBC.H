/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzLarsenDomainBC.H
  @brief  Declaration of a Larsen boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzLarsenDomainBC_H
#define CD_EBHelmholtzLarsenDomainBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_RtSpecies.H>
#include <CD_EBHelmholtzRobinDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Larsen domain boundary conditions into EBHelmholtzOp. 
  @details Larsen bcs are in the form A*phi + B*dphi/dn = C. We impose dphi/dn = (C - A*phi)/B where phi is extrapolated to the EB to first order. 
  @note The B-coefficient is not the same as in the Helmholtz equation...
  @note In our notation, n points into the computational domain. 
*/
class EBHelmholtzLarsenDomainBC : public EBHelmholtzRobinDomainBC {
public:

  /*!
    @brief Alias for source function
  */
  using SourceFunction = std::function<Real(const RealVect& a_position)>;

  /*!
    @brief Disallowed default constructor. 
  */
  EBHelmholtzLarsenDomainBC() = delete;

  /*!
    @brief Constructor. Source term is set to zero. 
    @param[in] a_species Radiative transfer species
    @param[in] a_r1      "First reflection coefficient"
    @param[in] a_r2      "Second reflection coefficient"
    @param[in] a_source   Source term
    @details Larsen coefficients (in Robin format) are set based on kappa and reflection coefficients. 
  */
  EBHelmholtzLarsenDomainBC(const RefCountedPtr<RtSpecies>& a_species,
			    const Real                      a_r1,
			    const Real                      a_r2,
			    const SourceFunction            a_source);

  /*!
    @brief Constructor. Source term is set to zero. 
    @param[in] a_species Radiative transfer species
    @param[in] a_r1      "First reflection coefficient"
    @param[in] a_r2      "Second reflection coefficient"
    @details Larsen coefficients (in Robin format) are set based on kappa and reflection coefficients. 
  */
  EBHelmholtzLarsenDomainBC(const RefCountedPtr<RtSpecies>& a_species,
			    const Real                      a_r1,
			    const Real                      a_r2);


  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzLarsenDomainBC();

protected:

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
    @brief Source term
  */
  SourceFunction m_source;

  /*!
    @brief Set Robin coefficients
  */
  void setRobinFunctions();
};

#include <CD_NamespaceFooter.H>

#endif
