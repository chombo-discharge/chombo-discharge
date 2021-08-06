/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBHelmholtzLarsenDomainBCFactory.H
  @brief  Declaration of a Electrostatic boundary condition class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzLarsenDomainBCFactory_H
#define CD_EBHelmholtzLarsenDomainBCFactory_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_EBHelmholtzRobinDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Larsen boundary conditions 
*/
class EBHelmholtzLarsenDomainBCFactory : public EBHelmholtzRobinDomainBCFactory {
public:

  /*!
    @brief Alias for source function
  */
  using SourceFunction = std::function<Real(const RealVect& a_position, const Real a_time)>;

  /*!
    @brief Disallowed weak constructor
  */
  EBHelmholtzLarsenDomainBCFactory() = delete;

  /*!
    @brief Constructor. Source term is set to zero. 
    @param[in] a_species Radiative transfer species
    @param[in] a_r1      "First reflection coefficient"
    @param[in] a_r2      "Second reflection coefficient"
    @param[in] a_source   Source term
    @details Larsen coefficients (in Robin format) are set based on kappa and reflection coefficients. 
  */
  EBHelmholtzLarsenDomainBCFactory(const RefCountedPtr<RtSpecies>& a_species,
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
  EBHelmholtzLarsenDomainBCFactory(const RefCountedPtr<RtSpecies>& a_species,
				   const Real                      a_r1,
				   const Real                      a_r2);

  /*!
    @brief No copy construction
  */
  EBHelmholtzLarsenDomainBCFactory(const EBHelmholtzLarsenDomainBCFactory& a_other) = delete;

  /*!
    @brief No move construction
  */
  EBHelmholtzLarsenDomainBCFactory(const EBHelmholtzLarsenDomainBCFactory&& a_other) = delete;

  /*!
    @brief Default constructor
  */
  ~EBHelmholtzLarsenDomainBCFactory();

  /*!
    @brief No copy assignment
  */
  EBHelmholtzLarsenDomainBCFactory& operator=(const EBHelmholtzLarsenDomainBCFactory& a_other) = delete;

  /*!
    @brief No move assignment
  */
  EBHelmholtzLarsenDomainBCFactory& operator=(const EBHelmholtzLarsenDomainBCFactory&& a_other) = delete;

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
  void setRobinCoefficients();
};

#include <CD_NamespaceFooter.H>

#endif
