/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBHelmholtzLarsenEBBCFactory.H
  @brief  Declaration of a Electrostatic boundary condition class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzLarsenEBBCFactory_H
#define CD_EBHelmholtzLarsenEBBCFactory_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_EBHelmholtzRobinEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Larsen boundary conditions 
*/
class EBHelmholtzLarsenEBBCFactory : public EBHelmholtzRobinEBBCFactory {
public:

  /*!
    @brief Disallowed weak constructor
  */
  EBHelmholtzLarsenEBBCFactory() = delete;

  /*!
    @brief Constructor. Source term is set to zero. 
    @param[in] a_species Radiative transfer species
    @param[in] a_r1      "First reflection coefficient"
    @param[in] a_r2      "Second reflection coefficient"
    @param[in] a_source  "Source term
    @details Larsen coefficients (in Robin format) are set based on kappa and reflection coefficients. 
  */
  EBHelmholtzLarsenEBBCFactory(const RefCountedPtr<RtSpecies>& a_species,
			       const Real                      a_r1,
			       const Real                      a_r2,
			       const Real                      a_source);

  /*!
    @brief No copy construction
  */
  EBHelmholtzLarsenEBBCFactory(const EBHelmholtzLarsenEBBCFactory& a_other) = delete;

  /*!
    @brief No move construction
  */
  EBHelmholtzLarsenEBBCFactory(const EBHelmholtzLarsenEBBCFactory&& a_other) = delete;

  /*!
    @brief Default constructor
  */
  ~EBHelmholtzLarsenEBBCFactory();

  /*!
    @brief No copy assignment
  */
  EBHelmholtzLarsenEBBCFactory& operator=(const EBHelmholtzLarsenEBBCFactory& a_other) = delete;

  /*!
    @brief No move assignment
  */
  EBHelmholtzLarsenEBBCFactory& operator=(const EBHelmholtzLarsenEBBCFactory&& a_other) = delete;

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
  Real m_source;

  /*!
    @brief Set Robin coefficients
  */
  void setRobinCoefficients();
};

#include <CD_NamespaceFooter.H>

#endif
