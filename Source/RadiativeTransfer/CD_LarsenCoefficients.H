/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LarsenCoefficients.H
  @brief  Declaration of a class that holds coefficients for Photon outflow boundary conditions
  @author Robert Marskar
*/

#ifndef CD_LarsenCoefficients_H
#define CD_LarsenCoefficients_H

// Std includes
#include <functional>

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_RobinCoefficients.H>
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that holds appropriate coefficients for Robin implementation of the Larsen boundary condition
*/
class LarsenCoefficients : public RobinCoefficients {
public:

  /*!
    @brief Alias for BC function
  */
  using BcFunction = std::function<Real(const RealVect a_pos, const Real a_time)>;

  /*!
    @brief Don't see why anyone would need this constructor
  */
  LarsenCoefficients() = delete;

  /*!
    @brief No copy construction
  */
  LarsenCoefficients(const LarsenCoefficients& a_other) = delete;

  /*!
    @brief No move construction
  */
  LarsenCoefficients(const LarsenCoefficients&& a_other) = delete;

  /*!
    @brief Constructor. Source term is set to zero. 
    @param[in] a_RtSpecies Radiative transfer species
    @param[in] a_r1 "First reflection coefficient"
    @param[in] a_r2 "Second reflection coefficient"
    @details Larsen coefficients (in Robin format) are set based on kappa and reflection coefficients
  */
  LarsenCoefficients(const RefCountedPtr<RtSpecies>& a_RtSpecies,
		     const Real                      a_r1,
		     const Real                      a_r2);

  /*!
    @brief Constructor. Source term is set to zero. 
    @param[in] a_RtSpecies Radiative transfer species
    @param[in] a_r1        "First reflection coefficient"
    @param[in] a_r2        "Second reflection coefficient"
    @param[in] a_bcFunction Source term on boundary
    @details Larsen coefficients (in Robin format) are set based on kappa and reflection coefficients. The last parameter provides
    an option to include sources on the boundaries. To include time-dependence, a_bcFunction should capture an external time by reference.
  */
  LarsenCoefficients(const RefCountedPtr<RtSpecies>& a_RtSpecies,
		     const Real                      a_r1,
		     const Real                      a_r2,
		     const BcFunction&               a_bcFunction);

  /*!
    @brief Destructor
  */
  ~LarsenCoefficients();

  /*!
    @brief No move assignment
  */
  LarsenCoefficients& operator=(const LarsenCoefficients& a_other) = delete;

  /*!
    @brief No copy assignment
  */
  LarsenCoefficients& operator=(const LarsenCoefficients&& a_other) = delete; 

  /*!
    @brief aco. Equal to 0.5*3*kappa*kappa*(1 + 3*r2)/(1 - 2*r1)
    @param[in] a_pos Position
  */
  virtual Real aco(const RealVect a_pos) const override;

  /*!
    @brief bco. Equal to -kappa
    @param[in] a_pos Position
  */
  virtual Real bco(const RealVect a_pos) const override;

  /*!
    @brief Surface source term. 
    @param[in] a_pos Position
  */
  virtual Real rhs(const RealVect a_pos) const override;
    
protected:

  /*!
    @brief Radiative transfer species
  */
  RefCountedPtr<RtSpecies> m_RtSpecies;

  /*!
    @brief "First reflection coefficient"
  */
  Real m_reflectionCoefficientOne;

  /*!
    @brief "Second reflection coefficient"
  */
  Real m_reflectionCoefficientTwo;

  /*!
    @brief Boundary source term
  */
  BcFunction m_bcFunction;
};

#include <CD_NamespaceFooter.H>

#endif
