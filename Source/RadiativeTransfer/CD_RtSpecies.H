/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtSpecies.H
  @brief  Declaration of a class which supplies a user interface to radiative transfer code. 
  @author Robert Marskar
*/

#ifndef CD_RtSpecies_H
#define CD_RtSpecies_H

// Std includes
#include <functional>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a class that defines an interface to radiative transfer solvers. 
  @details This is a lightweight class - it just holds a name for the group and the kappa for the radiative transfer species.
*/
class RtSpecies {
public:
  /*!
    @brief Weak constructor, must subsequently call define
  */
  RtSpecies();

  /*!
    @brief Constructor, constant kappa
    @param[in] a_name  Species name
    @param[in] a_kappa Kappa (inverse absorption length)
  */
  RtSpecies(const std::string a_name, const Real a_kappa);

  /*!
    @brief Constructor, spatially dependent kappa
    @param[in] a_name  Species name
    @param[in] a_kappa Kappa (inverse absorption length)
  */
  RtSpecies(const std::string a_name, const std::function<Real(const RealVect a_pos)>& a_kappa);

  /*!
    @brief Define function. This sets a constant absorption length.
    @param[in] a_name  Species name
    @param[in] a_kappa Kappa (inverse absorption length)
  */
  virtual void define(const std::string a_name, const Real a_kappa);

  /*!
    @brief Define function. This sets a spatially varying absorption length.
    @param[in] a_name  Species name
    @param[in] a_kappa Kappa (inverse absorption length)
  */
  virtual void define(const std::string a_name, const std::function<Real(const RealVect a_pos)>& a_kappa);

  /*!
    @brief Get species name.
  */
  virtual std::string getName() const;

  /*!
    @brief Get kappa at position. 
    @param[in] a_pos Physical coordinates. 
  */
  virtual Real getKappa(const RealVect a_pos) const;

  /*!
    @brief Get scattering coefficient. 
    @param[in] a_pos Physical coordinates. 
    @note This is not used anywhere (yet).
  */
  virtual Real getScatteringCoefficient(const RealVect a_pos) const;
	      
protected:

  /*!
    @brief Group name
  */
  std::string m_name;

  /*!
    @brief Switch for using constant value or function
  */
  bool m_constant;

  /*!
    @brief Constant kappa
  */
  Real m_kappa;

  /*!
    @brief Scatter coefficient
  */
  Real m_scatter;

  /*!
    @brief Variable kappa
  */
  std::function<Real(const RealVect a_pos)> m_varKappa;

  /*!
    @brief Variable kappa
  */
  std::function<Real(const RealVect a_pos)> m_varScatter;
};

#include <CD_NamespaceFooter.H>

#endif
  
