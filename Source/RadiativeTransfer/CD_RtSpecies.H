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

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a class that defines a group within the multigroup approximation
  @details This is a lightweight class - it just holds a name for the group and the kappa for the group
*/
class RtSpecies {
public:
  /*!
    @brief Weak constructor, must subsequently call define
  */
  RtSpecies();

  /*!
    @brief Constructor, constant kappa
  */
  RtSpecies(const std::string a_name, const Real a_kappa);

  /*!
    @brief Constructor, spatially dependent kappa
  */
  RtSpecies(const std::string a_name, Real (*a_kappa)(const RealVect a_pos));

  /*!
    @brief Define function
  */
  virtual void define(const std::string a_name, const Real a_kappa);

  /*!
    @brief Define function
  */
  virtual void define(const std::string a_name, Real (*a_kappa)(const RealVect a_pos));

  /*!
    @brief Get group name
  */
  virtual std::string getName() const;

  /*!
    @brief Constant kappa or not
  */
  virtual bool isKappaConstant() const;

  /*!
    @brief Get kappa at position. You really need to impement this, but it's not pure since this class may be useful for testing. 
  */
  virtual Real getKappa(const RealVect a_pos) const;

  /*!
    @brief Get scattering coefficient. 
  */
  virtual Real getScatteringCoefficient(const RealVect a_pos) const;
	      
protected:

  /*!
    @brief Group name
  */
  std::string m_name;
  
  /*!
    @brief Constant kappa or not
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
  Real (*m_varkappa)(const RealVect a_pos);

  /*!
    @brief Variable kappa
  */
  Real (*m_varscatter)(const RealVect a_pos);
};

#include <CD_NamespaceFooter.H>

#endif
  
