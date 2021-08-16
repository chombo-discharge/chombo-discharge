/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Dielectric.H
  @brief  Dielectric class
  @author Robert marskar
*/

#ifndef CD_Dielectric_H
#define CD_Dielectric_H

// Std includes
#include <functional>

// Chombo includes
#include <BaseIF.H>
#include <RefCountedPtr.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief dielectric class. Simple tuple.
*/
class Dielectric {
public:

  /*!
    @brief Weak constructor, must subsequently call one of the define methods.
  */
  Dielectric();
  
  /*!
    @brief Full constructor which uses constant permittivity.
    @param[in] a_baseIF       Implicit function
    @param[in] a_permittivity Constant permittivity
    @note Calls the define function for constant permittivity.
  */
  Dielectric(const RefCountedPtr<BaseIF>& a_baseIF, const Real a_permittivity);

  /*!
    @brief Full constructor which uses variable permittivity.
    @param[in] a_baseIF       Implicit function
    @param[in] a_permittivity Variable permittivity
    @note Calls the define function for variable permittivity.
  */
  Dielectric(const RefCountedPtr<BaseIF>& a_baseIF, const std::function<Real(const RealVect a_pos)>& a_permittivity);

  /*!
    @brief Destructor (does nothing).
  */
  virtual ~Dielectric();

  /*!
    @brief Define function
    @param[in] a_baseIF       Implicit function
    @param[in] a_permittivity Constant permittivity
  */
  virtual void define(const RefCountedPtr<BaseIF>& a_baseIF, const Real a_permittivity);

  /*!
    @brief Define function for variable permittivity
    @param[in] a_baseIF       Implicit function
    @param[in] a_permittivity Variable permittivity
  */
  virtual void define(const RefCountedPtr<BaseIF>& a_baseIF, const std::function<Real(const RealVect a_pos)>& a_permittivity);

  /*!
    @brief Get level-set function
    @return Returns m_baseIF
  */
  virtual const RefCountedPtr<BaseIF>& getImplicitFunction() const;

  /*!
    @brief Get permittivity
    @param[in] a_pos Position
    @return Returns permittivity for physical coordinates a_pos. 
  */
  virtual Real getPermittivity(const RealVect a_pos) const;

protected:

  /*!
    @brief Implicit function
  */
  RefCountedPtr<BaseIF> m_baseIF;

  /*!
    @brief Constant permittivity
  */
  Real m_constantPermittivity;

  /*!
    @brief Variable permittivity
  */
  std::function<Real(const RealVect a_pos) > m_variablePermittivity;

  /*!
    @brief Use constant permittivity or not. 
  */
  bool m_useConstant;

  /*!
    @brief Is defined or not.
  */
  bool m_isDefined;
};

#include <CD_NamespaceFooter.H>

#endif
