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
    @brief Weak constructor
  */
  Dielectric();
  
  /*!
    @brief Constructor. Calls define
  */
  Dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  /*!
    @brief Constructor, allows for spatially dependent permittivity
  */
  Dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos));

  /*!
    @brief 
  */
  virtual ~Dielectric();

  /*!
    @brief Define function
  */
  virtual void define(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  /*!
    @brief Define function for variable permittivity
  */
  virtual void define(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos));

  /*!
    @brief Get level-set function
  */
  virtual const RefCountedPtr<BaseIF>& getImplicitFunction() const;

  /*!
    @brief Get permittivity
  */
  virtual Real getPermittivity(const RealVect a_pos) const;

protected:

  RefCountedPtr<BaseIF> m_baseif;

  Real m_permittivity;

  Real (*m_variablepermittivity)(const RealVect a_pos);
  
  bool m_constant;
};

#include <CD_NamespaceFooter.H>

#endif
