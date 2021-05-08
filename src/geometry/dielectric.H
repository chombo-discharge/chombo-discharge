/*!
  @file dielectric.H
  @brief dielectric struct
  @author Robert marskar
  @date Nov. 2017
*/

#ifndef _DIELECTRIC_
#define _DIELECTRIC_

#include <BaseIF.H>
#include <RefCountedPtr.H>

#include "CD_NamespaceHeader.H"
/*!
  @brief dielectric class. Simple tuple.
*/
class dielectric {
public:

  /*!
    @brief Weak constructor
  */
  dielectric();
  
  /*!
    @brief Constructor. Calls define
  */
  dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  /*!
    @brief Constructor, allows for spatially dependent permittivity
  */
  dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos));

  /*!
    @brief 
  */
  virtual ~dielectric();

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
  virtual const RefCountedPtr<BaseIF>& get_function() const;

  /*!
    @brief Get permittivity
  */
  virtual Real get_permittivity(const RealVect a_pos) const;

protected:

  RefCountedPtr<BaseIF> m_baseif;

  Real m_permittivity;

  Real (*m_variablepermittivity)(const RealVect a_pos);
  
  bool m_constant;
};
#include "CD_NamespaceFooter.H"
#endif
