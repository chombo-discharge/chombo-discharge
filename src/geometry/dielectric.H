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

/*!
  @brief dielectric class. Simple tuple.
*/
class dielectric {
public:

  //
  dielectric();
  
  /*!
    @brief Constructor. Calls define
  */
  dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  /*!
    @brief 
  */
  ~dielectric();

  /*!
    @brief Define function
  */
  void define(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  /*!
    @brief Get level-set function
  */
  const RefCountedPtr<BaseIF>& get_function() const;

  /*!
    @brief Get permittivity
  */
  const Real& get_permittivity() const;

protected:
  std::pair<RefCountedPtr<BaseIF>, Real> m_tuple;
};
#endif
