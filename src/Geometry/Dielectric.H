/*!
  @file Dielectric.H
  @brief Dielectric struct
  @author Robert marskar
  @date Nov. 2017
*/

#ifndef _DIELECTRIC_
#define _DIELECTRIC_

#include <BaseIF.H>
#include <RefCountedPtr.H>

/*!
  @brief Dielectric class. Simple tuple.
*/
class Dielectric {
public:

  //
  Dielectric();
  
  /*!
    @brief Constructor. Calls define
  */
  Dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  /*!
    @brief 
  */
  ~Dielectric();

  /*!
    @brief Define function
  */
  void define(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  //
  const RefCountedPtr<BaseIF>& get_function() const;

  //
  const Real& get_permittivity() const;

protected:
  std::pair<RefCountedPtr<BaseIF>, Real> m_tuple;
};
#endif
