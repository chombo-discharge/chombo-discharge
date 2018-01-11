/*!
  @file electrode.H
  @brief electrode struct
  @author Robert marskar
  @date Nov. 2017
*/

#ifndef _ELECTRODE_
#define _ELECTRODE_

#include <BaseIF.H>
#include <RefCountedPtr.H>

/*!
  @brief electrode class. Simple tuple.
*/
class electrode {
public:

  //
  electrode();
  
  /*!
    @brief Constructor. Calls define
  */
  electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live);

  /*!
    @brief 
  */
  ~electrode();

  /*!
    @brief Define function
  */
  void define(RefCountedPtr<BaseIF> a_baseif, bool a_live);

  //
  const RefCountedPtr<BaseIF>& get_function() const;

  //
  const bool& is_live() const;

protected:
  std::pair<RefCountedPtr<BaseIF>, bool> m_tuple;
};
#endif
