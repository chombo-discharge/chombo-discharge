/*!
  @file Electrode.H
  @brief Electrode struct
  @author Robert marskar
  @date Nov. 2017
*/

#ifndef _ELECTRODE_
#define _ELECTRODE_

#include <BaseIF.H>
#include <RefCountedPtr.H>

/*!
  @brief Electrode class. Simple tuple.
*/
class Electrode {
public:

  //
  Electrode();
  
  /*!
    @brief Constructor. Calls define
  */
  Electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live);

  /*!
    @brief 
  */
  ~Electrode();

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
