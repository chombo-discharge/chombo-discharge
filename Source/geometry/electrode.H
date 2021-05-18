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

#include "CD_NamespaceHeader.H"
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
  electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

  /*!
    @brief 
  */
  virtual ~electrode();

  /*!
    @brief Define function
  */
  virtual void define(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

  /*!
    @brief Get implicit function
  */
  virtual const RefCountedPtr<BaseIF>& get_function() const;

  /*!
    @brief Get implicit function
  */
  virtual const bool& is_live() const;

  /*!
    @brief Get fraction
  */
  virtual const Real& get_fraction() const;

protected:
  std::pair<RefCountedPtr<BaseIF>, bool> m_tuple;

  /*!
    @brief Fraction of the live potential
  */
  Real m_fraction;
};
#include "CD_NamespaceFooter.H"

#endif
