/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_Electrode.H
  @brief  Declaring electrod
  @author Robert marskar
*/

#ifndef CD_Electrode_H
#define CD_Electrode_H

// Chombo includes
#include <BaseIF.H>
#include <RefCountedPtr.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief electrode class. Simple tuple.
*/
class Electrode {
public:

  //
  Electrode();
  
  /*!
    @brief Constructor. Calls define
  */
  Electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

  /*!
    @brief 
  */
  virtual ~Electrode();

  /*!
    @brief Define function
  */
  virtual void define(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

  /*!
    @brief Get implicit function
  */
  virtual const RefCountedPtr<BaseIF>& getImplicitFunction() const;

  /*!
    @brief Get implicit function
  */
  virtual const bool& isLive() const;

  /*!
    @brief Get fraction
  */
  virtual const Real& getFraction() const;

protected:
  std::pair<RefCountedPtr<BaseIF>, bool> m_tuple;

  /*!
    @brief Fraction of the live potential
  */
  Real m_fraction;
};

#include <CD_NamespaceFooter.H>

#endif
