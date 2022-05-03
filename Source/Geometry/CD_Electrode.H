/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Electrode.H
  @brief  Declaration of an electrode class
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
  @brief Electrode class for encapsulating the concept of an electrical conductor. 
*/
class Electrode
{
public:
  /*!
    @brief Default constructor. Must subsequently call define()
  */
  Electrode();

  /*!
    @brief Constructor. 
    @param[in] a_baseIF          Implicit function for the electrode
    @param[in] a_live            Live electrode or not
    @param[in] a_voltageFraction Fraction of live voltage (if live)
    @note Calls the define function
  */
  Electrode(const RefCountedPtr<BaseIF>& a_baseIF, const bool a_live, const Real a_voltageFraction = 1.0);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~Electrode();

  /*!
    @brief Define function
    @param[in] a_baseIF          Implicit function for the electrode
    @param[in] a_live            Live electrode or not
    @param[in] a_voltageFraction Fraction of live voltage (if live)
  */
  virtual void
  define(const RefCountedPtr<BaseIF>& a_baseIF, const bool a_live, const Real a_voltageFraction = 1.0);

  /*!
    @brief Get implicit function
    @return Returns m_baseIF
  */
  virtual const RefCountedPtr<BaseIF>&
  getImplicitFunction() const;

  /*!
    @brief Get live or not. 
    @return Returns m_isLive
  */
  virtual const bool&
  isLive() const;

  /*!
    @brief Get fraction of live voltage. 
    @return m_voltageFraction
  */
  virtual const Real&
  getFraction() const;

protected:
  /*!
    @brief Implicit function
  */
  RefCountedPtr<BaseIF> m_baseIF;

  /*!
    @brief Live or not
  */
  bool m_isLive;

  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Fraction of the live potential
  */
  Real m_voltageFraction;
};

#include <CD_NamespaceFooter.H>

#endif
