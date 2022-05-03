/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFMultigridInterpolator.H
  @brief  Declaration of a wrapper class for holding multifluid EBMultigridInterpolators
  @author Robert Marskar
*/

#ifndef CD_MFMultigridInterpolator_H
#define CD_MFMultigridInterpolator_H

// Std includes
#include <map>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Wrapper class for holding multifluid EBMultigridInterpolators. 
*/
class MFMultigridInterpolator
{
public:
  /*!
    @brief Weak constructor. Must subsquently call define.
  */
  MFMultigridInterpolator();

  /*!
    @brief Defining constructor.
    @param[in] a_interpolators Interpolators on each phase. 
  */
  MFMultigridInterpolator(const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators);

  /*!
    @brief Destructor. Does nothing.
  */
  ~MFMultigridInterpolator();

  /*!
    @brief Copy assignment
    @param[in] a_other Other interpolator
  */
  MFMultigridInterpolator&
  operator=(const MFMultigridInterpolator& a_other);

  /*!
    @brief Define function
    @param[in] a_interpolators Interpolators on each phase. 
  */
  void
  define(const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators);

  /*!
    @brief Get interpolator for specified phase
    @param[in] a_phase
  */
  RefCountedPtr<EBMultigridInterpolator>&
  getInterpolator(const int a_phase) const;

  /*!
    @brief Get number of filled ghost cells
  */
  int
  getGhostCF() const;

protected:
  /*!
    @brief Interpolators on each phase. 
  */
  mutable Vector<RefCountedPtr<EBMultigridInterpolator>> m_interpolators;
};

#include <CD_NamespaceFooter.H>

#endif
