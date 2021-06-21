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
#include <CD_Phases.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Wrapper class for holding multifluid EBMultigridInterpolators. 
*/
class MFMultigridInterpolator {
public:

  /*!
    @brief Weak constructor. Must subsquently call define.
  */
  MFMultigridInterpolator();

  /*!
    @brief Disallowed copy constructor. 
  */
  MFMultigridInterpolator(const MFMultigridInterpolator& a_other) = default;
  
  /*!
    @brief Full constructor. Calls define function
  */
  MFMultigridInterpolator(const std::map<Phase, RefCountedPtr<EBMultigridInterpolator> >& a_interpolators);

  /*!
    @brief Destructor. Does nothing.
  */
  ~MFMultigridInterpolator();

  /*!
    @brief Disallowed copy assignment.
  */
  MFMultigridInterpolator& operator=(const MFMultigridInterpolator& a_other) = default;

  /*!
    @brief Define function
  */
  void define(const std::map<Phase, RefCountedPtr<EBMultigridInterpolator> >& a_interpolators);

  /*!
    @brief Get interpolator for phase
  */
  const RefCountedPtr<EBMultigridInterpolator>& getInterpolator(const Phase a_phase) const;

protected:

  /*!
    @brief Interpolators on each phase. 
  */
  std::map<Phase, RefCountedPtr<EBMultigridInterpolator> > m_interpolators;
};

#include <CD_NamespaceFooter.H>

#endif
