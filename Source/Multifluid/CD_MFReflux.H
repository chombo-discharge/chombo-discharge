/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFReflux.H
  @brief  Declaration of a class for refluxing in a multiphase context. 
  @author Robert Marskar
*/

#ifndef CD_MFReflux_H
#define CD_MFReflux_H

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_EBReflux.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which wraps EBRefluxs in multiphase. 
*/
class MFReflux
{
public:
  /*!
    @brief Weak constructor. Must subsequently call define. 
  */
  MFReflux();

  /*!
    @brief Full constructor. Calls define
    @param[in] a_fluxRegs Flux registers on each phase
  */
  MFReflux(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs);

  /*!
    @brief Destructor (does nothing)
  */
  ~MFReflux();

  /*!
    @brief Define function
    @param[in] a_fluxRegs Flux registers on each phase
  */
  void
  define(const Vector<RefCountedPtr<EBReflux>>& a_fluxRegs);

  /*!
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
  */
  const RefCountedPtr<EBReflux>&
  getFluxRegPointer(const int a_phase) const;

  /*!
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
  */
  EBReflux&
  getFluxReg(const int a_phase);

  /*!
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
  */
  const EBReflux&
  getFluxReg(const int a_phase) const;

protected:
  /*!
    @brief Flux register for each phase. 
  */
  Vector<RefCountedPtr<EBReflux>> m_fluxRegs;
};

#include <CD_NamespaceFooter.H>

#endif
