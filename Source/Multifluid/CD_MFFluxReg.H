/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFFluxReg.H
  @brief  Declaration of a class for refluxing in a multiphase context. 
  @author Robert Marskar
*/

#ifndef CD_MFFluxReg_H
#define CD_MFFluxReg_H

// Chombo includes
#include <RefCountedPtr.H>
#include <EBFluxRegister.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Class which wraps EBFluxRegisters in multiphase. 
*/
class MFFluxReg {
public:

  /*!
    @brief Weak constructor. Must subsequently call define. 
  */
  MFFluxReg();
  
  /*!
    @brief Full constructor. Calls define
    @param[in] a_fluxRegs Flux registers on each phase
  */
  MFFluxReg(const Vector<RefCountedPtr<EBFluxRegister> >& a_fluxRegs);

  /*!
    @brief Destructor (does nothing)
  */
  ~MFFluxReg();

  /*!
    @brief Define function
    @param[in] a_fluxRegs Flux registers on each phase
  */
  void define(const Vector<RefCountedPtr<EBFluxRegister> >& a_fluxRegs);

  /*!
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
  */
  const RefCountedPtr<EBFluxRegister>& getFluxRegPointer(const int a_phase) const;

  /*!
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
  */
  EBFluxRegister& getFluxReg(const int a_phase);

  /*!
    @brief Get flux register for specified phase
    @param[in] a_phase Phase
  */
  const EBFluxRegister& getFluxReg(const int a_phase) const;
  
protected:

  /*!
    @brief Flux register for each phase. 
  */
  Vector<RefCountedPtr<EBFluxRegister> > m_fluxRegs;
  
};

#include <CD_NamespaceFooter.H>

#endif
