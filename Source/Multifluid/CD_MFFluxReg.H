/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFFluxReg.H
  @brief  Wrapper class for holding multifluid EBFastFluxReg
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
  @brief Class for filling ghost cells alà EBFasterFR for multiphase
*/
class MFFluxReg {
public:

  /*!
    @brief Constructor
  */
  MFFluxReg();
  
  /*!
    @brief Full constructor. Calls define
  */
  MFFluxReg(const Vector<RefCountedPtr<EBFluxRegister> >& a_fastfr);

  /*!
    @brief Destructor
  */
  ~MFFluxReg();

  /*!
    @brief Define function
  */
  void define(const Vector<RefCountedPtr<EBFluxRegister> >& a_fastfr);

  /*!
    @brief Get interpolation utility
  */
  const RefCountedPtr<EBFluxRegister>& getFluxRegPointer(const int a_phase) const;

  /*!
    @brief get interpolation utility for phase
  */
  EBFluxRegister& getFluxReg(const int a_phase);

  /*!
    @brief get interpolation utility for phase
  */
  const EBFluxRegister& getFluxReg(const int a_phase) const;
  
protected:

  /*!
    @brief Interpolation for each phase. 
  */
  Vector<RefCountedPtr<EBFluxRegister> > m_fastfr;
  
};

#include <CD_NamespaceFooter.H>

#endif
