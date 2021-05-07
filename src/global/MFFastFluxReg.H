/*!
  @file   MFFastFluxReg.H
  @brief  Wrapper class for holding multifluid EBFastFluxReg
  @author Robert Marskar
  @date   Aug. 2020
  @todo   Segregate implementation file. 
*/

#ifndef _MFFastFluxReg_H_
#define _MFFastFluxReg_H_

#include <RefCountedPtr.H>
#include <EBFluxRegister.H>

namespace ChomboDischarge {
  
  /*!
    @brief Class for filling ghost cells alà EBFasterFR for multiphase
  */
  class MFFastFluxReg {
  public:

    /*!
      @brief Constructor
    */
    MFFastFluxReg(){

    }

    /*!
      @brief Destructor
    */
    ~MFFastFluxReg(){

    }

    /*!
      @brief Full constructor. Calls define
    */
    MFFastFluxReg(const Vector<RefCountedPtr<EBFluxRegister> >& a_fastfr){
      this->define(a_fastfr);
    }

    /*!
      @brief Define function
    */
    void define(const Vector<RefCountedPtr<EBFluxRegister> >& a_fastfr){
      m_fastfr = a_fastfr;
    }

    /*!
      @brief Get interpolation utility
    */
    const RefCountedPtr<EBFluxRegister>& get_fastfr_ptr(const int a_phase) const {
      return m_fastfr[a_phase];
    }

    /*!
      @brief get interpolation utility for phase
    */
    EBFluxRegister& get_fastfr(const int a_phase) {
      return *m_fastfr[a_phase];
    }

    /*!
      @brief get interpolation utility for phase
    */
    const EBFluxRegister& get_fastfr(const int a_phase) const {
      return *m_fastfr[a_phase];
    }
  
  protected:


    /*!
      @brief Interpolation for each phase. 
    */
    Vector<RefCountedPtr<EBFluxRegister> > m_fastfr;
  
  };
}

#endif
