/*!
  @file   MFFastFR.H
  @brief  Wrapper class for holding multifluid EBFastFR's
  @author Robert Marskar
  @date   Aug. 2020
  @todo   Segregate implementation file. 
*/

#ifndef _MFFASTFR_H_
#define _MFFASTFR_H_

#include <RefCountedPtr.H>
#include <MFFastFR.H>

/*!
  @brief Class for filling ghost cells alà EBFastFR for multiphase
*/
class MFFastFR {
public:

  /*!
    @brief Constructor
  */
  MFFastFR(){

  }

  /*!
    @brief Destructor
  */
  ~MFFastFR(){

  }

  /*!
    @brief Full constructor. Calls define
  */
  MFFastFR(const Vector<RefCountedPtr<EBFastFR> >& a_fastfr){
    this->define(a_fastfr);
  }

  /*!
    @brief Define function
  */
  void define(const Vector<RefCountedPtr<EBFastFR> >& a_fastfr){
    m_fastfr = a_fastfr;
  }

  /*!
    @brief Get interpolation utility
  */
  const RefCountedPtr<EBFastFR>& get_fastfr_ptr(const int a_phase) const {
    return m_fastfr[a_phase];
  }

  /*!
    @brief get interpolation utility for phase
  */
  EBFastFR& get_fastfr(const int a_phase) {
    return *m_fastfr[a_phase];
  }

  /*!
    @brief get interpolation utility for phase
  */
  const EBFastFR& get_fastfr(const int a_phase) const {
    return *m_fastfr[a_phase];
  }
  
protected:


  /*!
    @brief Interpolation for each phase. 
  */
  Vector<RefCountedPtr<EBFastFR> > m_fastfr;
  
};

#endif
