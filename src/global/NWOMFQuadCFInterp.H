/*!
  @file NWOMFQuadCFInterp.H
  @brief Wrapper class for holding multifluid EBQuadCFInterps
  @author Robert Marskar
  @date Dec. 2017
  @todo Segregate implementation file. 
*/

#ifndef _NWOMFQUADCFINTERP_
#define _NWOMFQUADCFINTERP_

#include <RefCountedPtr.H>
#include <NWOEBQuadCFInterp.H>

/*!
  @brief Class for filling ghost cells alà NWOEBQuadCFInterp for multiphase
*/
class NWOMFQuadCFInterp {
public:

  /*!
    @brief Constructor
  */
  NWOMFQuadCFInterp(){

  }

  /*!
    @brief Destructor
  */
  ~NWOMFQuadCFInterp(){

  }

  /*!
    @brief Full constructor. Calls define
  */
  NWOMFQuadCFInterp(const Vector<RefCountedPtr<NWOEBQuadCFInterp> >& a_quadcfi){
    this->define(a_quadcfi);
  }

  /*!
    @brief Define function
  */
  void define(const Vector<RefCountedPtr<NWOEBQuadCFInterp> >& a_quadcfi){
    m_quadcfi = a_quadcfi;
  }

  /*!
    @brief Get interpolation utility
  */
  const RefCountedPtr<NWOEBQuadCFInterp>& get_quadcfi_ptr(const int a_phase) const {
    return m_quadcfi[a_phase];
  }

  /*!
    @brief get interpolation utility for phase
  */
  NWOEBQuadCFInterp& get_quadcfi(const int a_phase) {
    return *m_quadcfi[a_phase];
  }

  /*!
    @brief get interpolation utility for phase
  */
  const NWOEBQuadCFInterp& get_quadcfi(const int a_phase) const {
    return *m_quadcfi[a_phase];
  }
  
protected:


  /*!
    @brief Interpolation for each phase. 
  */
  Vector<RefCountedPtr<NWOEBQuadCFInterp> > m_quadcfi;
  
};

#endif
