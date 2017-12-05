/*!
  @file MFQuadCFInterp.H
  @brief Wrapper class for holding multifluid EBQuadCFInterps
  @author Robert Marskar
  @date Dec. 2017
  @todo Segregate implementation file. 
*/

#ifndef _MFQUADCFINTERP_
#define _MFQUADCFINTERP_

#include <RefCountedPtr.H>
#include <EBQuadCFInterp.H>

/*!
  @brief Class for filling ghost cells alà EBQuadCFInterp for multiphase
*/
class MFQuadCFInterp {
public:

  /*!
    @brief Constructor
  */
  MFQuadCFInterp(){

  }

  /*!
    @brief Destructor
  */
  ~MFQuadCFInterp(){

  }

  /*!
    @brief Full constructor. Calls define
  */
  MFQuadCFInterp(const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadcfi){
    this->define(a_quadcfi);
  }

  /*!
    @brief Define function
  */
  void define(const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadcfi){
    m_quadcfi = a_quadcfi;
  }

  /*!
    @brief get interpolation utility for phase
  */
  EBQuadCFInterp& get_quadcfi(const int a_phase) {
    return *m_quadcfi[a_phase];
  }

  /*!
    @brief get interpolation utility for phase
  */
  const EBQuadCFInterp& get_quadcfi(const int a_phase) const {
    return *m_quadcfi[a_phase];
  }
  
protected:


  /*!
    @brief Interpolation for each phase. 
  */
  Vector<RefCountedPtr<EBQuadCFInterp> > m_quadcfi;
  
};

#endif
