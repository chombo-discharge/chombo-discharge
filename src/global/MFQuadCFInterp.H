/*!
  @file   MFQuadCFInterp.H
  @brief  Wrapper class for holding multifluid EBQuadCFInterps
  @author Robert Marskar
  @date   Dec. 2017
*/

#ifndef _MFQUADCFINTERP_
#define _MFQUADCFINTERP_

#include <RefCountedPtr.H>
#include <EBQuadCFInterp.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Class for filling ghost cells alà EBQuadCFInterp for multiphase
*/
class MFQuadCFInterp {
public:

  /*!
    @brief Constructor
  */
  MFQuadCFInterp();

  /*!
    @brief Destructor
  */
  ~MFQuadCFInterp();

  /*!
    @brief Full constructor. Calls define
  */
  MFQuadCFInterp(const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadcfi);

  /*!
    @brief Define function
  */
  void define(const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadcfi);

  /*!
    @brief Get interpolation utility
  */
  const RefCountedPtr<EBQuadCFInterp>& getNWOEBQuadCFInterp_ptr(const int a_phase) const;

  /*!
    @brief get interpolation utility for phase
  */
  EBQuadCFInterp& getNWOEBQuadCFInterp(const int a_phase);

  /*!
    @brief get interpolation utility for phase
  */
  const EBQuadCFInterp& getNWOEBQuadCFInterp(const int a_phase) const;
  
protected:


  /*!
    @brief Interpolation for each phase. 
  */
  Vector<RefCountedPtr<EBQuadCFInterp> > m_quadcfi;
  
};
#include "CD_NamespaceFooter.H"

#endif
