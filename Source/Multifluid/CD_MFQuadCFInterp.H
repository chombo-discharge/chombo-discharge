/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFQuadCFInterp.H
  @brief  Wrapper class for holding multifluid EBQuadCFInterps
  @author Robert Marskar
*/

#ifndef CD_MFQuadCFInterp_H
#define CD_MFQuadCFInterp_H

// Chombo includes
#include <RefCountedPtr.H>
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
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
    @brief get interpolation utility for phase
  */
  const RefCountedPtr<EBQuadCFInterp>& getEBQuadCFInterpPointer(const int a_phase) const;

  /*!
    @brief get interpolation utility for phase
  */
  EBQuadCFInterp& getEBQuadCFInterp(const int a_phase);
  
protected:

  /*!
    @brief Interpolation for each phase. 
  */
  Vector<RefCountedPtr<EBQuadCFInterp> > m_quadcfi;
};

#include <CD_NamespaceFooter.H>

#endif
