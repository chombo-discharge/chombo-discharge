/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFQuadCFInterp.H
  @brief  Wrapper class for holding multifluid EBQuadCFInterps
  @author Robert Marskar
  @oto    Purge old code when multifluid has been reworked. 
*/

#ifndef CD_MFQuadCFInterp_H
#define CD_MFQuadCFInterp_H

// Chombo includes
#include <RefCountedPtr.H>
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_Phases.H>
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
    @brief Full constructor
    @param[in] Interpolators at various phases.
  */
  MFQuadCFInterp(const std::map<Phase, RefCountedPtr<EBQuadCFInterp> >& a_interpolators);

  /*!
    @brief Define function
    @param[in] a_interpolators Interpolators for various phases.
  */
  void define(const std::map<Phase, RefCountedPtr<EBQuadCFInterp> >& a_interpolators);

  /*!
    @brief Get interpolation utility for phase
  */
  const RefCountedPtr<EBQuadCFInterp>& getInterpolator(const Phase a_phase) const;

  // Old stuff below here. 

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

  /*!
    @brief Interpolation utility for each phase. 
  */
  std::map<Phase, RefCountedPtr<EBQuadCFInterp> > m_interpolators;
};

#include <CD_NamespaceFooter.H>

#endif
