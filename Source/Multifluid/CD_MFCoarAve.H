/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFCoarAve.H
  @brief  Wrapper class for holding multifluid EbCoarAves
  @author Robert Marskar
*/

#ifndef CD_MFCoarAve_H
#define CD_MFCoarAve_H

// Std includes
#include <map>

// Our includes
#include <CD_Phases.H>
#include <CD_EbCoarAve.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Class for filling ghost cells alà EbCoarAve for multiphase
*/
class MFCoarAve {
public:

  /*!
    @brief Constructor
  */
  MFCoarAve();

  /*!
    @brief Destructor
  */
  ~MFCoarAve();

  /*!
    @brief Full constructor. Calls define
  */
  MFCoarAve(const std::map<Phase, RefCountedPtr<EbCoarAve> >& a_aveOps);

  /*!
    @brief Define function
  */
  void define(const std::map<Phase, RefCountedPtr<EbCoarAve> >& a_aveOps);

  /*!
    @brief get interpolation utility for phase
  */
  EbCoarAve& getAveOp(const Phase a_phase);
  
protected:

  /*!
    @brief Interpolation for each phase. 
  */
  std::map<Phase, RefCountedPtr<EbCoarAve> > m_aveOps;
  
};

#include <CD_NamespaceFooter.H>

#endif
