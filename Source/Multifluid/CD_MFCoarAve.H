/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFCoarAve.H
  @brief  Wrapper class for holding multifluid EBCoarAves
  @author Robert Marskar
*/

#ifndef CD_MFCoarAve_H
#define CD_MFCoarAve_H

// Std includes
#include <map>

// Our includes
#include <CD_EBCoarAve.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for coarsening data in a multifluid context. 
*/
class MFCoarAve
{
public:
  /*!
    @brief Constructor. Must subsequently call defined. 
  */
  MFCoarAve();

  /*!
    @brief Destructor (does nothing). 
  */
  ~MFCoarAve();

  /*!
    @brief Full constructor. Calls define
    @param[in] aveOps Averaging operators on each phase. 
  */
  MFCoarAve(const Vector<RefCountedPtr<EBCoarAve>>& a_aveOps);

  /*!
    @brief Define function
    @param[in] aveOps Averaging operators on each phase. 
  */
  void
  define(const Vector<RefCountedPtr<EBCoarAve>>& a_aveOps);

  /*!
    @brief Get coarsening utility for specified phase.
    @param[in] a_phase Phase
  */
  const RefCountedPtr<EBCoarAve>&
  getAveOp(const int a_phase) const;

protected:
  /*!
    @brief Interpolation for each phase. 
  */
  Vector<RefCountedPtr<EBCoarAve>> m_aveOps;
};

#include <CD_NamespaceFooter.H>

#endif
