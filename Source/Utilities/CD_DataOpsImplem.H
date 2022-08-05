/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DataOpsImplem.H
  @brief  Implementation of CD_DataOps.H
  @author Robert Marskar
*/

#ifndef CD_DataOpsImplem_H
#define CD_DataOpsImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_NamespaceHeader.H>

template <typename T>
void
DataOps::setValue(LevelData<MFInterfaceFAB<T>>& a_lhs, const T& a_value)
{
  CH_TIME("DataOps::setValue");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].setVal(a_value);
  }
}

template <typename T>
int
DataOps::sgn(const T a_value)
{
  return ((a_value > 0.0) - (a_value < 0.0));
}

template <typename T>
void
DataOps::copy(EBAMRData<T>& a_dst, const EBAMRData<T>& a_src, const Interval a_dstInterv, const Interval a_srcInterv) {
  CH_TIME("DataOps::copy(EBAMRData)");

  CH_assert(a_dst.size() == a_src.size());
  CH_assert(a_dstInterv.size() == a_srcInterv.size());

  for (int lvl = 0; lvl < a_dst.size(); lvl++) {
    CH_assert(a_dst[lvl]->nComp() > a_dstInterv.end());
    CH_assert(a_src[lvl]->nComp() > a_srcInterv.end());    

    a_src[lvl]->localCopyTo(a_srcInterv, *a_dst[lvl], a_dstInterv);
  }
}

#include <CD_NamespaceFooter.H>

#endif
