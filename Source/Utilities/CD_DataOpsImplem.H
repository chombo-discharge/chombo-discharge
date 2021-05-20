/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file  CD_DataOpsImplem.H
  @brief Implementation of CD_DataOps.H
*/

#ifndef CD_DataOpsImplem_H
#define CD_DataOpsImplem_H

#include <CD_NamespaceHeader.H>

template <typename T> void DataOps::setValue(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(a_value);
  }
}

template <typename T> int DataOps::sgn(const T a_value){
  return ((a_value > 0.0) - (a_value < 0.0));
}

#include <CD_NamespaceFooter.H>

#endif
