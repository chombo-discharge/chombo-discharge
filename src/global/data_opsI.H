/*!
  @file data_opsI.H
  @brief Implementation of templated things in data_ops.H
*/

#ifndef _DATA_OPSI_
#define _DATA_OPSI_

#include "CD_NamespaceHeader.H"

template <typename T> void data_ops::set_value(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(a_value);
  }
}

template <typename T> int data_ops::sgn(const T a_value){
  return ((a_value > 0.0) - (a_value < 0.0));
}
#include "CD_NamespaceFooter.H"

#endif
