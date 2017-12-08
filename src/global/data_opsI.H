/*!
  @file data_opsI.H
  @brief Implementation of templated things in data_ops.H
*/

#ifndef _DATA_OPSI_
#define _DATA_OPSI_

template <typename T> void data_ops::set_value(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(a_value);
  }
}

#endif
