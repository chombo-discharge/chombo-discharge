/*!
  @file   CD_DataOpsImplem.H
  @brief  Implementation of CD_DataOps.H
  @author Robert Marskar
  @date   Nov 2017
*/

#ifndef _CD_DataOpsImplem_H_
#define _CD_DataOpsImplem_H_

template <typename T> void DataOps::setValue(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(a_value);
  }
}

template <typename T> int DataOps::sgn(const T a_value){
  return ((a_value > 0.0) - (a_value < 0.0));
}

#endif
