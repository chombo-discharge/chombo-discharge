/*!
  @file data_ops.cpp
  @brief Implementation of data_ops.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "data_ops.H"


void data_ops::incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    EBLevelDataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void data_ops::set_value(EBAMRCellData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}
