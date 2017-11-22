/*!
  @file data_ops.H
  @brief Agglomeration of useful data operations
  @author Robert Marskar
  @date Nov 2017
*/

#ifndef _data_ops_
#define _data_ops_

#include "type_definitions.H"

#include "EBLevelDataOps.H"

/*!
  @brief Agglomeration of useful data operations. Doing this because I use smart pointer all around
*/
class data_ops {
public:
  static void incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale);

  static void set_value(EBAMRCellData& a_lhs, const Real& a_value);
};
#endif
