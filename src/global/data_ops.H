/*!
  @file data_ops.H
  @brief Agglomeration of useful data operations
  @author Robert Marskar
  @date Nov 2017
*/

#ifndef _data_ops_
#define _data_ops_

#include "type_definitions.H"
#include "MFInterfaceFAB.H"

#include "EBLevelDataOps.H"

/*!
  @brief Agglomeration of useful data operations. Doing this because I use smart pointer all around
*/
class data_ops {
public:
  static void incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale);

  static void scale(EBAMRIVData& a_lhs, const Real& a_scale);

  static void scale(LevelData<BaseIVFAB<Real> >& a_lhs, const Real& a_scale);

  static void set_value(EBAMRCellData& a_lhs, const Real& a_value);

  static void set_value(EBAMRFluxData& a_lhs, const Real& a_value);

  static void set_value(EBAMRIVData& a_lhs, const Real& a_value);

  static void set_value(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value);

  template <typename T> static void set_value(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value);

};

#include "data_opsI.H"

#endif
