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

/*!
  @brief Agglomeration of useful data operations. Doing this because I use smart pointer all around
*/
class data_ops {
public:
  template <typename T> static void set_value(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value);
  
  static void incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale);

  static void copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src);

  static void scale(MFAMRCellData& a_lhs, const Real& a_scale);

  static void scale(EBAMRIVData& a_lhs, const Real& a_scale);

  static void scale(EBAMRCellData& a_lhs, const Real a_scale);

  static void scale(EBAMRFluxData& a_lhs, const Real a_scale);

  static void scale(LevelData<BaseIVFAB<Real> >& a_lhs, const Real& a_scale);

  static void set_value(EBAMRCellData& a_lhs, const Real& a_value);

  static void set_value(EBAMRFluxData& a_lhs, const Real& a_value);

  static void set_value(EBAMRIVData& a_lhs, const Real& a_value);

  static void set_value(MFAMRCellData& a_lhs, const Real& a_value);

  static void set_value(LevelData<MFCellFAB>& a_lhs, const Real& a_value);

  static void set_value(MFAMRFluxData& a_lhs, const Real& a_value);

  static void set_value(LevelData<MFFluxFAB>& a_lhs, const Real& a_value);

  static void set_value(MFAMRIVData& a_lhs, const Real& a_value);

  static void set_value(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value);

};

#include "data_opsI.H"

#endif
