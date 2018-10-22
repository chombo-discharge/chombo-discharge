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

  template <typename T> static int sgn(const T a_value);

  static void average_cell_to_face(EBAMRFluxData&               a_facedata,
				   const EBAMRCellData&         a_celldata,
				   const Vector<ProblemDomain>& a_domains);

  static void average_cell_to_face(LevelData<EBFluxFAB>&       a_facedata,
				   const LevelData<EBCellFAB>& a_celldata,
				   const ProblemDomain&        a_domain);

  static void average_cell_to_face_allcomps(EBAMRFluxData& a_face_data,
					    const EBAMRCellData& a_cell_data,
					    const Vector<ProblemDomain>& a_domains);

  static void average_cell_to_face_allcomps(LevelData<EBFluxFAB>&       a_facedata,
					    const LevelData<EBCellFAB>& a_celldata,
					    const ProblemDomain&        a_domain);

  static void dot_prod(EBAMRCellData& a_result, const EBAMRCellData& a_data1, const EBAMRCellData& a_data2);

  static void dot_prod(LevelData<EBCellFAB>& a_result, const LevelData<EBCellFAB>& a_data1, const LevelData<EBCellFAB>& a_data2);

  static void incr(MFAMRCellData& a_lhs, const MFAMRCellData& a_rhs, const Real a_scale);
  
  static void incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Real a_scale);
  
  static void incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale);
  
  static void incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real& a_scale);

  static void incr(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs, const Real& a_scale);

  static void incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real& a_scale);

  static void incr(EBAMRCellData& a_lhs, const EBAMRIVData& a_rhs, const Real a_scale);

  static void incr(LevelData<EBCellFAB>& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real a_scale);

  static void incr(EBAMRIVData& a_lhs, const EBAMRCellData& a_rhs, const Real a_scale);

  static void incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale);

  static void copy(MFAMRCellData& a_dst, const MFAMRCellData& a_src);

  static void copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src);

  static void copy(EBAMRIVData& a_dst, const EBAMRIVData& a_src);

  static void divide(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_lcomp, const int a_rcomp);

  static void divide(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_lcomp, const int a_rcomp);

  static void divide_scalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);
  
  static void divide_scalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void exponentiate(EBAMRCellData& a_lhs, const Real a_factor);

  static void exponentiate(LevelData<EBCellFAB>& a_lhs, const Real a_factor);

  static void floor(EBAMRCellData& a_lhs, const Real a_value);

  static void floor(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  static void get_max_min(Real& max, Real& min, EBAMRCellData& a_E, const int a_comp);

  static void get_max_min(Real& max, Real& min, LevelData<EBCellFAB>& a_E, const int a_comp);

  static void get_max_min(Vector<Real>& a_max, Vector<Real>& a_min, Vector<EBAMRCellData>& a_data);

  static void get_max_min_norm(Real& a_max, Real& a_min, EBAMRCellData& data);

  static void get_max_min_norm(Real& a_max, Real& a_min, LevelData<EBCellFAB>& data);

  static void get_max_min_norm(Real& a_max, Real& a_min, EBAMRIVData& data);

  static void get_max_min_norm(Real& a_max, Real& a_min, LevelData<BaseIVFAB<Real> >& data);

  static void kappa_sum(Real& a_mass, const LevelData<EBCellFAB>& a_lhs);

  static void kappa_scale(EBAMRCellData& a_data);

  static void kappa_scale(LevelData<EBCellFAB>& a_data);

  static void kappa_scale(MFAMRCellData& a_data);

  static void kappa_scale(LevelData<MFCellFAB>& a_data);

  static void ln(EBAMRCellData& a_lhs);

  static void ln(LevelData<EBCellFAB>& a_lhs);

  static void multiply(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  static void multiply(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void multiply_scalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void multiply_scalar(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs);
  
  static void multiply_scalar(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs);
  
  static void norm(Real& a_norm, const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p);

  static void scale(MFAMRCellData& a_lhs, const Real& a_scale);

  static void scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale);

  static void scale(EBAMRIVData& a_lhs, const Real& a_scale);

  static void scale(EBAMRCellData& a_lhs, const Real a_scale);

  static void scale(LevelData<EBCellFAB>& a_lhs, const Real a_scale);

  static void scale(EBAMRFluxData& a_lhs, const Real a_scale);

  static void scale(LevelData<BaseIVFAB<Real> >& a_lhs, const Real& a_scale);

  static void set_covered_value(EBAMRCellData& a_lhs, const int a_comp, const Real a_value);

  static void set_covered_value(LevelData<EBCellFAB>& a_lhs, const int a_comp, const Real a_value);

  static void set_value(EBAMRCellData& a_lhs, const Real& a_value);
  
  static void set_value(EBAMRCellData& a_lhs, const Real a_value, const int a_comp);

  static void set_value(LevelData<EBCellFAB>& a_lhs, const Real a_value, const int a_comp);

  static void set_value(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  static void set_value(LevelData<EBFluxFAB>& a_lhs, const Real a_value);
  
  static void set_value(LevelData<BaseIVFAB<Real> >& a_lhs, const Real a_value);

  static void set_value(LevelData<BaseIFFAB<Real> >& a_lhs, const Real a_value);

  static void set_value(EBAMRFluxData& a_lhs, const Real& a_value);

  static void set_value(EBAMRIVData& a_lhs, const Real& a_value);

  static void set_value(EBAMRIFData& a_lhs, const Real& a_value);

  static void set_value(MFAMRCellData& a_lhs, const Real& a_value);

  static void set_value(LevelData<MFCellFAB>& a_lhs, const Real& a_value);

  static void set_value(MFAMRFluxData& a_lhs, const Real& a_value);

  static void set_value(LevelData<MFFluxFAB>& a_lhs, const Real& a_value);

  static void set_value(MFAMRIVData& a_lhs, const Real& a_value);

  static void set_value(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value);

  static void sum(Real& a_value);

  static void vector_length(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  static void vector_length(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void vector_length(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& box);

};

#include "data_opsI.H"

#endif
