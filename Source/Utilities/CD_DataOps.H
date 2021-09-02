/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DataOps.H
  @brief  Agglomeration of useful data operations
  @author Robert Marskar
*/

#ifndef CD_DataOps_H
#define CD_DataOps_H

// Std includes

// Our includes
#include <CD_EBAMRData.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Agglomeration of useful data operations. 
*/
class DataOps {
public:
  
  template <typename T>
  static void setValue(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value);

  template <typename T>
  static int sgn(const T a_value);

  static void averageCellVectorToFaceScalar(EBAMRFluxData&               a_facedata,
					    const EBAMRCellData&         a_celldata,
					    const Vector<ProblemDomain>& a_domains);

  static void averageCellVectorToFaceScalar(LevelData<EBFluxFAB>&       a_facedata,
					    const LevelData<EBCellFAB>& a_celldata,
					    const ProblemDomain&        a_domain);

  static void averageCellToFace(EBAMRFluxData& a_face_data,
				const EBAMRCellData& a_cell_data,
				const Vector<ProblemDomain>& a_domains);

  static void averageCellToFace(LevelData<EBFluxFAB>&       a_facedata,
				const LevelData<EBCellFAB>& a_celldata,
				const ProblemDomain&        a_domain);

  static void averageFaceToCell(EBAMRCellData&               a_celldata,
				const EBAMRFluxData&         a_fluxdata,
				const Vector<ProblemDomain>& a_domains);

  static void averageFaceToCell(LevelData<EBCellFAB>&       a_celldata,
				const LevelData<EBFluxFAB>& a_fluxdata,
				const ProblemDomain&        a_domain);

  static void dotProduct(MFAMRCellData& a_result, const MFAMRCellData& a_data1, const MFAMRCellData& a_data2);

  static void dotProduct(LevelData<MFCellFAB>& a_result, const LevelData<MFCellFAB>& a_data1, const LevelData<MFCellFAB>& a_data2);
  
  static void dotProduct(EBAMRCellData& a_result, const EBAMRCellData& a_data1, const EBAMRCellData& a_data2);

  static void dotProduct(LevelData<EBCellFAB>& a_result, const LevelData<EBCellFAB>& a_data1, const LevelData<EBCellFAB>& a_data2);

  static void dotProduct(EBCellFAB& a_result, const EBCellFAB& a_data1, const EBCellFAB& a_data2, const Box& a_box);

  static void incr(MFAMRCellData& a_lhs, const MFAMRCellData& a_rhs, const Real a_scale);
  
  static void incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Real a_scale);
  
  static void incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale);
  
  static void incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real& a_scale);

  static void plus(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_srcComp, const int a_dstComp, const int a_numComp);

  static void plus(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_srcComp, const int a_dstComp, const int a_numComp);

  static void incr(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs, const Real& a_scale);

  static void incr(EBFluxFAB& a_lhs, const EBFluxFAB& a_rhs, const Real& a_scale);
  
  static void incr(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs, const Real& a_scale);

  static void incr(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs, const Real& a_scale);

  static void incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real& a_scale);

  static void incr(EBAMRIFData& a_lhs, const EBAMRIFData& a_rhs, const Real& a_scale);

  static void incr(LevelData<DomainFluxIFFAB>& a_lhs, const LevelData<DomainFluxIFFAB>& a_rhs, const Real& a_scale);

  static void incr(EBAMRCellData& a_lhs, const EBAMRIVData& a_rhs, const Real a_scale);

  static void incr(LevelData<EBCellFAB>& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real a_scale);

  static void incr(EBAMRIVData& a_lhs, const EBAMRCellData& a_rhs, const Real a_scale);

  static void incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale);

  static void copy(MFAMRCellData& a_dst, const MFAMRCellData& a_src);

  static void copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src);

  static void copy(EBAMRIVData& a_dst, const EBAMRIVData& a_src);

  static void divide(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_lcomp, const int a_rcomp);

  static void divide(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_lcomp, const int a_rcomp);

  static void divideByScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);
  
  static void divideByScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void exponentiate(EBAMRCellData& a_lhs, const Real a_factor);

  static void exponentiate(LevelData<EBCellFAB>& a_lhs, const Real a_factor);

  static void floor(EBAMRCellData& a_lhs, const Real a_value);

  static void floor(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  static void floor(EBAMRIVData& a_lhs, const Real a_value);

  static void floor(LevelData<BaseIVFAB<Real> >& a_lhs, const Real a_value); 

  static void getMaxMin(Real& max, Real& min, EBAMRCellData& a_E, const int a_comp);

  static void getMaxMin(Real& max, Real& min, LevelData<EBCellFAB>& a_E, const int a_comp);

  static void getMaxMin(Vector<Real>& a_max, Vector<Real>& a_min, Vector<EBAMRCellData>& a_data);

  static void getMaxMinNorm(Real& a_max, Real& a_min, EBAMRCellData& data);

  static void getMaxMinNorm(Real& a_max, Real& a_min, LevelData<EBCellFAB>& data);

  static void getMaxMinNorm(Real& a_max, Real& a_min, EBAMRIVData& data);

  static void getMaxMinNorm(Real& a_max, Real& a_min, LevelData<BaseIVFAB<Real> >& data);

  static void invert(EBAMRFluxData& a_data);

  static void invert(LevelData<EBFluxFAB>& a_data);

  static void invert(EBFluxFAB& a_data);

  static void kappaSum(Real& a_mass, const LevelData<EBCellFAB>& a_lhs);

  static void kappaScale(EBAMRCellData& a_data);

  static void kappaScale(LevelData<EBCellFAB>& a_data);

  static void kappaScale(MFAMRCellData& a_data);

  static void kappaScale(LevelData<MFCellFAB>& a_data);

  static void laplacian(EBAMRCellData& a_lapl, const EBAMRCellData& a_data);
  
  static void laplacian(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data);

  static void genLaplacian(EBAMRCellData& a_lapl, const EBAMRCellData& a_data);
  
  static void genLaplacian(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data);

  static void flashError(EBAMRCellData& a_lapl, const EBAMRCellData& a_data, const Real a_eps);
  
  static void flashError(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data, const Real a_eps);

  static void ln(EBAMRCellData& a_lhs);

  static void ln(LevelData<EBCellFAB>& a_lhs);

  static void multiply(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  static void multiply(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void multiply(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs);

  static void multiply(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs);

  static void multiply(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs);

  static void multiply(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs);
  
  static void multiplyScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  static void multiplyScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void multiplyScalar(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs);
  
  static void multiplyScalar(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs);
  
  static void norm(Real& a_norm, const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p);

  static void scale(MFAMRCellData& a_lhs, const Real& a_scale);

  static void scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale);

  static void scale(EBAMRIVData& a_lhs, const Real& a_scale);

  static void scale(EBAMRCellData& a_lhs, const Real a_scale);

  static void scale(LevelData<EBCellFAB>& a_lhs, const Real a_scale);

  static void scale(EBAMRFluxData& a_lhs, const Real a_scale);

  static void scale(LevelData<BaseIVFAB<Real> >& a_lhs, const Real& a_scale);

  static void setCoveredValue(EBAMRCellData& a_lhs, const int a_comp, const Real a_value);

  static void setCoveredValue(LevelData<EBCellFAB>& a_lhs, const int a_comp, const Real a_value);

  static void setValue(MFAMRCellData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp);

  static void setValue(LevelData<MFCellFAB>& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp);

  static void setValue(EBAMRCellData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp);

  static void setValue(LevelData<EBCellFAB>& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp);

  static void setValue(EBAMRFluxData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp);

  static void setValue(LevelData<EBFluxFAB>& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp);

  static void setValue(EBAMRIVData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp);

  static void setValue(LevelData<BaseIVFAB<Real> >& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp);

  static void setValue(EBAMRCellData& a_lhs, const std::function<RealVect(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx);

  static void setValue(LevelData<EBCellFAB>& a_lhs, const std::function<RealVect(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx);

  static void setValue(EBAMRCellData& a_lhs, const Real& a_value);
  
  static void setValue(EBAMRCellData& a_lhs, const Real a_value, const int a_comp);

  static void setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value, const int a_comp);

  static void setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  static void setValue(LevelData<EBFluxFAB>& a_lhs, const Real a_value);
  
  static void setValue(LevelData<BaseIVFAB<Real> >& a_lhs, const Real a_value);

  static void setValue(EBAMRFluxData& a_lhs, const Real& a_value);

  static void setValue(EBAMRIVData& a_lhs, const Real& a_value);

  static void setValue(MFAMRCellData& a_lhs, const Real& a_value);

  static void setValue(LevelData<MFCellFAB>& a_lhs, const Real& a_value);

  static void setValue(MFAMRFluxData& a_lhs, const Real& a_value);

  static void setValue(LevelData<MFFluxFAB>& a_lhs, const Real& a_value);

  static void setValue(MFAMRIVData& a_lhs, const Real& a_value);

  static void setValue(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value);

  static void setValue(EBAMRIFData& a_lhs, const Real a_value);

  static void setValue(LevelData<DomainFluxIFFAB>& a_lhs, const Real a_value);

  static void sum(Real& a_value);

  static void squareRoot(EBAMRFluxData& a_lhs);

  static void squareRoot(LevelData<EBFluxFAB>& a_lhs);

  static void vectorLength(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  static void vectorLength(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void vectorLength(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& box);

  static void vectorLength2(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  static void vectorLength2(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  static void vectorLength2(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& box);

  static void computeMinValidBox(RealVect& a_lo, RealVect& a_hi, const RealVect a_normal, const RealVect a_centroid);
  
  static bool allCornersInsideEb(const Vector<RealVect>& a_corners, const RealVect a_normal, const RealVect a_centroid);

  static void shiftCorners(Vector<RealVect>& a_corners, const RealVect& a_distance);

  static void filterSmooth(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_stride, const Real a_alpha);
  
  static void filterSmooth(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_stride, const Real a_alpha);

  static void computeParticleWeights(unsigned long long&      a_weight,
				     unsigned long long&      a_num,
				     unsigned long long&      a_remainder,
				     const unsigned long long a_numPhysicalParticles,
				     const int                a_ppc);

};

#include <CD_NamespaceFooter.H>

#include <CD_DataOpsImplem.H>

#endif
