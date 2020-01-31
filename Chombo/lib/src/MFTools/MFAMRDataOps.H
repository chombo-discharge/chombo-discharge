#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFAMRDATAOPS_H_
#define _MFAMRDATAOPS_H_

#include "MFLevelDataOps.H"
#include "MFIndexSpace.H"

#include "FaceIterator.H"
// #include "EBCoarseAverage.H"
// #include "EBCoarsen.H"

#include "EBIndexSpace.H"
#include "EBAMRIO.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

///
/**
   Repository of functions for MF Hierarchies.
 */
class MFAMRDataOps
{
public:
  MFAMRDataOps()
  {
  }

  ~MFAMRDataOps()
  {
  }

  ///
  /**
   */
  static void setCoveredAMRVal(Vector<LevelData<MFCellFAB>* >&      a_data,
                               const Vector< Vector<EBISLayout> >&  a_ebislv,
                               const Vector<int>&                   a_refRat,
                               const Real&                          a_value);

  static void setCoveredAMRVal(Vector<LevelData<MFCellFAB>* >&       a_data,
                               const Vector< Vector<EBLevelGrid> >&  a_eblg,
                               const Vector<int>&                    a_refRat,
                               const Real&                           a_value);

  ///
  /**
   */
  static void setCoveredVal(Vector<LevelData<MFCellFAB>* >& a_data,
                            const Real&                     a_value);

  ///
  /**
   */
  static void scale(Vector<LevelData<MFCellFAB>* >& a_lhs,
                    const Real&                     a_scale) ;


  ///
  /**
   */
  static void scale(Vector<LevelData<MFCellFAB>* >& a_lhs,
                    const Real&                     a_scale,
                    const int&                      a_comp);


  ///
  /**
   */
  static void setToZero(Vector<LevelData<MFCellFAB>* >& a_result);


  ///
  /**
   */
  static void setVal(Vector<LevelData<MFCellFAB>* >& a_result,
                     const Real&                     a_value);


  ///
  /**
   */
  static void setVal(Vector<LevelData<MFCellFAB>* >& a_result,
                     const Real&                     a_value,
                     const int&                      a_comp);

  ///
  /**
   */
  static void assign(Vector<LevelData<MFCellFAB>* >&       a_lhs,
                     const Vector<LevelData<MFCellFAB>* >& a_rhs);

  ///
  /**
   */
  static void assign(Vector<RefCountedPtr<LevelData<MFCellFAB> > >&       a_lhs,
                     const Vector<LevelData<MFCellFAB>* >&                a_rhs);

  ///
  /**
   */
  static void assign(Vector<LevelData<MFCellFAB>* >&       a_to,
                     const Vector<LevelData<MFCellFAB>* >& a_from,
                     const Interval&                       a_toInterval,
                     const Interval&                       a_fromInterval);

  ///
  /**
   */
  static void assign(Vector<RefCountedPtr<LevelData<MFCellFAB> > >& a_to,
                     const Vector<LevelData<MFCellFAB>* >&          a_from,
                     const Interval&                                a_toInterval,
                     const Interval&                                a_fromInterval);

  ///
  /**
   */
  static void incr(Vector<LevelData<MFCellFAB>* >& a_lhs,
                   const Vector<LevelData<MFCellFAB>* >&   a_rhs,
                   const Real& a_scale);

  ///
  /**
   */
  static void incr(Vector<LevelData<MFCellFAB>* >& a_lhs,
                   const Real& a_scale);

  ///
  /**
   */
  static void axby(Vector<LevelData<MFCellFAB>* >&       a_lhs,
                   const Vector<LevelData<MFCellFAB>* >& a_x,
                   const Vector<LevelData<MFCellFAB>* >& a_y,
                   const Real& a,
                   const Real& b);


  ///
  /**
   */
  static void axby(Vector<LevelData<MFCellFAB>* >&       a_lhs,
                   const Vector<LevelData<MFCellFAB>* >& a_x,
                   const Vector<LevelData<MFCellFAB>* >& a_y,
                   const Real& a_a,
                   const Real& a_b,
                   const int&  a_lhsComp,
                   const int&  a_xComp,
                   const int&  a_yComp);

  ///
  /**
   */
  static  void sum(Vector<LevelData<MFCellFAB>* >&       a_result,
                   const Vector<LevelData<MFCellFAB>* >& a_in1,
                   const Vector<LevelData<MFCellFAB>* >& a_in2);

  ///
  /**
   */
  static  void addConstant(Vector<LevelData<MFCellFAB>* >& a_data,
                           const Real&                     a_constant);

  ///
  /**
   */
  static  void product(Vector<LevelData<MFCellFAB>* >&       a_result,
                       const Vector<LevelData<MFCellFAB>* >& a_in1,
                       const Vector<LevelData<MFCellFAB>* >& a_in2);


  ///
  /**
   */
  static  void product(Vector<LevelData<MFCellFAB>* >&       a_result,
                       const Vector<LevelData<MFCellFAB>* >& a_in1,
                       const Vector<LevelData<MFCellFAB>* >& a_in2,
                       const int&                            a_rComp,
                       const int&                            a_1Comp,
                       const int&                            a_2Comp);

  ///
  /**
   */
  static void divideVectorByScalar(Vector<LevelData<MFCellFAB>* >&       a_vectorOut,
                                   const Vector<LevelData<MFCellFAB>* >& a_vectorIn,
                                   const Vector<LevelData<MFCellFAB>* >& a_scalar);

  ///
  /**
   */
  static  void divide(Vector<LevelData<MFCellFAB>* >&       a_result,
                      const Vector<LevelData<MFCellFAB>* >& a_in1,
                      const Vector<LevelData<MFCellFAB>* >& a_in2);


  ///
  /**
   */
  static  void divide(Vector<LevelData<MFCellFAB>* >&       a_result,
                      const Vector<LevelData<MFCellFAB>* >& a_in1,
                      const Vector<LevelData<MFCellFAB>* >& a_in2,
                      const int&                            a_rComp,
                      const int&                            a_1Comp,
                      const int&                            a_2Comp);

  ///
  /**
   */
  static  void kappaWeight(Vector<LevelData<MFCellFAB>* >& a_data);


  ///
  /**
   */
  static  void kappaScale(Vector<LevelData<MFCellFAB>* >& a_data,
                          const Real&                     a_scale);

  ///
  /**
   */
  static Real subtractOffMean(Vector<LevelData<MFCellFAB>* >&     a_data,
                              const Vector<DisjointBoxLayout>&    a_grids,
                              const Vector< Vector<EBISLayout> >& a_ebislv,
                              const Vector<int>&                  a_refRat);

  static Real subtractOffMean(Vector<LevelData<MFCellFAB>* >&      a_data,
                              const Vector< Vector<EBLevelGrid> >& a_eblg,
                              const Vector<int>&                   a_refRat);


protected:

};
#include "NamespaceFooter.H"
#endif
