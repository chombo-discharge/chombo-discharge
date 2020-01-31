#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFLEVELDATAOPS_H_
#define _MFLEVELDATAOPS_H_

#include "IntVect.H"
#include "LevelData.H"
#include "RefCountedPtr.H"

#include "MFCellFAB.H"

#include "EBCellFAB.H"
// #include "EBCellFactory.H"
#include "EBISLayout.H"
#include "NamespaceHeader.H"

///
/**
   Repository of useful MF functions over a level.  All functions are static so just
   use MFLevelDataOps::blah(...);
 */
class MFLevelDataOps
{
public:
  MFLevelDataOps()
  {
  }

  ~MFLevelDataOps()
  {
  }

  ///
  /**
   */
  static void setIrregVal(LevelData<MFCellFAB>&     a_data,
                          const DisjointBoxLayout&  a_dbl,
                          const Vector<EBISLayout>& a_ebisl,
                          const Real&               a_value);

  ///
  /**
   */
  static void setCoveredVal(LevelData<MFCellFAB>&   a_lData,
                            const Real&             a_value);


  ///
  /**
   */
  static void setCoveredVal(LevelData<MFCellFAB>&   a_lData,
                            const int&              a_comp,
                            const Real&             a_value);

  ///
  /**
   */
  static void averageMultiVofsToRegFAB(LevelData<MFCellFAB>&      a_data,
                                       const DisjointBoxLayout&   a_dbl,
                                       const Vector<EBISLayout>&  a_ebisl);

  ///
  /**
   */
  static void copyToMultiVofsFromRegFAB(LevelData<MFCellFAB>&     a_data,
                                        const DisjointBoxLayout&  a_dbl,
                                        const Vector<EBISLayout>& a_ebisl);

  ///
  /**
   */
  static void defineLevelData(LevelData<MFCellFAB>&      a_levelData,
                              const Vector<EBISLayout>&  a_ebisl,
                              const DisjointBoxLayout&   a_dbl,
                              const IntVect&             a_ghosts,
                              const int&                 a_nComp);

  ///
  /**
   */
  static void setToZero(LevelData<MFCellFAB>& a_result);


  ///
  /**
   */
  static void scale(LevelData<MFCellFAB>& a_lhs,
                    const Real&           a_scale) ;


  ///
  /**
   */
  static void scale(LevelData<MFCellFAB>& a_lhs,
                    const Real& a_scale,
                    const int&  a_comp);


  ///
  /**
   */
  static void setVal(LevelData<MFCellFAB>& a_result,
                     const Real&           a_value);

  ///
  /**
   */
  static void setVal(LevelData<MFCellFAB>& a_result,
                     const Real&           a_value,
                     const int&            a_comp);

  ///
  /**
   */
  static void assign(LevelData<MFCellFAB>&       a_to,
                     const LevelData<MFCellFAB>& a_from,
                     const Interval&             a_toInterva,
                     const Interval&             a_fromInterval);

  ///
  /**
   */
  static void assign(LevelData<MFCellFAB>& a_lhs,
                     const LevelData<MFCellFAB>& a_rhs);


  ///
  /**
   */
  static void clone(LevelData<MFCellFAB>& a_lhs,
                    const LevelData<MFCellFAB>& a_rhs);


  ///
  /**
   */
  static void incr( LevelData<MFCellFAB>& a_lhs,
                    const LevelData<MFCellFAB>&   a_rhs,
                    const Real& a_scale);


  ///
  /**
   */
  static void incr( LevelData<MFCellFAB>& a_lhs,
                    const Real& a_scale);


  ///
  /**
   */
  static void axby( LevelData<MFCellFAB>&       a_lhs,
                    const LevelData<MFCellFAB>& a_x,
                    const LevelData<MFCellFAB>& a_y,
                    const Real& a,
                    const Real& b);


  ///
  /**
   */
  static void axby( LevelData<MFCellFAB>&       a_lhs,
                    const LevelData<MFCellFAB>& a_x,
                    const LevelData<MFCellFAB>& a_y,
                    const Real& a,
                    const Real& b,
                    const int&  a_lhsComp,
                    const int&  a_xComp,
                    const int&  a_yComp);



  ///
  /**
   */
  static void sum(LevelData<MFCellFAB>&       a_result,
                  const LevelData<MFCellFAB>& a_in1,
                  const LevelData<MFCellFAB>& a_in2);


  ///
  /**
   */
  static void addConstant(LevelData<MFCellFAB>& a_data,
                          const Real&           a_constant);


  ///
  /**
   */
  static void power(LevelData<MFCellFAB>& a_result,
                    const Real&           a_exponent,
                    const int&            a_comp);

  ///
  /**
   */
  static void product(LevelData<MFCellFAB>&       a_result,
                      const LevelData<MFCellFAB>& a_in1,
                      const LevelData<MFCellFAB>& a_in2);


  ///
  /**
   */
  static void product(LevelData<MFCellFAB>&       a_result,
                      const LevelData<MFCellFAB>& a_in1,
                      const LevelData<MFCellFAB>& a_in2,
                      const int&                  a_rComp,
                      const int&                  a_1Comp,
                      const int&                  a_2Comp);


  ///
  /**
   */
  static void invert(LevelData<MFCellFAB>&       a_result,
                     const LevelData<MFCellFAB>& a_in1);


  ///
  /**
   */
  static  void divideVectorByScalar(LevelData<MFCellFAB>&       a_vectorOut,
                                    const LevelData<MFCellFAB>& a_vectorIn,
                                    const LevelData<MFCellFAB>& a_scalar);

  ///
  /**
   */
  static void vectorMagnitude(LevelData<MFCellFAB>&       a_scalarOut,
                              const LevelData<MFCellFAB>& a_vectorIn,
                              const int&                  a_pval);

  ///
  /**
   */
  static  void divide(LevelData<MFCellFAB>&       a_result,
                      const LevelData<MFCellFAB>& a_in1,
                      const LevelData<MFCellFAB>& a_in2);


  ///
  /**
   */
  static  void divide(LevelData<MFCellFAB>&       a_result,
                      const LevelData<MFCellFAB>& a_in1,
                      const LevelData<MFCellFAB>& a_in2,
                      const int&                  a_rComp,
                      const int&                  a_1Comp,
                      const int&                  a_2Comp);


  ///
  /**
   */
  static  void kappaWeight(LevelData<MFCellFAB>& a_data);


  ///
  /**
   */
  static  void kappaWeight(MFCellFAB& a_data);


  ///
  /**
   */
  static  void kappaScale(LevelData<MFCellFAB>& a_data,
                          const Real&           a_scale);


private:

};
#include "NamespaceFooter.H"
#endif
