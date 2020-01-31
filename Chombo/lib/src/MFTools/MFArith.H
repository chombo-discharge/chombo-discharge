#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFARITH_H_
#define _MFARITH_H_

#include "EBArith.H"

#include "MFCellFAB.H"
#include "EBISLayout.H"
#include "NamespaceHeader.H"
///
/**
   class to encapsulate the common arithmetic operations
   for mffabs
 */
class MFArith
{
public:
  ///
  /**
   */
  static int orderScript(int icomp, int inorm, int ncomp);
  ///
  /**
   */
  static void
  compareError(Vector<Real>&                            a_orders,
               const Vector< LevelData<MFCellFAB>* >&   a_errorFine,
               const Vector< LevelData<MFCellFAB>* >&   a_errorCoar,
               const Vector< DisjointBoxLayout >&       a_gridsFine,
               const Vector< DisjointBoxLayout >&       a_gridsCoar,
               const Vector< Vector<EBISLayout> >&      a_ebislvFine,
               const Vector< Vector<EBISLayout> >&      a_ebislvCoar,
               const Vector<int>&                       a_refRat,
               const Box&                               a_coarseDom,
               int a_testverbosity);

  ///
  /** return l-p norm of a_src.
      if p==0, v
      norm returned is Max(abs(a_src)) over
      uncovered regions.  otherwise,
      returns 1/vol(sum(volfrac*a_src(iv,comp)**p)^(1/p))
      of component comp of a_src
      weighted by local volume fraction and
      also returns volume of uncovered regions.
      Only uncovered regions count here.
      The data must have the same layout as a_layout with
      the possible exception of ghost cells.

      int pmode = -2;  //norm = (1/v)(sum(phi dv)) ---no absolute values and multiply kappa as you go
      int pmode = -1;  //norm = (1/v)(sum(phi dv)) ---no absolute values and assume kappa already multiplied in

   */
  static
  Real
  norm(const Vector< LevelData<MFCellFAB>* >& a_data,
       const Vector< DisjointBoxLayout > &    a_layout,
       const Vector< Vector<EBISLayout> >&    a_ebislv,
       const Vector< int >&                   a_refRatio,
       const int&                             a_comp,
       const int&                             a_p,
       EBNormType::NormMode = EBNormType::OverBoth);

  static
  Real
  normBoundaryExcluded(const Vector< LevelData<MFCellFAB>* >& a_src,
                       const Vector< DisjointBoxLayout >&     a_grids,
                       const Vector< Vector<EBISLayout> >&    a_ebislv,
                       const ProblemDomain&                   a_coarDomain,
                       const Vector<int>&                     a_refRatio,
                       const int&                             a_coarBoundaryExclusion,
                       const int&                             a_comp,
                       const int&                             a_pval,
                       EBNormType::NormMode = EBNormType::OverBoth);

  static
  void
  normByPhase(Vector<Real>&                          a_phaseNorms,
              const Vector< LevelData<MFCellFAB>* >& a_data,
              const Vector< DisjointBoxLayout > &    a_layout,
              const Vector< Vector<EBISLayout> >&    a_ebislv,
              const Vector< int >&                   a_refRatio,
              const int&                             a_comp,
              const int&                             a_p,
              EBNormType::NormMode = EBNormType::OverBoth);

};

#include "NamespaceFooter.H"
#endif
