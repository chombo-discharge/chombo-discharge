#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Fri, July 16, 1999

#ifndef _LEVELMG_H_
#define _LEVELMG_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "CoarseAverage.H"
#include "LevelOp.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Internal Chombo class to manage multigrid solve on a level
/**
    Multigrid solver on a single multigrid level.
    This class should be considered internal
    to AMRSolver and not a part of the Chombo API.
 */
class LevelMG
{
public:
  /// Has this object been fully defined?
  bool isDefined() const;

  /// Null Constructor
  LevelMG();

  /// Deprecated full constructor (use ProblemDomain instead of Box)
  LevelMG(const DisjointBoxLayout& a_ba,
          const DisjointBoxLayout* a_baseBaPtr,
          Real                     a_dxLevel,
          int                      a_refRatio,
          const  Box&              a_domain,
          int                      a_nCoarserLevels,
          const LevelOp* const     a_opin,
          int                      a_ncomp=1);

  /// Full constructor -- calls matching define function
  LevelMG(const DisjointBoxLayout& a_ba,
          const DisjointBoxLayout* a_baseBaPtr,
          Real                     a_dxLevel,
          int                      a_refRatio,
          const ProblemDomain&     a_domain,
          int                      a_nCoarserLevels,
          const LevelOp* const     a_opin,
          int                      a_ncomp=1);

  /// Deprecated full define function (use ProblemDomain instead of Box)
  void define(const DisjointBoxLayout& a_ba,
              const DisjointBoxLayout* a_baseBaPtr,
              Real                     a_dxLevel,
              int                      a_refRatio,
              const Box&               a_domain,
              int                      a_nCoarserLevels,
              const LevelOp* const     a_opin,
              int                      a_ncomp=1);

  /// Full Define Function
  void define(/// Grids at this multigrid level
              const DisjointBoxLayout& a_ba,
              /// Coarser grids
              const DisjointBoxLayout* a_baseBaPtr,
              /// Cell spacing at this level
              Real                     a_dxLevel,
              /// Refinement ratio
              int                      a_refRatio,
              /// Physical domain of problem
              const ProblemDomain&     a_domain,
              /// Number of multigrid levels coarser than this one
              int                      a_nCoarserLevels,
              /// LevelOp which encapsulates operator
              const LevelOp* const     a_opin,
              /// number of components
              int                      a_ncomp);

  /// Constructor for coarsened version of object.
  void define(const LevelMG&       a_level,
              int                  a_refCoarse,
              const LevelOp* const a_opin);

  /// Constructor for coarsened version of object.
  LevelMG(const LevelMG& a_level,
          int            a_refCoarse,
          const LevelOp* a_opin);

  /// Destructor
  ~LevelMG();

  /// Invoke multigrid relaxation step.
  /** Default is pure MG V-cycle, suitable
      for use in multilevel solver application; otherwise, use approximate
      solver such as CG at bottom level. It is assumed that the problem has
      already been put in residual-correction form. In particular, only
      the homogeneous form of the physical and coarse-fine boundary
      conditions need be invoked.
  */
  void mgRelax(LevelData<FArrayBox> & a_soln,
               const LevelData<FArrayBox> & a_rhs,
               bool bottomsolveflag);

  /// Set number of smoothing steps at bottom
  void setNumBottomGSRB(int a_numBottomGSRB)
  {
    m_numBottomGSRB = a_numBottomGSRB;
  }

  /// Set number of smoothing steps on the way up V-cycle
  void setNumSmoothUp(int a_numSmoothUp);

  /// Set number of smoothing steps on the way down V-cycle
  void setNumSmoothDown(int a_numSmoothDown);

  /// this is a dangerous access function that should not generally be used.
  LevelOp* levelOpPtr();

  /// this is another access function that should not be generally used.
  LevelMG* lCoarsePtr();

  /// sets convergence metric -- essentially a pass through to LevelOp
  void setConvergenceMetric(Real a_metric, int a_comp);

protected:
  void setDefaultValues();
  void clearMemory();

  bool m_isDefined;

  //
  LevelData<FArrayBox>  m_resid;

  //
  LevelData<FArrayBox>  m_crseResid;

  //
  LevelData<FArrayBox>  m_crseCorr;

  //
  DisjointBoxLayout m_ba;

  //
  DisjointBoxLayout m_baCoarsened;

  //
  const DisjointBoxLayout* m_baseBaPtr;

  //
  Real m_dxLevel;

  //
  ProblemDomain m_domain;

  //
  int m_nCoarserLevels;

  //
  LevelMG* m_lCoarsePtr;

  //
  LevelOp* m_levelopPtr;

  //
  int m_numBottomGSRB;

  //
  int m_numSmoothUp;

  //
  int m_numSmoothDown;

  //
  int m_refToCoar;

  //
  CoarseAverage m_averageOp;

private:
  /**
     correct fine on intersection with crse
     should only be called internally because this
     is not written for general grid layouts
  */
  void crseCorrect(LevelData<FArrayBox>&       a_fine,
                   const LevelData<FArrayBox>& a_crse,
                   int                         a_refRat);

  void operator=(const LevelMG& a_levmgin)
  {
  }

  LevelMG(const LevelMG& a_levmgin)
  {
  }
};

#include "NamespaceFooter.H"
#endif
