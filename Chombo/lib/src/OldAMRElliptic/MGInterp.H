#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MGINTERP_H_
#define _MGINTERP_H_

#include "REAL.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

class DisjointBoxLayout;

/// Internal Chombo class to interpolate coarse level data to a fine level
/**
    Class to interpolate coarse level data to a fine level
    data for multigrid.  This class should be considered internal
    to AMRSolver and not a part of the Chombo API.
*/
class MGInterp
{
public:
  /// Null constructor
  MGInterp();

  /// Destructor
  ~MGInterp();

  /// Full constructor (calls full define)
  MGInterp(const DisjointBoxLayout& a_fineDomain,
           int                      a_numcomps,
           int                      a_refRatio,
           const Box&               a_problemDomain);

  /// Full constructor (calls full define)
  MGInterp(const DisjointBoxLayout& a_fineDomain,
           int                      a_numcomps,
           int                      a_refRatio,
           const ProblemDomain&     a_problemDomain);

  /// Full define
  /**
      Specify the fine level domain and refinement ratio.
  */
  void define(const DisjointBoxLayout& a_fineDomain,     // the fine level domain
              int                      a_numcomps,       // the number of components
              int                      a_refRatio,       // the refinement ratio
              const Box&               a_problemDomain); // problem domain

  /// Full define
  /**
      Specify the fine level domain and refinement ratio.
  */
  void define(const DisjointBoxLayout& a_fineDomain,     // the fine level domain
              int                      a_numcomps,       // the number of components
              int                      a_refRatio,       // the refinement ratio
              const ProblemDomain&     a_problemDomain); // problem domain

  ///
  /**
      Interpolate coarse level data to fine level
  */
  void interpToFine(LevelData<FArrayBox>&       a_fineData,
                    const LevelData<FArrayBox>& a_coarseData);

  /// Is the object defined
  bool isDefined() const;

protected:
  bool m_isDefined;

  // the refinement ratio
  int m_refRatio;

  // work array for the coarse level data in a domain that is the
  // outline of the fine level domain on the coarse level
  //
  LevelData<FArrayBox> m_coarsenedFineData;
  DisjointBoxLayout m_grids;

  // problem domain
  ProblemDomain m_problemDomain;
};

#include "NamespaceFooter.H"
#endif
