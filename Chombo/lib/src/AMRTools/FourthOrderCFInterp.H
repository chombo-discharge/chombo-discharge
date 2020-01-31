#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERCFINTERP_H_
#define _FOURTHORDERCFINTERP_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "FourthOrderFineInterp.H"

#include "NamespaceHeader.H"

/// Fourth-order interpolation in time and space to ghost cells.

/**
 */
class FourthOrderCFInterp
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  FourthOrderCFInterp()
  {
    m_defined = false;
  }

  /// Full constructor
  FourthOrderCFInterp(/// layout at this level
                       const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                       /// layout at coarser level
                       const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                       /// number of variables
                       const int&                a_numStates,
                       /// problem domain on the coarser level
                       const ProblemDomain&      a_coarseDomain,
                       /// refinement ratio between this level and the coarser level
                       const int&                a_refineCoarse,
                       /// number of layers of ghost cells to fill by interpolation
                       const int&                a_interpRadius)
  {
    define(a_thisDisjointBoxLayout, a_coarserDisjointBoxLayout,
           a_numStates, a_coarseDomain, a_refineCoarse, a_interpRadius);
  }

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~FourthOrderCFInterp()
  {
  }

  ///fill ghostcells as if all the coarse data were zero (useful for multigrid)
  void 
  homogeneousCoarseFineInterp(/// interpolated solution on this level
                              LevelData<FArrayBox>&         a_fineData,
                              /// solution on coarser level
                              int                           a_srcComp,
                              /// starting fine data component
                              int                           a_dstComp,
                              /// number of data components to interpolate
                              int                           a_numComp);
  /// Actual constructor.
  /**
     Set up object.
   */
  void define(/// layout at this level
              const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              /// layout at coarser level
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              /// number of variables
              const int&                a_numStates,
              /// problem domain on the coarser level
              const ProblemDomain&      a_coarseDomain,
              /// refinement ratio between this level and the coarser level
              const int&                a_refineCoarse,
              /// number of layers of ghost cells to fill by interpolation
              const int&                a_interpRadius);



  
  /// Interpolate in space only.
  /**
     At a fixed time, interpolate in space to ghost cells of a_fine
     from a_coarse.
   */
  void coarseFineInterp(/// interpolated solution on this level
                        LevelData<FArrayBox>&         a_fineData,
                        /// solution on coarser level
                        const LevelData<FArrayBox>&   a_coarseData,
                        /// starting coarse data component
                        int                           a_srcComp,
                        /// starting fine data component
                        int                           a_dstComp,
                        /// number of data components to interpolate
                        int                           a_numComp);



protected:

  /// whether define() has been called
  bool m_defined;

  /// box layout for this level
  DisjointBoxLayout m_layout;

  /// this level's layout coarsened
  DisjointBoxLayout m_layoutCoarsened;

  /// box layout for the coarse level
  DisjointBoxLayout m_coarseLayout;

  /// number of layers of fine ghost cells to fill by interpolation
  int m_interpRadius;

  /// problem domain at the coarser level
  ProblemDomain m_coarseDomain;

  /// refinement ratio between this level and the next coarser
  int m_refineCoarse;

  /// 1 in m_fixedDims, m_refineCoarse in other dimensions
  IntVect m_refineVect;

  /// number of variables
  int m_numStates;

  /// data on ghosted coarsened fine grids at intermediate time in fillInterp
  LevelData<FArrayBox> m_coarsenedFineData;

  /// coarsened ghost cells of fine patches
  LayoutData<IntVectSet> m_coarsenedGhosts;


  /// interpolator in space
  FourthOrderFineInterp m_spaceInterpolator;

  /// Interpolate in space only, from coarsened fine grids.
  void fillInterpSpaceFromCoarsened(/// interpolated solution on this level
                                    LevelData<FArrayBox>&         a_fineData,
                                    /// starting coarse data component
                                    int                           a_srcComp,
                                    /// starting fine data component
                                    int                           a_dstComp,
                                    /// number of data components to interpolate
                                    int                           a_numComp);

private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderCFInterp&);
  FourthOrderCFInterp(const FourthOrderCFInterp&);
};

#include "NamespaceFooter.H"
#endif
