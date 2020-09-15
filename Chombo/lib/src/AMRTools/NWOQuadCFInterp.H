#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NWOQUADCFINTERP_H_
#define _NWOQUADCFINTERP_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "RealVect.H"
#include "ProblemDomain.H"
#include "CFIVS.H"
#include "NamespaceHeader.H"

/// Fourth-order interpolation in time and space to ghost cells.

/**
 */
class NWOQuadCFInterp
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  NWOQuadCFInterp()
  {
    m_defined = false;
  }

  /// Full constructor
  NWOQuadCFInterp(/// layout at this level
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
  ~NWOQuadCFInterp()
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
  DisjointBoxLayout m_layoutCoarsened;

  /// data on ghosted coarsened fine grids at intermediate time in fillInterp
  LevelData<FArrayBox> m_coarsenedFineData;

  /// coarsened ghost cells of fine patches
  LayoutData<CFIVS> m_cfivs;


  void 
  interpOnPatch(FArrayBox&         a_fineData,
                const FArrayBox&   a_coarseData,
                const DataIndex&   a_dit,
                int                a_srcComp,
                int                a_dstComp,
                int                a_numComp);

  void
  extrapolateValue(Real           & a_fineValue, 
                   const Real     & a_coarValue, 
                   const RealVect & a_firstDerivs, 
                   const RealVect & a_secondDerivs, 
                   const RealVect & a_mixedDerivs, 
                   const RealVect & a_distance);

  void
  getDerivs(RealVect& firstDerivs, 
            RealVect& secondDerivs,
            RealVect& mixedDerivs, 
            const FArrayBox & a_data, 
            const IntVect& a_ivCoar,
            const Real   & a_dx,
            const int    & a_icomp);


  int getMixedIndex(int a_derivDir1, int a_derivDir2)
  {
    int retval;
    //0 == xy, 1 = xz, 2 == yz
    if(((a_derivDir1 == 0) && (a_derivDir2 == 1)) || ((a_derivDir1 == 1) && (a_derivDir2 == 0)))
      {
        retval = 0;
      }
    else if (((a_derivDir1 == 0) && (a_derivDir2 == 2)) || ((a_derivDir1 == 2) && (a_derivDir2 == 0)))
      {
        retval = 1;
      }
    else if (((a_derivDir1 == 1) && (a_derivDir2 == 2)) || ((a_derivDir1 == 2) && (a_derivDir2 == 1)))
      {
        retval = 2;
      }
    else
      {
        MayDay::Error("bogus arguments to getmixedindex");
      }
    return retval;
  }

  void
  getMixedDerivDirections(int & a_derivDir1, int& a_derivDir2, const int& a_index)
  {
    //0 == xy, 1 = xz, 2 == yz
    if(a_index == 0)
      {
        a_derivDir1 = 0;
        a_derivDir2 = 1;
      }
    else if(a_index == 1)
      {
        a_derivDir1 = 0;
        a_derivDir2 = 2;
      }
    else if(a_index == 2)
      {
        a_derivDir1 = 1;
        a_derivDir2 = 2;
      }
    else
      {
        MayDay::Error("bogus arguments to getmixedderivdirections");
      }
  }
private:

  // Disallowed for all the usual reasons
  void operator=(const NWOQuadCFInterp&);
  NWOQuadCFInterp(const NWOQuadCFInterp&);
};

#include "NamespaceFooter.H"
#endif
