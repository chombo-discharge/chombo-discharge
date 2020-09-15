#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERMAPPEDCFINTERP_H_
#define _FOURTHORDERMAPPEDCFINTERP_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "FourthOrderMappedCFInterpStencil.H"

#include "NamespaceHeader.H"

/// Fourth-order interpolation in space

/**
 */
class FourthOrderMappedCFInterp
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  FourthOrderMappedCFInterp();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~FourthOrderMappedCFInterp();

  /// Actual constructor.
  /**
     Set up object.
   */
  void define(/// layout at finer level
              const DisjointBoxLayout&  a_layout,
              /// number of variables
              int                       a_numStates,
              /// refinement ratio between finer level and coarser level
              int                       a_refineCoarse,
              /// problem domain on finer level
              const ProblemDomain&      a_domain,
              /// number of coarse-level ghost cells
              int                       a_coarseGhostsFill = 0,
              /// dimensions that are fixed, not interpolated
              Interval                  a_fixedDims = Interval() );

  void setJ(
            /// <J> on fine grids
            const LevelData<FArrayBox>&  a_FnLevJ,
            /// <J> on coarsened fine grids
            const LevelData<FArrayBox>&  a_CrFnLevJ);

  void setCrFnLevUandJU(
                        /// <U> on coarsened fine grids
                        const LevelData<FArrayBox>&   a_CrFnLevU,
                        /// <JU> on coarsened fine grids
                        const LevelData<FArrayBox>&   a_CrFnLevJU);

  /**
     Interpolate to <JU> on fine level.
   */
  void interpToFineJU(/// interpolated fine <JU>
                      LevelData<FArrayBox>&         a_fineJU);

protected:

  void setStencil(BaseFab<IntVect>&  a_stencil,
                  const Box&         a_coarseBox);

  void getGradient(FArrayBox&        a_gradFab,
                   const FArrayBox&  a_srcFab,
                   const Box&        a_bx);

  void clearMemory();

  void freeJ();

  /// Problem domain - index space for this level
  ProblemDomain m_domain;

  // m_domain coarsened by m_refineVect
  ProblemDomain m_coarseDomain;

  /// Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  /// dimensions that are fixed, not interpolated
  Interval m_fixedDims;

  /// number of dimensions that are interpolated: SpaceDim - m_fixedDims.size()
  int m_nInterpDims;

  /// 1 in m_fixedDims, m_refineCoarse in other dimensions
  IntVect m_refineVect;
  
  /// 0 in m_fixedDims, 1 in other dimensions
  IntVect m_interpUnit;

  /// Number of variables
  int m_numStates;

  /// maximum distance of stencil from domain boundary
  int m_maxStencilDist;

  /// (-m_maxStencilDist*interpUnit, m_maxStencilDist*interpUnit)
  Box m_stencilBox;

  /// max degree of approximations
  int m_degree;

  /// number of nonzero powers
  int m_numNonzeroPowers;

  /// number of layers of ghost cells at the coarser level whose fine cells are to be filled in
  int m_coarseGhostsFill;

  /// Box layout for finer level
  DisjointBoxLayout m_layout;

  /// Box layout for coarser level:  m_layout coarsened by m_refineVect
  DisjointBoxLayout m_layoutCoarsened;

  /// all possible stencils, on (-m_maxStencilDist:+m_maxStencilDist)^SpaceDim
  BaseFab<FourthOrderMappedCFInterpStencil*> m_stencils;

  /// on m_layoutCoarsened; which stencil to use at each coarse point
  LevelData< BaseFab<IntVect> > m_stencilHere;

  /// <(z - z0)^p>k for each fine cell k.
  FArrayBox m_fineAvgPowers;

  /// gradient of m_fineAvgPowers
  FArrayBox m_gradFineAvgPowers;

  /// ghost vector on m_coarseData
  IntVect m_ghostVect;

  /// ghost vector on coarsened fine data
  IntVect m_coarseGhostVect;

  /// IntVect::Zero to m_refineVect-interpUnit
  Box m_fineBox;

  /// copied coarse data in interpToFine, lives on m_layoutCoarsened
  LevelData<FArrayBox> m_coarseData;

  /// <J*(z - z0)^p>k for each fine cell k.
  LevelData<FArrayBox> m_JfineAvgPowers;

  /// B(i) = sum_{k in F(i)} (<J*zvec**pvec>k); m_numNonzeroPowers components
  LevelData<FArrayBox> m_B;

  /// 2-norm of B(i); 1 component
  LevelData<FArrayBox> m_normB;

  LevelData<FArrayBox> m_alpha;

  LevelData<FArrayBox> m_beta;

  LevelData<FArrayBox> m_QRtau;

  LevelData<FArrayBox> m_FnLevJ;

  LevelData<FArrayBox> m_coeff0;

  LevelData<FArrayBox> m_coeffs;

  LayoutData< BaseFab< Vector<Real>* >* > m_QRmatrix;

  /// Vector<Real> *m_A1[dit](cell, 0) has a component for each stencil cell
  LayoutData < BaseFab< Vector<Real>* >* > m_A1;

  /// define() has been called
  bool m_defined;

  /// J has been set, together with all the arrays calculated from J
  bool m_haveJ;

  /// coefficients has been set, from call to setCrFnLevUandJU
  bool m_haveCoeffs;

private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderMappedCFInterp&);
  FourthOrderMappedCFInterp(const FourthOrderMappedCFInterp&);
};

#include "NamespaceFooter.H"

#endif
