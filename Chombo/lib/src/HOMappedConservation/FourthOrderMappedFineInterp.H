#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERMAPPEDFINEINTERP_H_
#define _FOURTHORDERMAPPEDFINEINTERP_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Class to interpolate from coarse to valid fine cells on mapped grids
 *
 *//*+*************************************************************************/

#include <map>

#include "IntVectSet.H"
#include "FArrayBox.H"
#include "CHArray.H"
#include "LevelData.H"
#include "Copier.H"
#include "MultiBlockLevelCopierAverage.H"

#include "NamespaceHeader.H"

//--Forward declarations

class RealVect;
class ProblemDomain;
class BoxIterator;
class IVSIterator;
class CoordTransform;
class LevelGridMetrics;


/*******************************************************************************
 */
///  Interpolate from coarse cells to valid fine cells on mapped grids.
/**
 *   This class is used to interpolate the solution from coarse cells to valid
 *   fine cells on mapped grids using a constrained least-squares procedure.
 *   We solve a least-squares system for the coefficients of a polynomial
 *   interpolant
 *   \f[
 *     \sum \limits_{\mathbf{p}} a_{\mathbf{p}} \langle \mathbf{\xi}^\mathbf{p}
 *     \rangle_{\mathbf{j}^l} = \langle U \rangle_{\mathbf{j}^l} \hbox{ , }
 *     \mathbf{j}^l \in \mathcal{I}(\mathbf{i}^l)
 *   \f]
 *   subject to the constraint
 *   \f[
 *     \sum \limits_{\mathbf{i}^{l+1} \in \mathcal{C}^{-1} {\mathbf{i}^l }}
 *     \sum \limits_{\mathbf{p}} a_\mathbf{p} \langle J \mathbf{\xi}^\mathbf{p}
 *     \rangle_{\mathbf{i}^{l+1}} = \langle J U \rangle_{\mathbf{i}^l}.
 *   \f]
 *   Given this interpolant, we can construct \f$\langle JU \rangle\f$ on the
 *   control volumes at level \f$l+1\f$.
 *   \f[
 *     \langle JU \rangle_{\mathbf{i}^{l+1}} = \sum \limits_\mathbf{p}
 *     a_\mathbf{p} \langle J\mathbf{\xi}^\mathbf{p} \rangle_{\mathbf{i}^{l+1}}
 *     \hbox{ , } \mathbf{i}^{l+1} \in \mathcal{C}^{-1}(\mathbf{i}^{l})
 *   \f]
 *   The details of the method are described in the Chombo document titled
 *   "Fourth-order conservation law example"
 *
 *   <B> Implementation </B><BR>
 *   To compute the coefficients, we solve the system
 *   \f[
 *     \min_{{B \mathbf{x}} = \mathbf{d}} \| A \mathbf{x} - \mathbf{b} \|_2
 *   \f]
 *   where \f$A\f$ is the lenM x lenN coefficient matrix of
 *   \f$\langle \mathbf{\xi}^\mathbf{p} \rangle_{\mathbf{j}^l}\f$.  lenM is
 *   the number of equations (number of cells in the stencil, the maximum value
 *   is for the interior stencil where \f$\mbox{lenM} =
 *   3^{\mbox{SpaceDim}} + 2\mbox{SpaceDim}\f$ and lenN is the
 *   number of powers of a displacement given by 0 = binomial(degree +
 *   SpaceDim, SpaceDim) where degree is the order of interpolation (i.e., 3
 *   for a cubic polynomial and fourth-order interpolation).  We solve for
 *   coeffients \f$\mathbf{x} = a_\mathbf{p}\f$ where \f$B \mathbf{x} =
 *   \mathbf{d}\f$ is the constraint equation shown above.  You will see this
 *   notation in the code, e.g, A, LSx, LSBu, LSbl, etc.  Implementation is
 *   broken down into two parts.  During construction, data is initialized that
 *   is valid for any boxes with the same mesh spacing and refinement ratios.
 *   Specifically, this is A and the displacements in the fine cells.  Then,
 *   the member function, 'interpolate', will then interpolate everywhere on a
 *   coarse box.
 *
 *   The solution can be achieved using the LAPACK [S|D]GGLSE routine or an
 *   internal method which solves the linear system using a Cholesky factoriza-
 *   tion.  Define statements in the Fortran routine 'maplsinterp' are used to
 *   select the method.  The internal method is faster by a factor of about 4
 *   for 2 components in the solution and should be even faster for more
 *   components.  If both are used, the results are compared to some tolerance.
 *
 *   <B> Stencil </B><BR>
 *   The standard interior stencil has an inner box with size 3 and outer
 *   cells added in each dimension.  See the Chombo doc mentioned above for
 *   more details and boundary stencils.  The maximum stencil extent is shown
 *   below when cell 'i' is adjacent to the boundary.  The encroachment of the
 *   wall in this example is by two cells relative to cell 'i' in the interior
 *   stencil.  We label this as a "stencil offset from interior" of 2 in this
 *   direction, (2, 0).  In this fashion, we can mark all the possible stencils.
 *                                                                     \verbatim
 *           +---+                                   +---+\\
 *           |   |                                   |   |\\
 *       +===========+                       +===========+\\
 *       H   |   |   H                       H   |   |   H\\
 *   +---H---+---+---H---+               +---H---+---+---H\\
 *   |   H   | i |   H   |               |   H   |   | i H\\
 *   +---H---+---+---H---+               +---H---+---+---H\\
 *       H   |   |   H                       H   |   |   H\\
 *       +===========+                       +===========+\\
 *           |   |                                   |   |\\
 *           +---+                                   +---+\\
 *              <---2--> stencil extent    <-----3----> max stencil extent
 *                                                                  \endverbatim
 *  The stencils along with the coefficient matrices are stored during
 *  construction.  To minimize the required memory, symmetric stencils are only
 *  stored once and coordinate transformations are used during actual
 *  interpolation.  In the above example, a wall on the left would produce an
 *  identical but reflected stencil.  If dx = dy, then walls on the top and
 *  bottom would also form the same stencil with some rotational similarity.
 *  Without any symmetry, the number of required stencils is given by
 *  \f$5^{\mbox{SpaceDim}}\f$.  By exploiting symmetry in a given direction
 *  (quadrant), the number of stencils is \f$3^{\mbox{SpaceDim}}\f$.  By
 *  exploiting symmetry between directions with the same mesh spacing, the
 *  number of stencils is binomial(SpaceDim + 2, 2).  This has a dramatic effect
 *  at higher dimensions.  At 6 dimensions, the memory required for 64 bit
 *  storage of matrix A are 7.246 GiB, 346 MiB, and 13 MiB, respectively.  For
 *  3 dimensions, this is 644 KiB, 139 KiB, and 52 KiB, respectively.  All
 *  possible symmetries are exploited herein, both within a direction and
 *  between directions with the same mesh spacing.
 *
 *  \note
 *   <ul>
 *     <li> Most of the hard-core work is written in Fortran in file
 *          FourthOrderMappedFineInerpF.ChF.
 *     <li> loadAvgXipA is written in both C++ and Fortran exploring ways of
 *          using the SpaceDim macros.  Both seem to have similar performance.
 *     <li> The min box size in any given direction is 4.  11/24/10 I don't
 *          think this is true -- proper nesting is all that is required.
 *   </ul>
 *
 ******************************************************************************/

class FourthOrderMappedFineInterp
{
protected:

  /// Status of internal coarsened-fine data
  enum CrFnStatus
  {
    CrFnStatusInvalid = 0,
    CrFnStatusHaveJU = 1,
    CrFnStatusHaveU  = 2
  };

public:


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  /// Non-defining construction
  FourthOrderMappedFineInterp(const int a_degree,
                              const int a_nComp);

  /// Constructor that initializes the stencils
  FourthOrderMappedFineInterp(const int       a_degree,
                              const int       a_nComp,
                              const IntVect&  a_refRatio,
                              const RealVect& a_h);

  /// Destructor
  ~FourthOrderMappedFineInterp();

private:

  //  Copy and assignment not allowed
  FourthOrderMappedFineInterp(const FourthOrderMappedFineInterp&);
  FourthOrderMappedFineInterp& operator=(const FourthOrderMappedFineInterp&);


/*==============================================================================
 * Member functions
 *============================================================================*/

public:

  /// Define that initializes the stencils
  virtual void defineStencils(const IntVect&  a_refRatio,
                              const RealVect& a_h);

  /// Define that allocates data for coarse representations of the fine grid
  virtual void defineGrids(
    const DisjointBoxLayout&         a_FnGrid,
    const DisjointBoxLayout&         a_CrGrid,
    const ProblemDomain&             a_CrProblemDomain,
    const IntVect&                   a_FnInterpRadVec = IntVect::Zero,
    const bool                       a_willFillGhostsWithCrFnJU = true,
    const MultiBlockLevelGeom *const a_mbgeoPtr = 0);

  /// Define that saves the coarse <J> on the coarsened-fine layout
  virtual void defineCrFnJ(const LevelData<FArrayBox>& a_CrJ);

  /// Define that initializes the coarse-fine interface for filling ghosts
  virtual void defineCFInterface(LevelGridMetrics *const a_CrMetricsPtr);

//--For interpolating to an entire fine patch during regrid

  /// Interpolate from a coarse to fine level on a mapped grid
  virtual void interpToFine(LevelData<FArrayBox>&       a_FnLevJU,
                    const LevelData<FArrayBox>& a_FnLevJ,
                    const LevelData<FArrayBox>& a_CrLevU,
                    const LevelData<FArrayBox>& a_CrLevJU,
                    LevelGridMetrics *const     a_CrMetricsPtr,
                    const Interval&             a_vectorIntv);

//--For filling ghost cells at the coarse-fine interface

  /// Fill the invalid ghost cells on a mapped grid
  void fillGhosts(LevelData<FArrayBox>&       a_FnLevU,
                  LevelData<FArrayBox>&       a_FnLevJU,
                  const LevelData<FArrayBox>& a_FnLevJ,
                  const LevelData<FArrayBox>& a_CrLevData,
                  const char* const           a_CrDataType,
                  const ProblemDomain&        a_CrProblemDomain);

  /// Fill invalid ghost cells from preset coarsened-fine data on a mapped grid
  virtual void fillGhosts(LevelData<FArrayBox>&       a_FnLevU,
                  LevelData<FArrayBox>&       a_FnLevJU,
                  const LevelData<FArrayBox>& a_FnLevJ,
                  LevelGridMetrics *const     a_CrMetricsPtr);

//--Support for both types of operations

  /// Interpolate from a coarsened-fine to fine level on a mapped grid
  virtual void interpCrFnToFine(LevelData<FArrayBox>&       a_FnLevJU,
                        const LevelData<FArrayBox>& a_FnLevJ,
                        const DisjointBoxLayout&    a_CrFnGrid,
                        const LevelData<FArrayBox>& a_CrFnLevU,
                        const LevelData<FArrayBox>& a_CrFnLevJU,
                        LevelGridMetrics *const     a_CrMetricsPtr) const;

  /// Interpolate from a coarse to fine box on a mapped grid
  int interpToFine(FArrayBox&           a_FnFabJU,
                   const FArrayBox&     a_FnFabJ,
                   const Box&           a_CrBox,
                   const FArrayBox&     a_CrFabU,
                   const FArrayBox&     a_CrFabJU,
                   const ProblemDomain& a_CrProblemDomain) const;

  /// Fill any ghost cells on a finer mesh
  void fillGhosts(FArrayBox&           a_FnFabU,
                  FArrayBox&           a_FnFabJU,
                  const FArrayBox&     a_FnFabJ,
                  const IntVectSet&    a_CrCells,
                  const BaseFab<int>&  a_CrCellsArIdx,
                  const FArrayBox&     a_CrFabU,
                  const ProblemDomain& a_CrProblemDomain) const;

  /// True if stencil defined
  bool isStencilDefined() const;

  /// True if coarsend-fine grid data defined
  bool isGridDefined() const;

  /// True if coarse-fine interface defined
  bool isCFInterfaceDefined() const;

  /// True if coarse-fine \<J\> defined
  bool isCrFnJDefined() const;

  /// Check for proper nesting of the fine level in the coarse level
  bool properNesting(const DisjointBoxLayout& a_CrBoxes,
                     const DisjointBoxLayout& a_CrFnBoxes,
                     LevelGridMetrics *const  a_CrMetricsPtr) const;

//--Access

  /// Preset the coarsened-fine \<U\> externally
  virtual LevelData<FArrayBox>& presetCrFnLevU();

  /// Preset the coarsened-fine \<JU\> externally
  virtual LevelData<FArrayBox>& presetCrFnLevJU();

  /// Preset the coarsened-fine \<U\> from coarse level data
  virtual void presetCrFnLevU(const LevelData<FArrayBox>& a_CrLevU,
                      const Interval&             a_vectorIntv);

// Should not be used
  // /// Preset the coarsened-fine \<JU\> from coarse level data
  // void presetCrFnLevJU(const LevelData<FArrayBox>& a_CrLevJU);

// Obsolete
  // /// Access coarsened-fine \<U\> to fill externally
  // const LevelData<FArrayBox>& getCrFnLevU() const;

  // /// Access coarsened-fine \<JU\> to fill externally
  // const LevelData<FArrayBox>& getCrFnLevJU() const;

  /// Invalidate the CrFn data, both \<U\> and \<JU\>
  virtual void invalidateCrFnLevData();

  /// Get the number of ghosts used to build the coarsened-fine data
  const IntVect& CrFnNumGhost() const;

  /// Return the number of components
  int nComp() const;

private:

  void selectStencil(const IntVect        a_center,
                     const ProblemDomain& a_CrProblemDomain,
                     int&                 a_iStencil,
                     CoordTransform&      a_ctStoP,
                     CoordTransform&      a_ctPtoS,
                     IntVect&             a_offset,
                     IntVect&             a_fineSmlEnd) const;

  /// Sort sets of dimensions with the same mesh spacing
  void sortEachSetDx(int *const a_a, int *const a_b) const;

  /// Get the stencil key from a list of offsets
  int offset2Key(const int *const a_offsets) const;

  /// Load \f$\langle \Delta\xi \rangle^p\f$ for all p
  void loadAvgXipA(Real *a_A, const RealVect& a_dxi, const RealVect a_h) const;

  /// Load \f$\langle \Delta\xi \rangle^p\f$ for all p, evaluated with K
  void loadAvgXiKpA(Real *a_A, const RealVect& a_dxi, const RealVect a_h) const;

  /// Load gradients of \f$\langle \Delta\xi \rangle^p\f$ for all p
  //  Applies to both of the above methods for finding AvgXip
  void loadAvgGradXip(Real*           a_gradXip,
                      const RealVect& a_dxi,
                      const RealVect  a_h) const;

protected:

  /// Is this a multiblock grid
  bool isMultiblock() const;


/*==============================================================================
 * Data members
 *============================================================================*/

public:

  const int m_degree;                 ///< Degree of interpolating polynomial
  const int m_stencilExtent;          ///< Maximum extent of interior stencil.
                                      ///< Hard coded as 2.
  const int m_nComp;                  ///< Number of components for \<U\> and
                                      ///< \<JU\>
  IntVect m_refRatio;                 ///< Refinement ratios between the levels
  int m_numStencil;                   ///< Number of possible stencils
  int m_numUniqueDx;                  ///< Number of different mesh spacing over
                                      ///< all dimensions
  int m_uniqueDx[SpaceDim];           ///< Marks sets of dimensions with the
                                      ///< same mesh spacing
  bool m_defined[4];                  ///< 0 - stencil, 1 - grids,
                                      ///< 2 - CF interface,
                                      ///< 3 - coarsened-fine\<J\>
  bool m_willFillGhostsWithCrFnJU;    ///< T - expecting to fill ghosts from
                                      ///<     from a coarsened-fine \<JU\> as
                                      ///<     specified in defineGrids().
  std::map<int, int> m_mapStencilIndex;
                                      ///< Maps a key based on the stencil
                                      ///< offset to a sequential index for a
                                      ///< stencil

  // Note: All CHArrays have column storage to avoid confusion with Fortran
  CHArray<IntVectSet, 1> m_ivsStencilMask;
                                      ///< Masks for the various stencils.
  CHArray<CHMatrix, 1, ArZeroCol, ArSp::ArrayOfMatrixAlloc> m_At;
                                      ///< The transpose of matrix A for solving
                                      ///< \f$A\mathbf{x} = \mathbf{b}\f$.
                                      ///< I.e., the powers of the displacements
                                      ///< \f$\delta\xi\f$ to the neighbour
                                      ///< cells.  The 'm_At' for all stencils
                                      ///< are stored in a continguous chunk of
                                      ///< memory by the ArrayOfMatrixAlloc
                                      ///< allocator
  CHArray<Real, SpaceDim+1, ArRangeCol>* m_fnXiLSx;
                                      ///< Fine displacements already
                                      ///< multiplied by the solution to the
                                      ///< constrained least-squares problem
                                      ///< (LSx).  Multiply by LSbl to get \<U\>
                                      ///< in a fine cell.
  CHArray<Real, SpaceDim+1, ArRangeCol>* m_fnGradXiLSx;
                                      ///< Gradients of fine displacements
                                      ///< already multiplied by the solution to
                                      ///< the constrained least-squares problem
                                      ///< (LSx).  Multiply by LSbl to get the
                                      ///< gradient of \<U\> in a fine cell.
  CHArray<Real, SpaceDim+1, ArRangeCol> m_avgXipFine;
                                      ///< Powers of displacements for fine
                                      ///< cells
  CHArray<Real, SpaceDim+1, ArRangeCol> m_avgXiKpFine;
                                      ///< Powers of displacements for fine
                                      ///< cells with constant added to impose
                                      ///< constraint when filling ghosts
  CHArray<Real, SpaceDim+1, ArRangeCol> m_avgGradXipFine;
                                      ///< Gradients of 'm_avgXipFine'.  This is
                                      ///< stored linearly in each cell as
                                      ///< (SpaceDim, p) with 'p' contiguous
                                      ///< since fortran only supports rank
                                      ///< 7 arrays (at least until F2008
                                      ///< becomes commonplace).

  // Same as FineInterp -- used once during a regrid but also when filling
  // ghost cells
  DisjointBoxLayout m_CrFnGrid;       ///< The coarsened-fine boxes
  LevelData<FArrayBox> m_CrFnLevU;    ///< \<U\> on a coarse representation of
                                      ///< the fine mesh
  LevelData<FArrayBox> m_CrFnLevJU;   ///< \<JU\> on a coarse representation of
                                      ///< the fine mesh
  LevelData<FArrayBox> m_CrFnLevJ;    ///< \<J\> on a coarse representation of
                                      ///< the fine mesh
  Copier m_copier;                    ///< Friendly neighbourhood copier for
                                      ///< copying from Cr to CrFn
  MultiBlockLevelCopierAverage m_mbcp;
                                      ///< To copy information between
                                      ///< diffferent layout across multiblock
                                      ///< boundaries

  // Data for filling interface ghost cells
  IntVect m_FnInterpRadVec;           ///< Number of fine ghost cells that must
                                      ///< be filled in each direction.
  IntVect m_CrInterpRadVec;           ///< Number of coarse cells that underlie
                                      ///< fine ghost cells in each direction
  IntVect m_CrFnNumGhost;             ///< Number of layers of ghost cells
                                      ///< required on the coarsened-fine mesh
  LayoutData<IntVectSet> m_CrFnFillInterp;
                                      ///< Coarse cells in which to interpolate
                                      ///< to fill fine ghost cells
  LayoutData<BaseFab<int> > m_ArIdxFillInterp;
                                      ///< Index into all the CHArray for an
                                      ///< IntVect
  CHArray<Box, 1, ArZeroCol> m_FnBoxFillInterp;
                                      ///< Box describing the fine cells that
                                      ///< need to be filled in a coarse cell.
                                      ///< Max size is n_{ref}^{l}
  int m_numCrFnFillGhost;             ///< Total number of coarse cells used to
                                      ///< fill fine ghost cells on this
                                      ///< processor (this is actually a count
                                      ///< of total cells, not layers like
                                      ///< numGhost usually is)
  int m_CrFnStatus;                   ///< Status of coarsened-fine data.  We
                                      ///< may know JU, U, or nothing.  This is
                                      ///< kept consistent with what
                                      ///< MappedLevelData knows about the
                                      ///< coarser data
  bool m_isMultiblock;                ///< T - A multiblock coordinates system
};


/*******************************************************************************
 *
 * Class FourthOrderMappedFineInterp: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//   True if stencil defined
/**  \return            T - stencil has been defined
 *//*-----------------------------------------------------------------*/

inline bool
FourthOrderMappedFineInterp::isStencilDefined() const
{
  return m_defined[0];
}

/*--------------------------------------------------------------------*/
//   True if coarsend-fine grid data defined
/**  \return            T - coarsend-fine grid data has been defined
 *//*-----------------------------------------------------------------*/

inline bool
FourthOrderMappedFineInterp::isGridDefined() const
{
  return m_defined[1];
}

/*--------------------------------------------------------------------*/
//   True if coarse-fine interface defined
/**  \return            T - coarse-fine interface has been defined
 *//*-----------------------------------------------------------------*/

inline bool
FourthOrderMappedFineInterp::isCFInterfaceDefined() const
{
  return m_defined[2];
}

/*--------------------------------------------------------------------*/
//   True if coarse-fine <J> defined
/**  Must be defined to fill ghosts using 'fillGhostsFromCrJU'
 *   \return            T - m_CrFnLevJ has been defined
 *//*-----------------------------------------------------------------*/

inline bool
FourthOrderMappedFineInterp::isCrFnJDefined() const
{
  return m_defined[3];
}

/*--------------------------------------------------------------------*/
//  Preset the coarsened-fine \<U\> externally
/** This returns the coarsened-fine representation and you have to do
 *  the copyTo externally
 *//*-----------------------------------------------------------------*/

inline LevelData<FArrayBox>&
FourthOrderMappedFineInterp::presetCrFnLevU()
{
  CH_assert(isGridDefined());
  m_CrFnStatus |= CrFnStatusHaveU;
  return m_CrFnLevU;
}

/*--------------------------------------------------------------------*/
//  Preset the coarsened-fine \<JU\> externally
/** This returns the coarsened-fine representation and you have to do
 *  the copyTo externally
 *//*-----------------------------------------------------------------*/

inline LevelData<FArrayBox>&
FourthOrderMappedFineInterp::presetCrFnLevJU()
{
  CH_assert(isGridDefined());
  m_CrFnStatus |= CrFnStatusHaveJU;
  return m_CrFnLevJU;
}

// Obsolete
// /*--------------------------------------------------------------------*/
// //  Access coarsened-fine \<U\>
// /*--------------------------------------------------------------------*/

// inline const LevelData<FArrayBox>&
// FourthOrderMappedFineInterp::getCrFnLevU() const
// {
//   CH_assert(isGridDefined());
//   return m_CrFnLevU;
// }

// /*--------------------------------------------------------------------*/
// //  Access coarsened-fine \<JU\>
// /*--------------------------------------------------------------------*/

// inline const LevelData<FArrayBox>&
// FourthOrderMappedFineInterp::getCrFnLevJU() const
// {
//   CH_assert(isGridDefined());
//   return m_CrFnLevJU;
// }

/*--------------------------------------------------------------------*/
//  Invalidate the CrFn data, both <U> and <JU>
/*--------------------------------------------------------------------*/

inline void
FourthOrderMappedFineInterp::invalidateCrFnLevData()
{
  m_CrFnStatus = CrFnStatusInvalid;
}

/*--------------------------------------------------------------------*/
//  Get the number of ghosts used to build the coarsened-fine data
/*--------------------------------------------------------------------*/

inline const IntVect&
FourthOrderMappedFineInterp::CrFnNumGhost() const
{
  CH_assert(isGridDefined());
  return m_CrFnNumGhost;
}

/*--------------------------------------------------------------------*/
//  Return the number of components
/*--------------------------------------------------------------------*/

inline int
FourthOrderMappedFineInterp::nComp() const
{
  return m_nComp;
}

/*--------------------------------------------------------------------*/
//  Is this a multiblock grid?
/** \return             T - Yes
 *//*-----------------------------------------------------------------*/

inline bool
FourthOrderMappedFineInterp::isMultiblock() const
{
  return m_isMultiblock;
}

#include "NamespaceFooter.H"
#endif
