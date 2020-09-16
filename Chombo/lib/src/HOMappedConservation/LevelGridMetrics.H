#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELGRIDMETRICS_H_
#define _LEVELGRIDMETRICS_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Grid metrics on and between levels.
 *
 *//*+*************************************************************************/

#include "AMRLevel.H"
#include "CodimBox.H"
#include "CoarseAverageFace.H"
#include "CoarseAverageCodim.H"
#include "MultiBlockCoordSys.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderMappedFineInterp.H"
#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "MultiBlockFaceRegister.H"

#include "NamespaceHeader.H"

//--Forward declarations

class TimeInterpolatorRK4;  // (part of timeIntermediate hack)

// //--Debug support

// extern int LGMDebug;


/*******************************************************************************
 */
///  Average from fine intersections
/**
 *   The objective of this class is to average down information from A^(l+1)
 *   to C^(l) but only on the intersection of A^(l+1) and B^(l+1).  This class
 *   is designed to support class 'LevelGridMetrics' where A is the snapback
 *   flux on the old fine mesh, B is the new mesh, and C (= A^l) is the snapback
 *   flux on the old coarse mesh.  These are the steps (Cr is coarsen operator):
 *   <ol>
 *     <li> Average A^(l+1) to D^(l) on boxes of Cr(A^(l+1)).  This step is
 *          performed at level 'l+1'.  D and the coarsened box layout of A are
 *          saved.
 *     <li> Now at level 'l', copy to E^(l) on boxes of C^(l).  These first two
 *          steps are done in the averager
 *     <li> Compute mask for each box of C^(l).  The mask is the
 *          intersection of Cr(A^(l+1)), Cr(B^(l+1)), and C^(l)
 *     <li> Copy from E^(l) to C^(l) using the mask.
 *   </ol>
 *
 ******************************************************************************/

class IntersectionAverageFace
{
public:

  /// Constructor
  IntersectionAverageFace(const LevelData<FluxBox>& a_A,
                          const int                 a_nRef);

  // Use synthesized destructor

  /// Copy the average using a mask
  void copyTo(LevelData<FluxBox>&      a_C,
              const DisjointBoxLayout& a_dblB);

private:

  // Copy and assignment not allowed
  IntersectionAverageFace(const IntersectionAverageFace&);
  IntersectionAverageFace& operator=(const IntersectionAverageFace&);

protected:

  const int m_nComp;                  ///< Number of components in A (same as C)
  const int m_nRef;                   ///< Refinement between A & C
  DisjointBoxLayout m_crDblA;         ///< Coarsened layout for A
  CoarseAverageFace m_averager;       ///< The averager itself stores the
                                      ///< averaged data between averaging and
                                      ///< copying
};


/*******************************************************************************
 */
///  Grid metrics for a level
/**
 *   This class adds infrastructure for mapped grids to the fourth-order
 *   Cartesian AMRLevel class.  The interface closely follows that of
 *   AMRLevel.  Typically this class is a member of the user's derived
 *   AMRLevel and the member functions of LevelGridData should be called
 *   from member functions of the derived AMRLevel that have the same name.
 *
 *   <p> This class provides m_N, the metric terms \f$\mathbf(N)\f$ on each
 *   face, and m_J, the determinant (physical cell volume divided by
 *   computational cell volume).  These quantities are appropriately maintained
 *   during all aspects of an AMR solution.
 *
 *   <p> GHOSTS --- The determinants, (J), are defined on all ghost cells at
 *   the domain interior and 1 layer of ghost cells outside a non-periodic
 *   domain.  Hence, gradients of J are available everywhere an application may
 *   store data.  To achieve this, it is assumed that the metrics can be
 *   evaluated everywhere in \f$\mathcal(G)(\Gamma,2)\f$.
 *   HOWEVER, metrics that are consistent between levels are only maintained in
 *   the valid region of the finer level.  A valid cell on a coarse level will
 *   have metrics computed separately from any overlapping ghost cells from a
 *   finer level (i.e., the volumes of the fine ghost cells may not sum to the
 *   volume of the valid coarse cell).  Consistency between levels everywhere
 *   could be achieved by treating all ghost cells as part of the valid region
 *   --- this would simply (probably simply) require some modifications to the
 *   "averagers".  However we do not do this because:
 *   <ul>
 *     <li> it breaks a consistent definition of valid/invalid
 *     <li> there is no requirement for conservation in ghost cells since they
 *          are only used for gradient reconstruction, limiting, etc.
 *   </ul>
 *   Any ghost cell that is covered by a valid cell on that level is consistent
 *   with the metrics from any finer level by means of the copiers (used in
 *   the averagers) copying to both the invalid and valid regions of a
 *   destination.
 *   <p>
 *   Transverse terms in N are determined only on the faces of grow(N.box(), -1)
 *   if obtained by averaging, i.e., on one less layer than where the normal
 *   components are known.  If obtained by analytic methods, they are available
 *   everywhere.
 *
 *   <p> Normally, N is calculated using script N to ensure freestream
 *   preservation.  But in some cases, e.g., a 2-D solution on the surface of
 *   a sphere, there is no definition of script N.  In those cases,
 *   m_useScriptN should be set to false and N will be calculated and averaged
 *   directly.  The averaging ensures conservation but freestream-preservation
 *   is sacrificed.  The value of m_useScriptN is obtained from the
 *   MultiBlockCoordSys and should be set to 'false' in the constructor of
 *   coordinate systems derived from class MultiBlockCoordSys.  Otherwise, by
 *   default, it is set to true.
 *
 *   <p>
 *   The following routines should be called from the corresponding routines
 *   with the same name derived AMRLevel:
 *   <ul>
 *     <li> initialGrid
 *     <li> postInitialGrid
 *     <li> postTimeStep
 *     <li> preRegrid
 *     <li> regrid
 *   </ul>
 *
 *   <p> Setup for the class is as follows:
 *   <ol>
 *     <li> Call define to setup the class members -- metrics are still
 *          undefined.
 *     <li> Optionally call initialGrid to set the new box layout on a level.
 *          This can be deferred until postInitialGrid if desired.
 *     <li> Call postInitialGrid to set the new box layout (if not done in
 *          previous step) and define the metrics.
 *   </ol>
 *
 *   <p> During a run, be sure to:
 *   <ol>
 *     <li> Call postTimeStep to average down the metric terms from finer
 *          levels.  This only has an effect if they changed on the finer
 *          level
 *   </ol>
 *
 *   <p> For a regrid do the following:
 *   <ol>
 *     <li> Call preRegrid to compute the new metrics and correct the solution
 *          on the old mesh to the metrics of the new mesh (snapback).
 *     <li> Optionally call regrid to interpolate from the coarse mesh to
 *          the new fine mesh.  This can be done is user code but must
 *          use FourthOrderMappedFineInterp to ensure a conservative and
 *          free-stream preserving interpolation on the mapped grid.
 *   </ol>
 *
 *   \note
 *   <ul>
 *     <li> To facilitate products that involve gradients, the metrics (N) are
 *          determined on the faces of boxes with size ghostBox+2.  This is
 *          labeled an 'NReqGhostBox' herein.  At domain boundaries, the boxes
 *          for 'N' extend by 2 cells outside the physical domain and it is
 *          assumed that the mappings exist there. The determinants (J) are
 *          determined on the cells of boxes with size ghostBox+1.  This is
 *          labeled a 'JReqGhostBox' herein.  At domain boundaries, the boxes
 *          for 'J' extends by 1 cell outside the physical domain
 *     <li> Averages of fluxes.  Most of the fluxes, \<F>\ are stored as average
 *          quantities to be applied to a face and do not include area
 *          information.  I.e., the actual flux across a 2-D face is
 *          \f$\langle F\rangle h^2\f$  Hence when correcting from fine to
 *          coarse, you want to average the fine \<F\> and then muliply by the
 *          area of the coarse face (although in most equations, the area
 *          cancels out with another term).  However, the \f$\mathcal(N)^s\f$,
 *          which are a type of flux, are stored differently in that they have
 *          already been multiplied by the length term (equivalent to
 *          \f$\langle F\rangle h^2\f$ in the example above).  Hence, these
 *          fluxes needed to be summed when correcting from coarse to fine.
 *     <li> WARNING!  If transverse components of <N> are computed analytically,
 *          they are only averaged down when finer levels have changed before
 *          entering a postTimeStep (along with normal components).
 *          Consequently, it is currently difficult to know if they have been
 *          computed or averaged from a finer level at any given point.  But
 *          note that if transverse components of <N> are computed by averaging
 *          from normal components, then they should always be consistent.
 *          **FIXME -- probably the best fix is to only average the normal
 *          components if transverse are computed analytically.
 *   </ul>
 *
 *//*+*************************************************************************/

class LevelGridMetrics
{

public:

  /// Age of metric terms during a regrid
  enum AgeType
  {
    AgeOld = 0,
    AgeNew = 1,
    AgeNum = 2
  };

  /// Type of operation for computing transverse terms in N
  enum TransverseNOpType
  {
    TransverseNOpNone,
    TransverseNOpAnalytic,              // Analytic calculation
    TransverseNOpAverage                // Average from normals on nearby
                                        // transverse faces
  };


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  LevelGridMetrics(const int a_nComp,
                   const int a_spaceOrder);

  /// Destructor
  virtual ~LevelGridMetrics();

private:

  // Copy and assignment not allowed
  LevelGridMetrics(const LevelGridMetrics&);
  LevelGridMetrics& operator=(const LevelGridMetrics&);


/*==============================================================================
 * Member functions
 *============================================================================*/

public:

  /// Free memory
  void clear();

  /// Does a finer level exist?
  bool hasFiner() const;

  /// Does a coarser level exist?
  bool hasCoarser() const;

  /// Define the class.  No metrics terms are yet available.
  void define(
    AMRLevel *const                  a_AMRLevel,
    MultiBlockCoordSysFactory *const a_coordSysFact,
    LevelGridMetrics *const          a_coarserLevelGridMetrics,
    const RealVect&                  a_dxVect,
    const IntVect&                   a_ghostVect,
    const TransverseNOpType          a_transverseNOpType = TransverseNOpNone,
    const bool                       a_haveMultiBlockVectorData = true,
//**FIXME do we need this?
    const int                        a_numBlocks = 1);

  /// Define the class.  No metrics terms are yet available.
  void define(
    AMRLevel *const                  a_AMRLevel,
    MultiBlockCoordSys *const        a_coordSys,
    LevelGridMetrics *const          a_coarserLevelGridMetrics,
    const RealVect&                  a_dxVect,
    const IntVect&                   a_ghostVect,
    const TransverseNOpType          a_transverseNOpType = TransverseNOpNone,
    const bool                       a_haveMultiBlockVectorData = true,
//**FIXME do we need this?
    const int                        a_numBlocks = 1);

  /// Are the metrics defined on this level?
  bool metricsDefined() const;

  /// Did the metrics change on the finer mesh?
  bool didFinerChange() const;

  /// Retrieve the multi-block coordinate system.
  const MultiBlockCoordSys& getCoordSys() const;

  /// Get the coordinate system
  /* Are we ever going to use a second-order system?  If so, I assume the
   * CoordSys routines will be virtualized such that getCoordSys2 can be
   * changed to getCoordSys and the cast will be unnecessary.
   */
  const NewCoordSys* getCoordSys2(const Box& box) const;

  ///
  const NewFourthOrderCoordSys* getCoordSys(const Box& box) const;

  /// Get the block index
  int getBlock(const Box& a_box) const;

  /// Get a problem domain representing a block
  const ProblemDomain& blockDomain(const Box& a_box, const int a_numGhost = 0);

  /// Get the computational mesh spacing
  const RealVect& dxVect() const;

//--AMRLevel extensions

  /// After a time step, average down \f$N\f$ and \f$J\f$
  void postTimeStep();

  /// Set up initial grid
  virtual void initialGrid(const DisjointBoxLayout *const a_newGrid);

  /// Set up the initial metric terms.
  void postInitialGrid(const DisjointBoxLayout *const a_newGrid);

  /// Compute new metrics and correct the solution on the coarser meshes
  void preRegrid(const int                         a_baseLevel,
                 const DisjointBoxLayout&          a_newGrid,
                 const LevelData<FArrayBox> *const a_coarseUOldPtr,
                 LevelData<FArrayBox>&             a_UOld,
                 LevelData<FArrayBox>&             a_JUOld);

  /// Regrid operations -- performs a coarse to fine interpolation
  virtual void regrid(LevelData<FArrayBox>&       a_JU,
                      const LevelData<FArrayBox>& a_CrU,
                      const LevelData<FArrayBox>& a_CrJU,
                      const Interval&             a_vectorIntv);

//   /// Post regrid operations -- currently unused.
//   void postRegrid(int a_base_level);

//--Convenience

  /// Fill invalid fine ghost cells (with preset CrFnU)
  virtual void fillFineGhostCells(LevelData<FArrayBox>& a_U,
                                  LevelData<FArrayBox>& a_JU,
                                  Real a_time = -1);

  /// Computes \<U\> in valid cells
  void computeValidU(LevelData<FArrayBox>&       a_U,
                     const LevelData<FArrayBox>& a_JU);

  /// Compute the minimum grid buffer size for a fourth-order interpolation
  static int bufferSize4thO(const std::vector<int>& a_refRatio,
                            const int               a_maxLevel,
                            const int               a_numGhost);

  /// Exchange \<U\> across multiblock boundaries
  void multiblockExchangeU(LevelData<FArrayBox>& a_U,
                           const Interval&       a_vectorIntv) const;

  /// Supports a call to TimeInterpolatorRK4::intermediate for multi-block grids
  void timeIntermediate(const TimeInterpolatorRK4& a_timeInterpolator,
                        const Real&                a_timeInterpCoeff,
                        const int&                 a_stage,
                        const Interval&            a_vectorIntv);

//--Consistent access to structural boxes

  /// Get the boxes
  const DisjointBoxLayout& getBoxes() const;

  /// Get a level iterator
  DataIterator getDataIterator() const;

  /// Resize a box to include the ghost cells
  Box resizeWithGhosts(const Box& box) const;

  /// Determine the ghost vector required for \f$J\f$.
  virtual IntVect getJReqGhostVect() const;

  /// Determine the ghost vector required for \f$N\f$.
  virtual IntVect getNReqGhostVect() const;

  /// Resize a box to the size required for \f$J\f$
  Box resizeWithJReqGhosts(const Box& box) const;

  /// Resize a box to the size required for \f$N\f$
  Box resizeWithNReqGhosts(const Box& box) const;

//--Access to the interpolator's temporary structures

  /*
   * The following return a coarsened representation of 'this' level, i.e.,
   * 'this' is the fine level.
   */

  /// Set the coarsened-fine \<U\> used by the interpolator (CrFn access)
  LevelData<FArrayBox>& presetCr2ThisInterpolatorCrFnLevU();

  /// Set the coarsened-fine \<JU\> used by the interpolator (CrFn access)
  LevelData<FArrayBox>& presetCr2ThisInterpolatorCrFnLevJU();

// Obsolete
  // /// Get the coarsened-fine \<U\> used by the interpolator
  // const LevelData<FArrayBox>& getCr2ThisInterpolatorCrFnLevU() const;

  // /// Get the coarsened-fine \<JU\> used by the interpolator
  // const LevelData<FArrayBox>& getCr2ThisInterpolatorCrFnLevJU() const;

  /*
   * The following make use of the coarse data on 'this' level, i.e., 'this' is
   * the coarse level.
   */

  /// Set the coarsened-fine \<U\> used by the interpolator
  void presetThis2FnInterpolatorCrFnLevU(
    const LevelData<FArrayBox>& a_CrLevU,
    const Interval&             a_vectorIntv);

// Should not be used
  // /// Set the coarsened-fine \<JU\> used by the interpolator
  // void presetThis2FnInterpolatorCrFnLevJU(
  //   const LevelData<FArrayBox>& a_CrLevJU);

  /// Invalidate the CrFn data, both \<U\> and \<JU>\, used by the interpolator
  void invalidateCr2ThisInterpolatorCrFnLevData();

  /// Invalidate the CrFn data, both \<U\> and \<JU>\, used by the interpolator
  void invalidateThis2FnInterpolatorCrFnLevData();

  /// Number of ghosts used to build the coarsened-fine data in the interpolator
  const IntVect& interpolatorCrFnNumGhost(const bool a_force = false) const;

  /// Is this a multiblock grid
  bool isMultiblock() const;

  ///
  const ProblemDomain& problemDomain() const;

protected:

//--Easier access to AMRLevel data (these are also defined in AMRLevel.H).

  ///
  int level() const;

  ///
  int refRatio() const;

  ///
  static int verbosity();

//--Internal

  /// Completely define the metrics for this level
  void defineMetrics();

  /// Label this and all finer levels as having undefined metrics
  void undefineMetrics();

  /// Compute \f$N\f$ on the hyperfaces from \f$\mathcal{N}^s\f$
  //  Only used if m_useScriptN == true
  void faceNormalN(LevelData<FluxBox>&        a_N,
                   const LevelData<CodimBox <FArrayBox> >& a_scrN);

  /// Compute \f$N\f$ on the hyperfaces from \f$\mathcal{N}^s\f$
  //  Only used if m_useScriptN == false
  void getN(LevelData<FluxBox>& a_N);

  /// Add transverse components to \f$N\f$
  void faceTransverseN(LevelData<FluxBox>& a_N);

  /// Integrate \f$\mathcal{N}^s\f$ on each hyperedge
  void getScriptN(LevelData<CodimBox<FArrayBox> >& a_scrN);

  /// Average \f$\mathcal{N}^s\f$
  //  Only used if m_useScriptN == true
  void averageScriptN(const bool           a_isBase,
                      LevelData<CodimBox<FArrayBox> >& a_scrN);

  /// Average \f$N\f$ directly if \f$\mathcal{N}^s\f$ doesn't exist.
  //  Only used if m_useScriptN == false
  void averageN(const bool          a_isBase,
                LevelData<FluxBox>& a_N);

  /// Compute the volume flux
  void getVolFlux(LevelData<FluxBox>&        a_NtX,
                  const LevelData<FluxBox>&  a_N);

  /// Take divergence of \f$N^T X\f$ to set m_J.
  virtual void setAvgJfromNtX(LevelData<FluxBox>& a_NtX);

  /// Average the volume flux \f$N^T X\f$.
  void averageVolFlux(const bool          a_isBase,
                      LevelData<FluxBox>& a_NtX,
                      const AgeType       a_age);

  /// Compute \f$N^T X\f$ on all finer levels and average down
  void volFluxCascade();

  /// Average the snapback solution flux, \f$N^T F\f$.
  void averageSnapbackSolFlux(const bool          a_isBase,
                              LevelData<FluxBox>& a_NtFDiff);

  /// Store the DisjointBoxLayout and define multiblock structures
  virtual void defineGrids(const DisjointBoxLayout& a_grids);

  /// Create a non-periodic DisjointBoxLayout
  void defineGridsNP(const DisjointBoxLayout& a_grids);


/*==============================================================================
 * Data members
 *============================================================================*/

public:

  LevelData<FluxBox> m_N;             ///< \f$\mathbf(N)\f$ on each face
  LevelData<FArrayBox> m_J;           ///< Physical cell volume divided by
                                      ///< computational cell volume

protected:

//--Single-block

  // The following averagers are not defined on the finest level
  CoarseAverageCodim* m_scrNAverager; ///< Averager for \f$\mathcal(N)^s\f$.
                                      ///< Only used if m_useScriptN == true
  CoarseAverageFace* m_NAverager;     ///< Averager for \f$N\f$
                                      ///< Only used if m_useScriptN == false
  CoarseAverageFace* m_NtXAverager[AgeNum];
                                      ///< This averages the volume flux from a
                                      ///< finer level.  Averagers exist for old
                                      ///< "AgeOld" and new "AgeNew" meshes
                                      ///< during a regrid
  IntersectionAverageFace* m_NtFDiffAverager;
                                      ///< Averager for the snapback solution
                                      ///< flux.  This averager only operates
                                      ///< on the intersections of three
                                      ///< layouts.

  MultiBlockCoordSys* m_coordSys;     ///< The multi-block coordinate system.

  LevelGridMetrics* m_finerLevelGridMetrics;
                                      ///< Grid metrics on the next finer level
  LevelGridMetrics* m_coarserLevelGridMetrics;
                                      ///< Grid metrics on the next coarser
                                      ///< level

  AMRLevel* m_parentAMRLevel;         ///< The AMR level we are attached to

  IntVect m_ghostVect;                ///< Ghost cell vector
  DisjointBoxLayout m_grids;          ///< Grid layout
  DisjointBoxLayout m_gridsNP;        ///< Non-periodic grid layout.  This is a
                                      ///< shallow copy of m_grids, with a non-
                                      ///< periodic problem domain and a new
                                      ///< neighbor iterator if problemDomain()
                                      ///< is periodic.
  RealVect m_dxVect;                  ///< Mesh spacing in each direction
  const int m_spaceOrder;             ///< Spatial-order of accuracy - must be 4

//--AMR tools

  FourthOrderMappedFineInterp* m_interpolator;
                                      ///< Constrained least-squares
                                      ///< interpolator for interpolating \<JU\>
                                      ///< from coarse to fine meshes.  'This'
                                      ///< mesh is the coarse mesh.

//--Multiblock tools

  ProblemDomain m_cachedBlockDomain;  ///< Representing only the block (it is
                                      ///< assumed periodic boundaries are not
                                      ///< used for multiblock)
  int m_cachedBlockIdx;               ///< Index of the cached block
  int m_cachedBlockNG;                ///< Number of ghosts that the cached
                                      ///< block was grown by.
  MultiBlockLevelGeom m_mbgeo;        ///< Topology of the multiblock layout
  MultiBlockLevelExchangeAverage m_mbex;
                                      ///< To exchange information across
                                      ///< multiblock boundaries
  MultiBlockFaceRegister* m_Nmbreg;   ///< Flux-like register for resolving
                                      ///< \f$\mathbf{N}\f$ on faces at coarse-
                                      ///< fine interfaces that overlap
                                      ///< multiblock boundaries
  MultiBlockFaceRegister* m_NtXmbreg[AgeNum];
                                      ///< Flux-like register for resolving
                                      ///< \f$\mathbf{N}^T\mathbf{X}\f$ on faces
                                      ///< at coarse-fine interfaces that
                                      ///< overlap multiblock boundaries

//--Miscellaneous/flags

  TransverseNOpType m_transverseNOpType;
                                      ///< Type of operatation for computing
                                      ///< transverse terms in N
  bool m_metricsChanged;              ///< T - The metrics on this level have
                                      ///<     changed
  bool m_metricsDefined;              ///< T - The metrics have been defined for
                                      ///<     this level
  bool m_useScriptN;                  ///< For some mappings, it is not possible
                                      ///< to define scriptN.  In those cases,
                                      ///< N is calculated and averaged
                                      ///< directly.  The solution is
                                      ///< conservative but freestream-
                                      ///< preservation is lost across coarse-
                                      ///< fine interfaces.
                                      ///< T - Script N is used (default).
                                      ///< F - Script N is not used.  N is
                                      ///<     calculated directly and
                                      ///<     freestream preservation is lost
  bool m_isMultiblock;                ///< T - A multiblock coordinates system
  bool m_haveMultiBlockVectorData;    ///< T - Space vector data will be
                                      ///<     exchanged across multiblock
                                      ///<     boundaries

  static int s_verbosity;
};


/*******************************************************************************
 *
 * Class LevelGridMetrics: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Are the metrics defined on this level?
/** After preRegrid, this also indicates if a level is used.
 *  \return             T - yes.
 *//*-----------------------------------------------------------------*/

inline bool
LevelGridMetrics::metricsDefined() const
{
  return m_metricsDefined;
}

/*--------------------------------------------------------------------*/
//  Did the metrics change on the finer mesh?
/** \return             T - yes.
 *//*-----------------------------------------------------------------*/

inline bool
LevelGridMetrics::didFinerChange() const
{
  return m_finerLevelGridMetrics->m_metricsChanged;
}

/*--------------------------------------------------------------------*/
//  Get the coordinate system
/** \param[in]  box     Assigned to this box
 *//*-----------------------------------------------------------------*/

inline const NewCoordSys*
LevelGridMetrics::getCoordSys2(const Box& box) const
{
  return m_coordSys->getCoordSys(box);
}

inline const MultiBlockCoordSys&
LevelGridMetrics::getCoordSys() const
{
  return *m_coordSys;
}

#if 0
/*--------------------------------------------------------------------*/
//  Get the coordinate system
/** \param[in]  dit     Assigned to this data iterator
 *//*-----------------------------------------------------------------*/

inline const NewCoordSys*
LevelGridMetrics::getCoordSys2(const DataIterator& dit) const
{
  return m_coordSys->getCoordSys(m_grids[dit()]);
}
#endif

/*--------------------------------------------------------------------*/
//  Get the "fourth-order" coordinate system
/** \param[in]  box     Assigned to this box
 *//*-----------------------------------------------------------------*/

inline const NewFourthOrderCoordSys*
LevelGridMetrics::getCoordSys(const Box& box) const
{
  return static_cast<const NewFourthOrderCoordSys*>(m_coordSys->getCoordSys(box));
}

#if 0
/*--------------------------------------------------------------------*/
//  Get the "fourth-order" coordinate system
/** \param[in]  dit     Assigned to this data iterator
 *//*-----------------------------------------------------------------*/

inline const NewFourthOrderCoordSys*
LevelGridMetrics::getCoordSys(const DataIterator& dit) const
{
  return static_cast<const NewFourthOrderCoordSys*>(m_coordSys->getCoordSys(m_grids[dit()]));
}
#endif

/*--------------------------------------------------------------------*/
//  Get the block index
/** \param[in]  a_box   Assigned to this interior box
 *  \return             Block index
 *//*-----------------------------------------------------------------*/

inline int
LevelGridMetrics::getBlock(const Box& a_box) const
{
  return m_coordSys->whichBlock(a_box);
}

/*--------------------------------------------------------------------*/
//  Get a problem domain representing a block
/** \param[in]  a_box   Assigned to this interior box
 *  \param[in]  a_numGhost
 *                      The number of ghosts that the box representing
 *                      the block should be grown by before cropping.
 *                      Cropping *only* occurs at boundaries defined
 *                      by the physical domain
 *  \return             If singleblock, the problem domain.  If
 *                      multiblock, non-periodic problem domain
 *                      representing the block
 *  \note
 *  <ul>
 *    <li> If single block, the problem domain must be cached when
 *         the coordSys is defined
 *  </ul>
 *//*-----------------------------------------------------------------*/

inline const ProblemDomain&
LevelGridMetrics::blockDomain(const Box& a_box, const int a_numGhost)
{
  if (isMultiblock())
    {
      CH_assert(a_numGhost >= 0);
      const int iBlock = getBlock(a_box);
      CH_assert(iBlock >= 0);
      if (iBlock != m_cachedBlockIdx || a_numGhost != m_cachedBlockNG)
        {
          m_cachedBlockIdx = iBlock;
          m_cachedBlockNG  = a_numGhost;
          // Find the problem domain representing the block
          Box grownBlock =
            grow(m_coordSys->mappingBlocks()[iBlock], a_numGhost);
          if (a_numGhost > 0)
            {
              m_coordSys->keepInDomain(grownBlock, iBlock);
            }
          m_cachedBlockDomain.define(grownBlock);
        }
    }
  return m_cachedBlockDomain;
}

/*--------------------------------------------------------------------*/
//  Get the computational mesh spacing
/** \return             Mesh spacing in each direction
 *//*-----------------------------------------------------------------*/

inline const RealVect&
LevelGridMetrics::dxVect() const
{
  return m_dxVect;
}


/*====================================================================*/
/** \name Consistent access to structural boxes
 *
 *  There are several places to get the boxes and interators from but
 *  all should be the same.  These functions provide consistent
 *  access with sizes relevant to the metrics.
 *//*=================================================================*/
//@{

//  Get the boxes
inline const DisjointBoxLayout&
LevelGridMetrics::getBoxes() const
{
  return m_grids;
}

//  Get a level iterator
inline DataIterator
LevelGridMetrics::getDataIterator() const
{
  return getBoxes().dataIterator();
}

//  Resize a box to include the ghost cells (defines J)
inline Box
LevelGridMetrics::resizeWithGhosts(const Box& box) const
{
  return grow(box, m_ghostVect);
}

//  Determine the ghost vector required for \f$J\f$.
/** This is +1 over the ghost box to allow for derivatives needed
 *  to manipulate products
 */
inline IntVect
LevelGridMetrics::getJReqGhostVect() const
{
  return IntVect::Unit;
}

//  Determine the ghost vector required for \f$N\f$.
/** This is +2 over the ghost box to allow for tangential derivatives needed
 *  to compute the volume flux
 */
inline IntVect
LevelGridMetrics::getNReqGhostVect() const
{
  return 2*IntVect::Unit;
}

//  Resize a box to the size required for \f$J\f$
inline Box
LevelGridMetrics::resizeWithJReqGhosts(const Box& box) const
{
  return grow(box, getJReqGhostVect());
}

//  Resize a box to the size required for \f$N\f$
inline Box
LevelGridMetrics::resizeWithNReqGhosts(const Box& box) const
{
  return grow(box, getNReqGhostVect());
}
//@}


/*====================================================================*/
/** \name Easier access to AMRLevel data
 *//*=================================================================*/
//@{

inline int
LevelGridMetrics::level() const
{
  return m_parentAMRLevel->level();
}

inline int
LevelGridMetrics::refRatio() const
{
  return m_parentAMRLevel->refRatio();
}

inline const ProblemDomain&
LevelGridMetrics::problemDomain() const
{
  return m_parentAMRLevel->problemDomain();
}
//@}

/*--------------------------------------------------------------------*/
//  Set the coarsened-fine <U> used by the interpolator (CrFn access)
/** Grids must have been defined for the interpolator.  Returns a
 *  coarsening of 'this' level to be filled externally using a copyTo.
 *  'this' is the fine level.  I.e., 'this' is the fine level.
 *//*-----------------------------------------------------------------*/

inline LevelData<FArrayBox>&
LevelGridMetrics::presetCr2ThisInterpolatorCrFnLevU()
{
  CH_assert(hasCoarser());
  return m_coarserLevelGridMetrics->m_interpolator->presetCrFnLevU();
}

/*--------------------------------------------------------------------*/
//  Set the coarsened-fine <JU> used by the interpolator (CrFn access)
/** Grids must have been defined for the interpolator.  Returns a
 *  coarsening of 'this' level to be filled externally using a copyTo.
 *  'this' is the fine level.  I.e., 'this' is the fine level.
 *//*-----------------------------------------------------------------*/

inline LevelData<FArrayBox>&
LevelGridMetrics::presetCr2ThisInterpolatorCrFnLevJU()
{
  CH_assert(hasCoarser());
  return m_coarserLevelGridMetrics->m_interpolator->presetCrFnLevJU();
}

// Obsolete
// /*--------------------------------------------------------------------*/
// //  Get the coarsened-fine <U> used by the interpolator
// /** Grids must have been defined for the interpolator.  Returns a
//  *  coarsening of 'this' level.  I.e., 'this' is the fine level.
//  *//*-----------------------------------------------------------------*/

// inline const LevelData<FArrayBox>&
// LevelGridMetrics::getCr2ThisInterpolatorCrFnLevU() const
// {
//   CH_assert(hasCoarser());
//   return m_coarserLevelGridMetrics->m_interpolator.getCrFnLevU();
// }

// /*--------------------------------------------------------------------*/
// //  Get the coarsened-fine <JU> used by the interpolator
// /** Grids must have been defined for the interpolator.  Returns a
//  *  coarsening of 'this' level.  I.e., 'this' is the fine level.
//  *//*-----------------------------------------------------------------*/

// inline const LevelData<FArrayBox>&
// LevelGridMetrics::getCr2ThisInterpolatorCrFnLevJU() const
// {
//   CH_assert(hasCoarser());
//   return m_coarserLevelGridMetrics->m_interpolator.getCrFnLevJU();
// }

/*--------------------------------------------------------------------*/
//  Set the coarsened-fine <U> used by the interpolator
/** Grids must have been defined for the interpolator.  Copies the
 *  coarse data from 'this' level to the coarsened-fine representation
 *  in the interpolator.  'this' is the coarse level.
 *  \param[in]  a_CrU  \<U\> on the coarse level
 *  \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *//*-----------------------------------------------------------------*/

inline void
LevelGridMetrics::presetThis2FnInterpolatorCrFnLevU(
  const LevelData<FArrayBox>& a_CrLevU,
  const Interval&             a_vectorIntv)
{
  CH_assert(hasFiner());
  m_interpolator->presetCrFnLevU(a_CrLevU, a_vectorIntv);
}

// Should not be used
// /*--------------------------------------------------------------------*/
// //  Set the coarsened-fine <JU> used by the interpolator
// /** Grids must have been defined for the interpolator.  Copies the
//  *  coarse data from 'this' level to the coarsened-fine representation
//  *  in the interpolator.  'this' is the coarse level.
//  *//*-----------------------------------------------------------------*/

// inline void
// LevelGridMetrics::presetThis2FnInterpolatorCrFnLevJU(
//   const LevelData<FArrayBox>& a_CrLevJU)
// {
//   CH_assert(hasFiner());
//   m_interpolator.presetCrFnLevJU(a_CrLevJU);
// }

/*--------------------------------------------------------------------*/
//  Invalidate the CrFn data, both <U> and <JU>, used by the
//  interpolator
/** 'this' is the fine level.  If data goes out-of-date on the coarser
 *  level, then the interpolator, for interpolating to 'this' level,
 *  has out-of-date coarsened-fine data.  There must be a coarser
 *  level
 *//*-----------------------------------------------------------------*/

inline void
LevelGridMetrics::invalidateCr2ThisInterpolatorCrFnLevData()
{
  CH_assert(hasCoarser());
  m_coarserLevelGridMetrics->m_interpolator->invalidateCrFnLevData();
}

/*--------------------------------------------------------------------*/
//  Invalidate the CrFn data, both <U> and <JU>, used by the
//  interpolator
/** 'this' is the coarse level.  If data goes out-of-date on 'this'
 *  level, then the interpolator, for interpolating to the next finer
 *  level, has out-of-date coarsened-fine data.  You can invalidate
 *  even if there is no finer level.
 *//*-----------------------------------------------------------------*/

inline void
LevelGridMetrics::invalidateThis2FnInterpolatorCrFnLevData()
{
  m_interpolator->invalidateCrFnLevData();
}

/*--------------------------------------------------------------------*/
//  Number of ghosts used to build the coarsened-fine data in the
//  interpolator
/** \param[in]  a_force T - avoid assertion checking.
 *                      F - (default)
 *  \note
 *  <ul>
 *    <li> The caller should ensure that a finer level exists if the
 *         assertion is disabled.  Commonly, the return value is
 *         required inbetween initialGrid() and postInitialGrid()
 *         where the metrics have not yet been defined (so the
 *         assertion would fail) but the grid has been defined (so
 *         this value is known).
 *  </ul>
 *//*-----------------------------------------------------------------*/

inline const IntVect&
LevelGridMetrics::interpolatorCrFnNumGhost(const bool a_force) const
{
  if (!a_force)
    {
      CH_assert(hasFiner());
    }
  return m_interpolator->CrFnNumGhost();
}

/*--------------------------------------------------------------------*/
//  Is this a multiblock grid?
/** \return             T - Yes
 *//*-----------------------------------------------------------------*/

inline bool
LevelGridMetrics::isMultiblock() const
{
  return m_isMultiblock;
}

#include "NamespaceFooter.H"
#endif
