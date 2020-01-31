#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 10 Apr 2008

#ifndef _MULTIBLOCKLEVELDATA_H_
#define _MULTIBLOCKLEVELDATA_H_

#include "MappedDomain.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "IVSFAB.H"

#include "NamespaceHeader.H"

/// Class to hold data on a mapped multiblock domain.
class MultiBlockLevelData
{
public:

  /**
     \name Constructor, destructor, and defines.
  */
  /*@{*/

  /// null constructor
  MultiBlockLevelData();

  /// null constructor
  virtual ~MultiBlockLevelData();

  ///
  /** Constructor.
   */
  MultiBlockLevelData(/// multiblocked mapped domain
                      const MappedDomain&        a_mappedDomain,
                      /// CELL-centered grids
                      const DisjointBoxLayout&   a_grids,
                      /// ghost width
                      const IntVect&             a_ghost = IntVect::Zero,
                      /// radius of neighborhood from which to interpolate
                      int                        a_radius = 1,
                      /// max degree of moments of displacement
                      int                        a_degree = 1);


  ///
  /** Define MultiBlockLevelData object.
  */
  virtual void define(/// multiblocked mapped domain
                      const MappedDomain&        a_mappedDomain,
                      /// CELL-centered grids
                      const DisjointBoxLayout&   a_grids,
                      /// ghost width
                      const IntVect&             a_ghost = IntVect::Zero,
                      /// radius of neighborhood from which to interpolate
                      int                        a_radius = 1,
                      /// max degree of moments of displacement
                      int                        a_degree = 1);

  /*@}*/
  /**
     \name Settings
  */
  /*@{*/

  ///
  /** Set verbosity.  Default is 0.
   */
  virtual void setVerbose(int   a_verbose);

  ///
  /** Set order, 2 or 4.  Default is 2.
   */
  void setOrder(int   a_order);

  ///
  /** Set factor by which to multiply m_ghost to get m_ghostMore.
      Normally radius + 1.
   */
  void setGhostFactor(int   a_ghostFactor);

  ///
  /** Set whether to use averages in calculation.
   */
  void setUseAverage(bool   a_useAverage);

  ///
  /** Set whether to find averages for ghosts.
   */
  void setAverageGhost(bool   a_averageGhost);

  /*@}*/
  /**
     \name Functions for integrals
  */
  /*@{*/

  ///
  /** set up weights and points, Gauss quadrature
   */
  virtual void setGaussQuadrature(int   a_integralLength);

  ///
  /** set up weights and points, Newton-Cotes quadrature
   */
  virtual void setNewtonCotesQuadrature(int   a_integralLength);

  ///
  /** integrate
   */
  virtual void integrateOnCells(/// CELL-centered integral on each cell, predefined and filled in by this function
                                BoxLayoutData<FArrayBox>&         a_integral,
                                /// function values on all fine integration points
                                const BoxLayoutData<FArrayBox>&   a_function);

  /*@}*/
  /**
     \name Interpolations
  */
  /*@{*/

  ///
  /** Using average function values, interpolate function values on centers
      of valid cells, and store in a_funCenter.

      If m_order == 2 then merely copy a_funAvg to a_funCenter on valid cells.
      If m_order == 4 then get some derivatives and set a_funCenter on
      valid cells interior to each block.  (Valid cells next to block boundaries
      stay undefined, to be filled in by fillAvgGhosts().)
   */
  virtual void interpCenterFromAvg(/// function values on CELL centers, predefined and filled in by this function
                                   LevelData<FArrayBox>&         a_funCenter,
                                   /// averaged function values over CELLs
                                   const LevelData<FArrayBox>&   a_funAvg);

  /*@}*/
  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Retrieve cell centers.
   */
  const LevelData<FArrayBox>& cellCenters();

  ///
  /** Retrieve integration points.
   */
  const LevelData<FArrayBox>& physIntegrationPointsValidAll();

  ///
  /** Retrieve integration points.
   */
  const LevelData<FArrayBox>& physIntegrationPointsGhostedAll();

  ///
  /** Retrieve integration points.
   */
  const LevelData<FArrayBox>& cellVolumesValid();

  ///
  /** Retrieve integration points.
   */
  const LevelData<FArrayBox>& cellVolumesGhosted();

  ///
  /** Retrieve block number.
   */
  int block(const DataIterator&   a_dit) const
  {
    return m_block[a_dit];
  }

  ///
  /** Retrieve block number.
   */
  int block(const DataIndex&   a_dind) const
  {
    return m_block[a_dind];
  }

  /*@}*/
  /**
     \name Exchange functions
  */
  /*@{*/

  ///
  /** Exchange a_calcCenter and *a_calcAvgPtr, on cell centers or averages,
      depending on m_useAverage.
   */
  virtual void exchangeAll(LevelData<FArrayBox>&   a_calcAvg,
                           LevelData<FArrayBox>*   a_condPtr = NULL,
                           LevelData<FArrayBox>*   a_calcCenterPtr = NULL);

  ///
  /** If m_useAverage:

      Entering this function, a_calcCenter holds function values
      on centers of valid cells and internal ghost cells of blocks;
      on ghost cells of blocks, it holds the same obtained from
      exchangeCopyOnly().

      This function then interpolates to obtain function values
      on centers of ghost cells of blocks.
      At each ghost cell, solve a least squares problem,
      and store the condition number in *a_condPtr.

      If a_calcAvgPtr is non-NULL, then also fill in interpolated
      averages on ghost cells of blocks.

      If !m_useAverage:

      Entering this function, a_calcAvg holds average function values
      over valid cells and over internal ghost cells of blocks;
      on ghost cells of blocks, it holds the same obtained from
      exchangeCopyOnly().

      This function then interpolates to obtain average function values
      over ghost cells of blocks.
      At each ghost cell, solve a least squares problem,
      and store the condition number in *a_condPtr.

      If a_calcCenterPtr is non-NULL, then also fill in interpolated
      averages over ghost cells of blocks;
      and if also m_order == 4 then do this not only on ghost cells of blocks,
      but also on valid cells next to block boundaries.
      The interior cell centers should have been filled in already
      by interpCenterFromAvg.
   */
  virtual void fillGhosts(/// function values at cell centers
                          LevelData<FArrayBox>&   a_calcCenter,
                          /// condition numbers from least squares
                          LevelData<FArrayBox>*   a_condPtr,
                          /// interpolated averages on ghost cells (optional)
                          LevelData<FArrayBox>*   a_calcAvgPtr);

  ///
  /** Set moments for filling ghosts.
   */
  virtual void setMoments(Vector<Real>&      a_moments,
                          Vector<Real>&      a_ghostMoments,
                          Vector<IntVect>&   a_neighborhood,
                          int                a_block,
                          const IntVect&     a_ghostCell,
                          const RealVect&    a_ghostCellCenter,
                          const FArrayBox&   a_cellCentersValidFab,
                          const FArrayBox&   a_cellVolumesValidFab,
                          const FArrayBox*   a_cellVolumesGhostedFabPtr,
                          const FArrayBox&   a_weightedJacobianValidFab,
                          const FArrayBox*   a_weightedJacobianGhostedFabPtr,
                          FArrayBox&         a_physIntegrationPointsFab);

  /*@}*/

protected:

  ///
  /** Exchange given components of a_data, including through topological
      neighbors given in MappedDomain.
   */
  virtual void exchangeCopyOnly(LevelData<FArrayBox>&   a_data,
                                const Interval&         a_intvl,
                                bool                    a_refine = false);

  ///
  /** Exchange all components of a_data, including through topological
      neighbors given in MappedDomain.
   */
  virtual void exchangeCopyOnly(LevelData<FArrayBox>&   a_data,
                                bool                    a_refine = false);


  ///
  /** Set m_physIntegrationPointsValidAll and m_weightedJacobianAll.
   */
  virtual void setIntegralPointsAndWeights(/// cartesian coordinates within unit cell
                                           const FArrayBox&   a_cartesianOriginCoordsFab,
                                           /// weights of points
                                           const FArrayBox&   a_integWeightsFab);

  ///
  /** Set gradients of Jacobian at cell centers.
   */
  virtual void setGradJacobianCellCenters();

  ///
  /** Compute cell volumes, and store in m_cellVolumes.
      Called by setIntegralPointsAndWeights().
   */
  void setCellVolumes();

  ///
  /** Fin all moments.
   */
  void setAllMoments();

  ///
  /** Returns whether a_iv is in the forbidden zone of a_bx:
      that is, more than one component of a_iv is outside a_bx.
   */
  virtual bool forbiddenZone(const IntVect&   a_iv,
                             const Box&       a_bx);

  /** verbosity level
   */
  int m_verbose;

  /** Do we have m_mappedDomain?
   */
  bool m_isDefined;

  /** grids on which data lie
   */
  DisjointBoxLayout m_grids;

  /** m_grids grown by m_ghost except on external boundaries
   */
  //   BoxLayout m_gridsGhostedInternal;

  /** ghost layer of m_data
   */
  IntVect m_ghost;

  /** expanded ghost layer needed in topological exchange
   */
  IntVect m_ghostMore;

  /** radius of neighborhood from which to interpolate
   */
  int m_radius;

  /** max degree of moments of displacement
   */
  int m_degree;

  /** size of neighborhood with radius m_radius:
      (2*radius + 1)^SpaceDim
   */
  int m_neqns;

  /** number of coefficients up to degree m_degree:
      C(m_degree + SpaceDim, SpaceDim)
      = (m_degree + SpaceDim)! / (SpaceDim! * m_degree!)
   */
  int m_nvars;

  /** order of approximation:  2 or 4
   */
  int m_order;

  /** factor by which to multiply m_ghost to get m_ghostMore;
      normally m_radius + 1.
   */
  int m_ghostFactor;

  /** degrees from 0 to m_degree in each dimension
   */
  Box m_degreeBox;

  /** length of m_workArray
   */
  int m_lwork;

  /** used for least squares
   */
  Real* m_workArray;

  /** block in which each grid is located
   */
  LayoutData<int> m_block;

  /** SpaceDim components of physical coordinates of cell centers
      on m_grids, with ghost layer m_ghost.
      At ghost cells, copy from box with valid cells.
   */
  LevelData<FArrayBox> m_cellCentersValid;

  /** SpaceDim components of physical coordinates of cell centers
      on m_grids, with ghost layer m_ghost.
      Do not do exchange.  If m_ghost == IntVect::Zero then this is unused.
   */
  LevelData<FArrayBox> m_cellCentersGhosted;

  /** Volumes of cells.  On m_grids, 1 component.
      We'll call exchangeCopyOnly() on it.
   */
  LevelData<FArrayBox> m_cellVolumesValid;

  /** Volumes of cells, used iff m_averageGhost.  On m_grids, 1 component.
      We will NOT call exchangeCopyOnly() on it.
   */
  LevelData<FArrayBox> m_cellVolumesGhosted;

  /** multiblock mapped domain
   */
  MappedDomain m_mappedDomain;

  /** multiblock mapped domain refined for quadrature
   */
  MappedDomain m_mappedDomainRefinedForQuadrature;

  /** ghost cells of m_grids
   */
  LayoutData< IntVectSet > m_ghostBetweenBlockCells;

  /** number of points in neighborhood of ghost cell; 1 component
   */
  LayoutData< IVSFAB<int>* > m_neighborhoodSize;

  /** moments; m_neqns * m_nvars components
   */
  LayoutData< IVSFAB<Real>* > m_moments;

  /** ghost moments, if m_averageGhost; m_nvars components
   */
  LayoutData< IVSFAB<Real>* > m_ghostMoments;

  /** indices of neighborhood of each ghost cell; m_neqns components
   */
  LayoutData< IVSFAB<IntVect>* > m_neighborhood;

  /** gradients of Jacobians at cell centers; SpaceDim components
   */
  LevelData<FArrayBox> m_gradJacobianCellCenters;

  /** have we set gradients of Jacobians at cell centers?
   */
  bool m_gotGradJacobianCellCenters;

  /** have we set integration points and weights?
   */
  bool m_gotIntegralStuff;

  /** lengths of integrals (in each dimension) in each cell
   */
  int m_integralLength;

  /** CELL-centered grids on which Gauss quadrature points lie
   */
  DisjointBoxLayout m_gridsRefinedForQuadrature;

  /** grids on which integration points lie
   */
  BoxLayout m_gridsIntegrationAll;

  /** which quadrature method?
   */
  enum QuadratureMethod
  {
    NEWTONCOTES,
    GAUSS
  };

  QuadratureMethod m_quadratureMethod;

  /** integration points on each cell
   */
  Box m_integralCellBox;

  /** whether to use averages over cells
   */
  bool m_useAverage;

  /** whether to find averages over ghost cells
   */
  bool m_averageGhost;

  /** on m_integralCellBox, with SpaceDim*(m_degree + 1) components;
      used for storing powers of displacement
   */
  FArrayBox m_physPowFab;

  /** integration points, on m_gridsIntegrationAll, SpaceDim components,
      with ghost layer m_integralLength*m_ghostMore
   */
  LevelData<FArrayBox> m_physIntegrationPointsValidAll;

  /** integration points, on m_gridsIntegrationAll, SpaceDim components,
      with ghost layer m_integralLength*m_ghostMore.
      NO exchange.
   */
  LevelData<FArrayBox> m_physIntegrationPointsGhostedAll;

  /** integration weights, on m_gridsRefinedForQuadrature,
      with ghost layer m_integralLength*m_ghostMore
   */
  LevelData<FArrayBox> m_weightedJacobianValidAll;

  /** integration weights, on m_gridsRefinedForQuadrature,
      with layer m_integralLength*m_ghostMore
   */
  LevelData<FArrayBox> m_weightedJacobianGhostedAll;
};

#include "NamespaceFooter.H"

#endif // include guard
