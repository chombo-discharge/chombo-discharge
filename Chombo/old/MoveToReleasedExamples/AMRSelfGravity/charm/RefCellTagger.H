#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REFCELLTAGGER_H_
#define _REFCELLTAGGER_H_

#include <iostream>
#include <string>

#include "LevelData.H"
#include "FArrayBox.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "REAL.H"
#include "RealVect.H"

///
/**
   Class that tag cells that need to be refined based on various criteria:
     a_useRefineShocks       - whether shocks         ...
     a_useRefineJeansLength -   ...   Jeans length   ...
     a_useRefineOverDensity -   ...   overdensity    ...
     a_useRefineGradient    -   ...   gradient  should be refined
   Theresholds:
     a_presJumpThreshold
     a_cellMassThreshold
     a_jeansResolThreshold
     a_gradientThreshold
   Parameters:
     a_gradVarInterval     - interval for which to compute the grad
 */

// define how to enforce refinement given a member refine_box: refine
// this whole box alone (ONLY), refine within this box alone when any
// other criterion is fulfilled (AND), refine within this box in
// addition to any other cells (OR)

enum RefineMode
{
  FIX, AND, OR
};

class RefCellTagger
{
public:
  /// Null constructor
  RefCellTagger();

  /// Constructor which defines no parameters
  RefCellTagger(const ProblemDomain& a_problemDomain,
                const Real&          a_dx,
                const int&           a_refRatio);

  /// Destructor
  virtual ~RefCellTagger()
  {
  }

  ///
  void define();

  ///
  void define(const ProblemDomain& a_problemDomain,
              const Real&          a_dx,
              const int&           a_refRatio);
  ///
  void setParameters(const bool           a_useRefineGradient,
                     const Real&          a_gradientThreshold,
                     const Interval&      a_gradVarInterval,
                     const bool           a_useRefineOverDensity,
                     const Real&          a_cellMassThreshold,
                     const bool           a_useRefineShocks,
                     const Real&          a_presJumpThreshold,
                     const bool           a_useRefineVorticity,
                     const Real&          a_vorticityThreshold,
                     const bool           a_useRefineJeansLength,
                     const Real&          a_jeansResolThreshold,
                     const bool           a_useRefineRegion,
                     const Vector<Box>&   a_refineRegion,
                     const Vector<RefineMode>& a_refineMode);

  /// Factory method - this object is its own factory
  virtual RefCellTagger* new_refCellTagger() const;

  /// set the pressure jump threshold for refinement
  void setRefineShocks(const bool  a_useRefineShocks,
                       const Real& a_presJumpThreshold);

  /// set the pressure jump threshold for refinement
  void setRefineVorticity(const bool  a_useRefineVorticity,
                          const Real& a_vortexThreshold);

  /// set the variables for gradient
  void setRefineGradient(const bool      a_useRefineGradient,
                         const Real&     a_gradientThreshold,
                         const Interval& a_gradVarInterval);

  /// set the cell mass threshold for refinement
  void setRefineOverdense(const bool  a_useRefineOverDensity,
                          const Real& a_cellMassThreshold);

  /// set the jeans resol. threshold for refinement
  void setRefineJeans(const bool  a_useRefineJeans,
                      const Real& a_jeansResolThreshold);

  /// setup for refinement of a region: a mode can be union or intersection for
  void setRefineRegion(const bool                a_useRefineRegion,
                       const Vector<Box>&        a_refineRegion,
                       const Vector<RefineMode>& a_refMode);

  ///
  bool refineShocks()
  {
    return m_useRefineShocks;
  }

  ///
  bool refineVorticity()
  {
    return m_useRefineVorticity;
  }

  ///
  bool refineOverdense()
  {
    return m_useRefineOverDensity;
  }

  ///
  bool refineGradient(const std::string a_var);

  ///
  bool refineJeans()
  {
    return m_useRefineJeansLength;
  }

  ///
  bool refineRegion()
  {
    return m_useRefineRegion;
  }

  /// whether or not need to interopolate/copy ghost cells
  virtual bool needGhostCells()
  {
    return (m_useRefineGradient || m_useRefineShocks || m_useRefineVorticity);
  }

  ///
  virtual IntVectSet tagShocks(const LevelData<FArrayBox>& a_U);

  ///
  virtual IntVectSet tagVorticity(const LevelData<FArrayBox>& a_U);

  ///
  virtual IntVectSet tagJeans(const LevelData<FArrayBox>& a_U);

  ///
  virtual IntVectSet tagOverdense(const LevelData<FArrayBox>& a_rho);

  ///
  virtual IntVectSet tagGradient(const LevelData<FArrayBox>& a_U);

  /// Tag all cells in the specified region of space
  virtual IntVectSet tagRegion(const RealVect&          a_low,
                               const RealVect&          a_high,
                               const DisjointBoxLayout& a_grids);

  /// Tag all cells in the region defined in the setRefinedRegion function.
  virtual IntVectSet tagRegion(const IntVectSet& a_tags,
                               const int         a_level,
                               const DisjointBoxLayout& a_grids);

protected:
  //
  void relativeGradient(FArrayBox&           a_gradMagFab,
                        const FArrayBox&     a_varFab,
                        const int            a_comp,
                        const Box&           a_box);

protected:
  // has define been called, so that dx, problemDomain and refRatio have been defined
  bool m_isDefined;

  // grid mesh size
  Real m_dx;

  // grid mesh size
  ProblemDomain m_problemDomain;

  // refinement ratio
  Real m_refRatio;

  // pressure jump threshold for refinement, use this criter., is all defined
  Real m_presJumpThreshold;
  bool m_useRefineShocks;
  bool m_isRefineShocksSet;

  // vorticity threshold for refinement, use this criter., is all defined
  Real m_vorticityThreshold;
  bool m_useRefineVorticity;
  bool m_isRefineVorticitySet;

  // cell density threshold for refinement, use this criter., is all defined
  Real m_cellMassThreshold;
  bool m_useRefineOverDensity;
  bool m_isRefineOverDensitySet;

  // Jeans resol. threshold for refinement, use this criter., is all defined
  Real m_jeansResolThreshold;
  bool m_useRefineJeansLength;
  bool m_isRefineJeansLengthSet;

  // gradient threshold for refinement, use this criter., is all defined
  Real m_gradientThreshold;
  bool m_useRefineGradient;
  Interval m_gradVarInterval;
  bool m_isRefineGradientSet;

  // refine within this box according to a specified refinement mode;
  // this will be applied internally
  Vector<Box>     m_refineRegion;
  Vector<RefineMode> m_refineMode;
  bool m_useRefineRegion;
  bool m_isRefineRegionSet;
};

#endif
