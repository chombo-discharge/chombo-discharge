#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELHELMHOLTZSOLVER_H_
#define _LEVELHELMHOLTZSOLVER_H_

#include "LevelSolver.H"
#include "BaseHelmholtzOp.H"
#include "GhostBC.H"
#include "NamespaceHeader.H"

/// special case LevelSolver to solve Helmholtz eqn
/** LevelHelmholtzSolver is a special case of LevelSolver
    for solving the Helmholtz equation on a level.  It
    allows the coefficients to be changed without allocating
    an entire new LevelSolver. */
class LevelHelmholtzSolver: public LevelSolver
{

public:
  /// Default constructor
  LevelHelmholtzSolver();

  virtual ~LevelHelmholtzSolver();

  ///
  LevelHelmholtzSolver(const DisjointBoxLayout& a_grids,
                       const DisjointBoxLayout* a_baseGrids,
                       const ProblemDomain& a_domain,
                       Real a_dxLevel,
                       int a_nRefCrse,
                       const BaseHelmholtzOp* const a_opin,
                       int a_ncomp = 1);

  /// deprecated interface using Box instead of ProblemDomain
  LevelHelmholtzSolver(const DisjointBoxLayout& a_grids,
                       const DisjointBoxLayout* a_baseGrids,
                       const Box& a_domain,
                       Real a_dxLevel,
                       int a_nRefCrse,
                       const BaseHelmholtzOp* const a_opin,
                       int a_ncomp = 1);

  ///
  virtual void define(const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_baseGrids,
                      const ProblemDomain& a_domain,
                      Real a_dxLevel,
                      int a_nRefCrse,
                      const BaseHelmholtzOp* const a_opin,
                      int a_ncomp = 1);


  /// (deprecated) interface using Box instead of ProblemDomain
  virtual void define(const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_baseGrids,
                      const Box& a_domain,
                      Real a_dxLevel,
                      int a_nRefCrse,
                      const BaseHelmholtzOp* const a_opin,
                      int a_ncomp = 1);


  /// this function hides the LevelSolver version -- should never be called
  virtual void define(const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_baseGrids,
                      const ProblemDomain& a_domain,
                      Real a_dxLevel,
                      int a_nRefCrse,
                      const LevelOp* const a_opin,
                      int a_ncomp = 1);

  /// this function hides the LevelSolver version -- should never be called
  virtual void define(const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_baseGrids,
                      const Box& a_domain,
                      Real a_dxLevel,
                      int a_nRefCrse,
                      const LevelOp* const a_opin,
                      int a_ncomp = 1);



  /// multiply beta by a_scale
  void scaleBeta(Real a_scale);

  ///
  void setBC(const DomainGhostBC& a_BC);



};

#include "NamespaceFooter.H"
#endif


