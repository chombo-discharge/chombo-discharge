#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEHEATSOLVER_H_
#define _BASEHEATSOLVER_H_

#include <Box.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <LevelFluxRegister.H>
#include <FArrayBox.H>
#include <ProblemDomain.H>
#include "BaseHelmholtzOp.H"
#include "NamespaceHeader.H"

/// Pure virtual base class to encapsulate heat equation solution.
/** The baseHeatSolver class is a virtual base class for implementing
    different approaches to solving the heat equation using BaseHelmholtzOp
    derived operators (note that we use BaseHelmholtzOp as a base class rather
    than LevelOp),
*/

class BaseHeatSolver
{

public:

  /// default constructor
  BaseHeatSolver()
  {
  }

  /// destructor
  virtual ~BaseHeatSolver()
  {
  }

  /// complete constuctor
  /** Note that the BaseHelmholtzOp passed in must already have the
      diffusion coefficient (nu) set.
  */
  BaseHeatSolver(const DisjointBoxLayout& a_levelGrids,
                 const DisjointBoxLayout* a_crseGrids,
                 const ProblemDomain& a_domain,
                 Real a_dxLevel,
                 int a_nRefCrse,
                 const BaseHelmholtzOp* a_helmop,
                 int a_ncomp = 1)
  {
    MayDay::Error("BaseHeatSolver is a virtual base class!");
  }

  /// define statement
  /** Note that the BaseHelmholtzOp passed in must already have the
      diffusion coefficient (nu) set.
  */
  virtual void define(const DisjointBoxLayout& a_levelGrids,
                      const DisjointBoxLayout* a_crseGrids,
                      const ProblemDomain& a_domain,
                      Real a_dxLevel,
                      int a_nRefCrse,
                      const BaseHelmholtzOp* a_helmop,
                      int ncomp = 1) =0;

  /// do solve
  /** update phi according to  d(phi)/dt = L(phi) + src
   */
  virtual void updateSoln(LevelData<FArrayBox>& a_phiNew,
                          LevelData<FArrayBox>& a_phiOld,
                          LevelFluxRegister* a_fineFluxRegPtr,
                          LevelFluxRegister* a_crseFluxRegPtr,
                          const LevelData<FArrayBox>* a_crsePhiOldPtr,
                          Real crseOldTime,
                          const LevelData<FArrayBox>* a_crsePhiNewPtr,
                          Real crseNewTime,
                          // note that a_src can't be const because of BCs
                          LevelData<FArrayBox>& a_src,
                          Real oldTime,
                          Real dt) = 0;

  /// computes time-centered diffusive term for explicit updating
  /** compute time-centered L(phi) for use in subsequent update
      operations
  */
  virtual void computeDiffusion(LevelData<FArrayBox>& a_DiffusiveTerm,
                                LevelData<FArrayBox>& a_phiOld,
                                LevelFluxRegister* a_fineFluxRegPtr,
                                LevelFluxRegister* a_crseFluxRegPtr,
                                const LevelData<FArrayBox>* a_crsePhiOldPtr,
                                Real crseOldTime,
                                const LevelData<FArrayBox>* a_crsePhiNewPtr,
                                Real crseNewTime,
                                LevelData<FArrayBox>& a_src,
                                Real oldTime,
                                Real dt) = 0;

  /// "virtual constructor"
  virtual BaseHeatSolver* new_heatSolver() const = 0;

  /// set solver tolerance
  virtual void setSolverTolerance(Real a_solverTol) = 0;

  /// accessor
  virtual Real solverTolerance() const = 0;

};

#include "NamespaceFooter.H"
#endif
