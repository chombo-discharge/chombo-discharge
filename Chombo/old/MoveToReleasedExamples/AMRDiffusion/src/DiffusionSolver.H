#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIFFUSIONSOLVER_H_
#define _DIFFUSIONSOLVER_H_

#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BiCGStabSolver.H"

#include "BackwardEuler.H"

#include "DiffusionParams.H"

#include "UsingNamespace.H"

/// A solver for the diffusion equation with source and sink terms
class DiffusionSolver
{
  public:
    /// Construct with the solver parameters
    DiffusionSolver(const DiffusionParams & a_params);

    /// Destructor
    ~DiffusionSolver();

    /// Initialize the domain, geometry, and solver
    void init();

    /// Run the solver to completion
    void run();

  protected:
    // Initialize the index space
    void initIndexSpace();

    // Initialize the data
    void initData();

    // Initialize the solver
    void initSolver();

    // Set the current source term
    void setSource();

    // Write an output file
    void writeOutput(const int  & a_step,
                     const Real & a_time);

    // All the solver parameters
    DiffusionParams m_params;

    // All the grids at all AMR levels
    Vector<DisjointBoxLayout> m_grids;

    // The solution at the old and new times
    Vector<LevelData<FArrayBox>* > m_solnOld;
    Vector<LevelData<FArrayBox>* > m_solnNew;
    Vector<LevelData<FArrayBox>* > m_source;

    // The backward Euler time integrator
    BackwardEuler* m_integrator;

    // The bottom solver used by multigrid when it can no longer coarsen
    BiCGStabSolver<LevelData<FArrayBox> > m_bottomSolver;

  private:
};

#endif
