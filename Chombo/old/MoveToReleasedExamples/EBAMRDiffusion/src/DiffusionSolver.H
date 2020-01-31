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

#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBBackwardEuler.H"

#include "DiffusionParams.H"
#include "EBConductivityOpFactory.H"

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

    void getEBLGAndQuadCFI(Vector<EBLevelGrid>                   & a_ebLevelGrids,
                           Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadCFInterp);

    void getVariableCoeffOpFactory(RefCountedPtr<EBConductivityOpFactory>& a_factory);

    void defineVariableCoeffs(Vector<RefCountedPtr<LevelData<EBCellFAB> > >        & a_aco,
                              Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        & a_bco,
                              Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > & a_bcoIrreg);

    void setConductivityCoefs(LevelData<EBCellFAB>        & a_aco,
                              LevelData<EBFluxFAB>        & a_bco,
                              LevelData<BaseIVFAB<Real> > & a_bcoIrreg,
                              const Real                  & a_dx);


    // Initialize the geometry
    void initGeometry();

    // Initialize the index space
    void initIndexSpace();

    // Initialize the data
    void initData();

    // Initialize the solver
    void initSolver();

    // Initialize the solver in the case of constant coeff
    void initSolverConstantCoeff();

    // Initialize the solver in the case of variable coeff
    void initSolverVariableCoeff();

    // Set the current source term
    void setSource();

    // Write an output file
    void writeOutput(const int  & a_step,
                     const Real & a_time);

    // All the solver parameters
    DiffusionParams m_params;

    // All the grids at all AMR levels
    Vector<DisjointBoxLayout> m_grids;

    // EBISLayouts for all AMR levles
    Vector<EBISLayout> m_ebisl;

    // The solution at the old and new times
    Vector<LevelData<EBCellFAB>* > m_solnOld;
    Vector<LevelData<EBCellFAB>* > m_solnNew;

    // The source/sink terms (combined into one "source")
    Vector<LevelData<EBCellFAB>* > m_source;

    // The backward Euler time integrator
    EBBackwardEuler* m_integrator;

    // Current time
    Real m_time;

    // The bottom solver used by multigrid when it can no longer coarsen
    BiCGStabSolver<LevelData<EBCellFAB> > m_bottomSolver;

  private:
};

#endif
