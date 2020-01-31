#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Tues, July 6, 1999

#ifndef _AMRLEVELMG_H_
#define _AMRLEVELMG_H_

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelFluxRegister.H"
#include "LevelMG.H"
#include "MGInterp.H"
#include "Copier.H"
#include "NamespaceHeader.H"

// forward declaration of AMRSolver
class AMRSolver;

/// Class which manages data and operations on a level for AMRSolver
/**
   AMRLevelMG manages the operations on a level for AMRSolver.
   AMRSolver contains a vector of AMRLevelMG's (one for each level).
   This class should be considered internal to AMRSolver and should
   not be considered part of the Chombo API.
 */
class AMRLevelMG
{
public:
  friend class AMRSolver;

  /// Null constructor
  AMRLevelMG();

  /// Destructor
  ~AMRLevelMG();

  /// Define AMRLevelMG from the inputs.
  AMRLevelMG(const AMRSolver* const a_parent,
             int                    a_level,
             const LevelOp* const   a_opin,
             int                    a_ncomp=1);

  /// Define AMRLevelMG from the inputs.
  /**
     The AMRLevelMG will only work for the
     input AMRSolver.  A LevelMG is created.
  */
  void define(const AMRSolver* const a_parent,
              int                    a_level,
              const LevelOp* const   a_opin,
              int                    a_ncomp=1);

  /// Set number of multigrid smoothings on way up V-cycle;  Default == 4.
  void setNumSmoothUp(int a_numSmoothUp);

  ///  Set number of multigrid smoothings on way down V-cycle;  Default == 4.
  void setNumSmoothDown(int a_numSmoothDown);

  /// Complete, 3-level operator.
  /**
      Applies multilevel operator with
      inhomogeneous coarse-fine boundary conditions (and
      inhomogeneous domain boundary conditions.  Refluxing
      is used to compute the operator at coarse-fine boundaries.
  */
  void applyAMROperator(Vector<LevelData<FArrayBox> *>& a_phiLevel,
                        LevelData<FArrayBox>&           a_Lofphi);

  /// Complete, 3-level operator with homogeneous physicalboundary conditions.
  /** Applies multilevel operator with inhomogeneous
      physical boundary conditions.  Refluxing is used to compute
      the operator at coarse-fine boundaries.
  */
  void applyAMROperatorHphys(Vector<LevelData<FArrayBox> *>& a_phiLevel,
                             LevelData<FArrayBox>&           a_LofPhi);

  /// Compute complete, 3-level residual and put it into local data m_resid
  void computeAMRResidual(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                          const Vector<LevelData<FArrayBox> *>& a_rhsLevel);

  /// Compute complete, 3-level residual with homogeneous physical BCs
  /**
      Compute complete, 3-level residual with homogeneous physical
      BCs and inhomogeneous C/F BCs and put it into local data m_resid
  */
  void computeAMRResidualHphys(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                               const Vector<LevelData<FArrayBox> *>& a_rhsLevel);

  /// Compute complete, 3-level residual and return it in a_resid
  void computeAMRResidual(LevelData<FArrayBox>&                 a_resid,
                          Vector<LevelData<FArrayBox>* >&       a_phiLevel,
                          const Vector<LevelData<FArrayBox>* >& a_rhsLevel);

  ///  Sweep up multigrid v-cycle.
  void upSweep(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
               const Vector<LevelData<FArrayBox> *>& a_rhsLevel);

  ///  Sweep down v-cycle
  void downSweep(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                 const Vector<LevelData<FArrayBox> *>& a_rhsLevel);

  /// Compute and return norm of internal data m_resid.
  Vector<Real> computeResidualNorm(int normType) const;

  /// Return pointer to the LevelOp used to define this AMRLevelMG
  LevelOp* levelOpPtr() const;

  /// Returns true if this object has been defined.
  bool isDefined() const;

  /// Returns normType norm of mfab
  Vector<Real> computeNorm(const LevelData<FArrayBox>& a_mfinput,
                           int                         a_normType) const;

  /// sets convergences metrics (if appropriate) for this problem.
  /** If we're using this class in the context of a solve, we may
      wish to set a metric against which to definte convergence (for
      example, the norm of the RHS).  As currently implemented, this
      is simply a pass-through to the LevelOp classes from this class.
  */
  void setConvergenceMetric(Real a_metric, int a_comp);

protected:
  ///
  /** Use short v-cycle or levelmg's smoother to reduce high wave
      numbers of residual
  */
  void smooth(LevelData<FArrayBox>&       a_soln,
              const LevelData<FArrayBox>& a_rhs);

  ///
  /**
      Delete internal memory and set pointers to null
  */
  void clearMemory();

  ///
  /**
      Set default values.  does not deal with memory
  */
  void setDefaultValues();

  ///
  /**
     Reflux enforces flux-matching from finer levels.
     Steps involved:
     - a) Initialize flux register with coarse flux
     - b) Increment flux register with sum of fine fluxes
     - b) Modify solution with flux diffrence
  */
  void reflux(Vector<LevelData<FArrayBox> *>& a_phiLevel,
              LevelData<FArrayBox>&           a_Lofphi);

  ///
  /**
     Fill fluxregisters with coarse flux
  */
  void initFRCoarse(Vector<LevelData<FArrayBox> *>& a_phiLevel);

  ///
  /**
     Increment fluxregisters with fine flux
  */
  void incrementFRFine(Vector<LevelData<FArrayBox> *>& a_phiLevel);

  ///
  /**
     The parent AMRSolver which manages the overall multilevel solve
  */
  const AMRSolver* m_parent;

  ///
  /**
     LevelMG level solver object to relax on this level
  */
  LevelMG m_levelMG;

  ///
  /**
     For interpolation
  */
  MGInterp m_mginterp;

  ///
  /**
      Flux register with next coarser level
  */
  LevelFluxRegister m_levfluxreg;

  ///
  /**
      Has this AMRLevelMG been defined yet?
  */
  bool m_isDefined;

  ///
  /**
      Coarsened version of the grids at this level
  */
  DisjointBoxLayout m_coarsenedGrids;

  ///
  /**
     Residual
  */
  LevelData<FArrayBox> m_resid;

  ///
  /**
     Correction
  */
  LevelData<FArrayBox> m_corr;

  ///
  /**
     The infamous correction to the correction
  */
  LevelData<FArrayBox> m_dcorr;

  ///
  /**
      Storage for L(phi)
  */
  LevelData<FArrayBox> m_lofPhi;

  ///
  /**
     Extra copy of solution
  */
  LevelData<FArrayBox> m_phiSave;

  ///
  /**
      Storage for coarsened residual
  */
  LevelData<FArrayBox> m_resC;

  ///
  /**
      Copier that operates for moving  coarsened residual to
      coarse level grid space.
  */
  Copier m_residualCopier;

  /**
     Copier that's used for exchange on the fine grid data
     before fine flux register data is incremented.
  */
  Copier m_fineExchangeCopier;

  ///
  /**
     Which AMR Level we are on
  */
  int m_level;

  ///
  /**
     The LevelOp which encapsulates which operator we're using
     in this solve.
   */
  LevelOp* m_levelopPtr;

  /**
     Used for averaging residual to coarser level
  */
  CoarseAverage m_averageOp;

private:
  AMRLevelMG& operator = (const AMRLevelMG&);
  AMRLevelMG(const AMRLevelMG&);
};

#include "NamespaceFooter.H"
#endif
