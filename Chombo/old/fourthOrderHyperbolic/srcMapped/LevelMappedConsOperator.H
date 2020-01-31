#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELMAPPEDCONSOPERATOR_H_
#define _LEVELMAPPEDCONSOPERATOR_H_

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FourthOrderFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "GodunovPhysics.H"
#include "PatchMappedConsOperator.H"
#include "TimeInterpolatorRK4.H"
#include "PhysIBC.H"
#include "FourthOrderCoordSys.H"
#include "Tuple.H"
#include "UsingNamespace.H"

/// Level Cons Operator - a class for evaluating the Laplacian.

/**
 */
class LevelMappedConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelMappedConsOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelMappedConsOperator();

  /// Actual constructor.
  /**
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const
     (coding standards).  a_refine is the refinement ratio between
     this level and the next coarser level.  For the coarsest level, an
     empty DisjointBoxLayout is passed in for coarserDisjointBoxLayout.
   */
  void define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              const ProblemDomain&      a_domain,
              const int&                a_refineCoarse,
              const Real&               a_dx,
              const GodunovPhysics* const  a_gdnvPhysics,
              const int&                a_numStates,
              const bool&               a_hasCoarser,
              const bool&               a_hasFiner);

  /// Evaluate the operator (div(flux) ) at a given time.
  /**
     For the coarsest level, UcoarseOld and UcoarseNew are empty
     LevelData<FArrayBox> objects.  Also, either UcoarseOld or
     UcoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  U must be defined on the same
     DisjointBoxLayouts as were used in define().  Coarse flux register
     is flux register with the next coarser level.  Fine flux register is
     the flux register with the next finer level.  To the finer level FR,
     this level is the coarse level.  To the coarser level FR, this level
     is the fine level.
   */
  void evalRHS(LevelData<FArrayBox>&       a_LofU,
               LevelData<FArrayBox>&       a_U,
               LevelFluxRegister&          a_finerFluxRegister,
               LevelFluxRegister&          a_coarserFluxRegister,
               const LevelData<FArrayBox>& a_UcoarseOld,
               const Real&                 a_timeCoarseOld,
               const LevelData<FArrayBox>& a_UcoarseNew,
               const Real&                 a_timeCoarseNew,
               Real                        a_time,
               Real                        a_fluxweight);


  /// update solution -- soln += dt*rhs (required by LevelRK4)
  void updateODE(LevelData<FArrayBox>& a_soln,
                 const LevelData<FArrayBox>& a_rhs,
                 Real a_dt);


  /// reset fluxes contained in this object to zero
  void resetFluxes();


  /// returns reference to fluxes computed and accumulated by this operator.
  LevelData<FluxBox>& getFluxes();


  /// define newSoln to match existingSoln, including ghost cells
  /** (required by LevelRK4)
   */
  void defineSolnData(LevelData<FArrayBox>& a_newSoln,
                      const LevelData<FArrayBox>& a_existingSoln);


  /// define RHS data based on existingSoln (in this case, w/o ghost cells)
  /** (required by LevelRK4)
   */
  void defineRHSData(LevelData<FArrayBox>& a_newRHS,
                     const LevelData<FArrayBox>& a_existingSoln);


  /// copy data from src->dest
  /** (required by LevelRK4)
   */
  void copySolnData(LevelData<FArrayBox>& a_dest,
                    const LevelData<FArrayBox>& a_src);

  /// set spatial order of accuracy
  /*
    Can be 2 or 4 (default)
   */
  void spaceOrder(int a_spaceOrder);


  /// set coordinate system object
  void setCoordSys(FourthOrderCoordSys* a_coordSysPtr);

  /// if true, limit face values
  void limitFaceValues(bool a_limitFaceValues);

  /// sets whether to flatten extrapolations to faces
  void useFlattening(bool a_useFlattening);

  /// sets m_flattening
  void defineFlattening();

  /// sets m_noPPM
  void noPPM(bool a_noPPM);

  /// sets m_doDeconvolution
  void doDeconvolution(bool a_doDeconvolution);

  /// sets m_doFaceDeconvolution
  void doFaceDeconvolution(bool a_doFaceDeconvolution);

  /// sets whether to use artificial viscosity (from divergence)
  void useArtificialViscosity(bool a_useArtificialViscosity);

  /// sets coefficient of artificial viscosity (from divergence)
  void artificialViscosity(Real a_artificialViscosity);

  /// set number of calls that will be made to evalRHS
  void evalCountMax(int a_evalCountMax);

  /// reset m_evalCount to 0
  void resetEvalCount();

  /// add artificial viscosity to a_Unew
  void addArtificialViscosity(LevelData<FArrayBox>&        a_Unew,
                              const LevelData<FArrayBox>&  a_Uold,
                              LevelFluxRegister&   a_finerFluxRegister,
                              LevelFluxRegister&   a_coarserFluxRegister,
                              Real  a_oldTime,
                              Real  a_weight);

  /// update a_finerFluxRegister and a_coarserFluxRegister
  void updateFluxTotalsAndRegisters(// used to update flux registers:  it is shifted and shifted back, but data remain unchanged
                                    FluxBox&             a_Fface,
                                    LevelFluxRegister&   a_finerFluxRegister,
                                    LevelFluxRegister&   a_coarserFluxRegister,
                                    const DataIndex&     a_dataIndex,
                                    Real                 a_weight);

  /// Return time interpolator object within patch-fill object.
  TimeInterpolatorRK4& getTimeInterpolator();

  /// fill in ghost cells of a_U
  void fillGhosts(LevelData<FArrayBox>&   a_U,
                  const Real&             a_time,
                  const Real&             a_timeCoarseOld,
                  const Real&             a_timeCoarseNew);

  void cellUJToCellU(
                     LevelData<FArrayBox>& a_Uavg,
                     const LevelData<FArrayBox>& a_UJavg,
                     const LevelData<FArrayBox>& a_J) const;

protected:

  // put this in PatchMappedConsOperator
  // as well as computeFaceAreas computeFaceGradU computeNJInvGradUOnBox
  // computeNTNJInvGradUOnBox computeNormalDiffusiveFluxVector
  // secondOrderMappedGridDivergence divideOutCellVolume secondOrderLaplacian

  void computeCompFaceFluxes(LevelData<FluxBox>& uTimesV,
                             const LevelData<FluxBox>& u,
                             const LevelData<FluxBox>& v) const;

  void computePhysFaceFluxes(LevelData<FluxBox>&,
                             const LevelData<FluxBox>&);

  void getPhysicalCellVolumes(LevelData<FArrayBox>&) const;

  void computeFaceAverages(LevelData<FluxBox>& face_data,
                           const LevelData<FArrayBox>& cell_data ) const;

  void secondOrderLaplacian( LevelData<FArrayBox>&       a_lap_u,
                             const LevelData<FArrayBox>& a_u,
                             const bool a_div_cell_volume = true );

  void secondOrderMappedGridDivergence( LevelData<FArrayBox>&      a_divF,
                                        const LevelData<FluxBox>&  a_F);

  void computeNormalDiffusiveFluxVector( LevelData<FluxBox>& a_flux,
                                         const LevelData<FArrayBox>& a_u );

  void computeNJInvGradUOnBox( FluxBox&                   a_NJInvGradU,
                               const FluxBox&             a_NJInv,
                               const FluxBox&             a_gradU,
                               const int                  a_comp );

  void computeNTNJInvGradUOnBox( FluxBox&                   a_fluxNormal,
                                 const FluxBox&             a_N,
                                 const FluxBox&             a_NJInvGradU,
                                 const int                  a_comp );

  void computeFaceGradU( LevelData<FluxBox>&         a_grad_u,
                         const LevelData<FArrayBox>& a_u );

  // Box layout for this level
  DisjointBoxLayout m_grids;

  // Patch integrator
  PatchMappedConsOperator m_patchConsOperator;

  // Number of ghost cells need locally for this level
  int m_numGhost;

  // Exchange copier
  Copier m_exchangeCopier;

  // pointer to advection velocities
  //  LevelData<FluxBox>* m_advVelPtr;

  LevelData<FluxBox> m_fluxes;

  FourthOrderCoordSys* m_coordSysPtr;

  // Interpolator for filling in ghost cells from the next coarser level
  // QuadCFInterp m_patcher;
  // PiecewiseLinearFillPatch m_patcher;
  FourthOrderFillPatch m_patcher;

  // Grid spacing
  Real m_dx;

  // Problem domain - index space for this level
  ProblemDomain m_domain;

  // Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  // Number of conserved variables
  int m_numCons;

  // Number of conservative fluxes and other face centered quantities
  int m_numFluxes;

  // Flags for the existence of coarser and finer levels, respectively
  bool m_hasCoarser;
  bool m_hasFiner;

  // flattening coefficients, if m_useFlattening
  BoxLayoutData<FArrayBox> m_flattening;

  // spatial order of accuracy
  int m_spaceOrder;

  // if true, limit face values
  bool m_limitFaceValues;

  // if true, use flattening
  bool m_useFlattening;

  int m_numFields;

  // Avoid PPM?  Default false.
  bool m_noPPM;

  // Do deconvolution?  Default yes.
  bool m_doDeconvolution;

  // Do face deconvolution?  Default yes.
  bool m_doFaceDeconvolution;

  // number of calls that will be made to evalRHS
  int m_evalCountMax;

  // number of calls that have been made to evalRHS so far
  int m_evalCount;

  // use artificial viscosity (from divergence)?
  bool m_useArtificialViscosity;

  // coefficient of artificial viscosity (from divergence)
  Real m_artificialViscosity;

  // [face direction][component]
  Tuple<IntVect, SpaceDim> m_metricTermComponents;

  // define() has been called
  bool m_defined;

private:
  // Disallowed for all the usual reasons
  void operator=(const LevelMappedConsOperator& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelMappedConsOperator(const LevelMappedConsOperator& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
