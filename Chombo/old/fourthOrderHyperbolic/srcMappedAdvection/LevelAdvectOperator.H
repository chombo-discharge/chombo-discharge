#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELADVECTOPERATOR_H_
#define _LEVELADVECTOPERATOR_H_

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "QuadCFInterp.H"
#include "ExtrapFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "AdvectOpF_F.H"
#include "BasicIBC.H"
//#include "SecondOrderCoordSys.H"
#include "FourthOrderCoordSys.H"
#include "UsingNamespace.H"

/// Level Advect Operator - a class for evaluating the Laplacian.

/**
 */
class LevelAdvectOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelAdvectOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelAdvectOperator();

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
              const int&                a_numStates,
              const RealVect&           a_dx,
              const bool&               a_hasCoarser,
              const bool&               a_hasFiner);

  /// Evaluate the operator (div(flux) ) at a given time.
  /**
     For the coarsest level, phiCoarseOld and phiCoarseNew are empty
     LevelData<FArrayBox> objects.  Also, either phiCoarseOld or
     phiCoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  phi must be defined on the same
     DisjointBoxLayouts as were used in define().  Coarse flux register
     is flux register with the next coarser level.  Fine flux register is
     the flux register with the next finer level.  To the finer level FR,
     this level is the coarse level.  To the coarser level FR, this level
     is the fine level.
   */
  void evalRHS(LevelData<FArrayBox>&       a_LOfPhi,
               LevelData<FArrayBox>&       a_Phi,
               LevelFluxRegister&          a_finerFluxRegister,
               LevelFluxRegister&          a_coarserFluxRegister,
               const LevelData<FArrayBox>& a_phiCoarseOld,
               const Real&                 a_TCoarseOld,
               const LevelData<FArrayBox>& a_phiCoarseNew,
               const Real&                 a_TCoarseNew,
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

  /// returns reference to N^T*vel computed and accumulated by this operator.
  LevelData<FluxBox>& getNTvel();

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

  /// set face-centered advection velocities
  void setAdvectionVel(LevelData<FArrayBox>* a_advVelPtr);

  /// accumulate N^T*vel here
  LevelData<FluxBox> m_NTvel;

  /// accumulate fluxes
  LevelData<FluxBox> m_fluxes;

  Real m_sumWeights;

  /// set boundary condition object
  void setBC(BasicIBC* a_basicIBCPtr);

  /// set coordinate system object
  void setCoordSys(CoordSys<FArrayBox,FluxBox>* a_coordSysPtr);

  /// Average down from the fine grid to the coarse grid.
  /**
     This operator averages the solution from the fine grid onto the
     grid at the next coarser AMR level. Our averaging operator given as

     phiCoarse = <phi> - h^2 <L(phi)>_2 / 24 .

     Here <.> denotes the arithmetic average over the fine cells covered
     by each coarse cell, L(phi) is the (2*SpaceDim + 1)-point
     discretization of the Laplacian, with ghost cells set by QuadCFInterp,
     and <.>_2 is the average over the fine subblock of size 2 located in the
     center of the coarse cell.

  */

  void avgdown(LevelData<FArrayBox>&       a_Phi,
               LevelData<FArrayBox>&       a_phiCoarse);


  /// stabilizes the scheme using hyperviscosity
  void addHyperviscousFluxDivergence( LevelData<FArrayBox>& a_divF,
                                      const LevelData<FArrayBox>& a_u );
  void addHyperviscousFluxes( LevelData<FArrayBox>& a_divF,
                              const LevelData<FArrayBox>& a_u );


  /// set spatial order of accuracy
  /*
    Can be 2 or 4 (default)
   */
  void spaceOrder(int a_spaceOrder);

  /// if true, limit face values
  void limitFaceValues(bool a_limitFaceValues);

  bool limitFaceValues() const
  {
    return m_limitFaceValues;
  }

  void useHyperviscosity(bool a_useHyperviscosity);

  bool useHyperviscosity() const
  {
    return m_useHyperviscosity;
  }

  void hyperviscosity(Real a_hyperviscosity);

  Real hyperviscosity() const
  {
    return m_hyperviscosity;
  }

  bool useFourthOrder() const
  {
    return m_spaceOrder==4;
  }

  bool useSecondOrder() const
  {
    return m_spaceOrder==2;
  }

  void estimateEigenvalues( const LevelData<FArrayBox>& a_Phi,
                            const Real& a_time,
                            const Real& a_dt );

protected:
  // define() has been called
  bool m_defined;

  // Box layout for this level
  DisjointBoxLayout m_grids;

  // pointer to advection velocities
  LevelData<FArrayBox>* m_advVelPtr;

  // boundary condition object
  BasicIBC* m_basicIBCPtr;

  // boundary condition object
  CoordSys<FArrayBox,FluxBox>* m_coordSysPtr;
  FourthOrderCoordSys* m_FOCS;

  // Interpolator for filling in ghost cells from the next coarser level
  QuadCFInterp m_patcher;

  // Grid spacing
  RealVect m_dx;

  // Problem domain - index space for this level
  ProblemDomain m_domain;

  // spatial order of accuracy
  int m_spaceOrder;

  // if true, limit face values
  bool m_limitFaceValues;

  // if true, add hyperviscous fluxes to each flux for stabilization
  bool m_useHyperviscosity;
  Real m_hyperviscosity;

  // Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  // Number of conservative fluxes and other face centered quantities
  int m_numFluxes;

  int m_numGhost;
  int m_numFields;

  // Flags for the existence of coarser and finer levels, respectively
  bool m_hasCoarser;
  bool m_hasFiner;

private:
  // Disallowed for all the usual reasons
  void operator=(const LevelAdvectOperator& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelAdvectOperator(const LevelAdvectOperator& a_input)
  {
    MayDay::Error("invalid operator");
  }

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

  /// estimate eigenvalues of RHS matrix
  void mtrxVectMult( Real ax[],
                     const Real v[],
                     LevelData<FArrayBox>& a_Phi,
                     LevelData<FArrayBox>& a_LOfPhi,
                     const Real& a_time,
                     const int& a_ido);

};

#endif
