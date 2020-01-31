#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "LevelAdvectOperator.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilF_F.H"
#include "LevelFluxRegister.H"
#include "QuadCFInterp.H"
#include "timeInterp.H"
#include "AMRIO.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "BoxIterator.H"

#include "mappedLimiter.H"
#include "mappedAdvectionFlux.H"

#include "LoHiSide.H"
#include <sstream>
#include <fstream>

#include <complex>
#include <algorithm>

#define NEWHVFD 1

extern "C"
{
   void FORTRAN_NAME( DNAUPD , dnaupd )
      (int *ido, const char *bmat, int *n, const char *which,
       int *nev, double *tol, double *resid,
       int *ncv, double *V, int *ldv,
       int *iparam, int *ipntr, double *workd,
       double *workl, int *lworkl, int *info);

   void FORTRAN_NAME( DNEUPD, dneupd )
      (int *rvec, char *HowMny, int *select,
       double *dr, double *di, double *Z,
       int *ldz, double *sigmar,
       double *sigmai, double *workev,
       const char *bmat, int *n, const char *which,
       int *nev, double *tol, double *resid,
       int *ncv, double *V, int *ldv,
       int *iparam, int *ipntr,
       double *workd, double *workl,
       int *lworkl, int *info);

   void FORTRAN_NAME( DAXPY, daxpy )
      (int*, double*, double*, int*, double*, int*);

   double FORTRAN_NAME( DNRM2, dnrm2 )(int*, double*, int*);

   double FORTRAN_NAME( DLAPY2, dlapy2 )(double*, double*);

}

// Constructor - set up some defaults
LevelAdvectOperator::LevelAdvectOperator()
{
  m_defined = false;
  m_refineCoarse = 0;
}

// Destructor - free up storage
LevelAdvectOperator::~LevelAdvectOperator()
{
}

// Define the object so that time stepping can begin
void LevelAdvectOperator::define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                                 const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                                 const ProblemDomain&      a_domain,
                                 const int&                a_refineCoarse,
                                 const int&                a_numFields,
                                 const RealVect&           a_dx,
                                 const bool&               a_hasCoarser,
                                 const bool&               a_hasFiner)
{
  // Sanity checks
 CH_assert(a_refineCoarse > 0);
 CH_assert(a_dx > RealVect::Zero);

  // Make a copy of the current grids
  m_grids  = a_thisDisjointBoxLayout;

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;
  m_numFields = a_numFields;

  // set up storage
  m_fluxes.define(m_grids, m_numFields);

  // NTvel nees one layer of ghost cells
  m_NTvel.define(m_grids, 1, IntVect::Unit);

  // Set up the interpolator if there is a coarser level
  if (m_hasCoarser)
    {
      m_patcher.define(a_thisDisjointBoxLayout,
                       &a_coarserDisjointBoxLayout,
                       m_dx[0],
                       a_refineCoarse,
                       m_numFields,
                       a_domain);
    }

  // Everything is defined
  m_defined = true;
}

inline void
LevelAdvectOperator::computeFaceAverages(
   LevelData<FluxBox>& a_face_data,
   const LevelData<FArrayBox>& a_cell_data ) const
{
   if (useFourthOrder())
   {
      fourthOrderCellToFace( a_face_data, a_cell_data );
   }
   else if (useSecondOrder())
   {
      CellToEdge( a_cell_data, a_face_data );
   }
   else
   {
      MayDay::Error("Bad Space Order in LevelAdvectOperator");
   }
}


inline Real ipow( const Real x, const int p )
{
   Real result = 1.0;
   for (int i=0; i<p; i++) result *= x;
   return result;
}


inline void
LevelAdvectOperator::getPhysicalCellVolumes(
   LevelData<FArrayBox>& a_volumes ) const
{
//   if (useSecondOrder())
//   {
//      SecondOrderCoordSys* coordSysPtr
//         = dynamic_cast<SecondOrderCoordSys*>(m_coordSysPtr);
      a_volumes.define( m_coordSysPtr->getCellVolumes() );
//   }
//   else if (useFourthOrder())
//   {
//      FourthOrderCoordSys* coordSysPtr
//         = dynamic_cast<FourthOrderCoordSys*>(m_coordSysPtr);
//      a_volumes.define( coordSysPtr->getCellVolumes() );
//   }
//   else {
//     MayDay::Error("Bad Space Order in LevelAdvectOperator");
//   }
}

// Evaluate the advect operator (-div(vel*phi) ) at time a_time.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// The flux registers are incremented with the normal derivatives of a_phi
// at the boundary, times a fraction of the time step corresponding to
// which stage of an explicit time-stepping scheme (e.g. Runge-Kutta)
// is being invoked.
void LevelAdvectOperator::evalRHS(
                                  LevelData<FArrayBox>&       a_LOfPhi,
                                  LevelData<FArrayBox>&       a_phi,
                                  LevelFluxRegister&          a_finerFluxRegister,
                                  LevelFluxRegister&          a_coarserFluxRegister,
                                  const LevelData<FArrayBox>& a_phiCoarseOld,
                                  const Real&                 a_TCoarseOld,
                                  const LevelData<FArrayBox>& a_phiCoarseNew,
                                  const Real&                 a_TCoarseNew,
                                  Real                        a_time,
                                  Real                        a_weight)
{
  // Make sure everything is defined
  CH_assert(m_defined);

  m_sumWeights += a_weight;


  // Create temporary storage with a layer of "m_numGhost" ghost cells
  IntVect ivGhost = m_numGhost*IntVect::Unit;

  if (m_hasCoarser)
  {
    // Check that current fine-level time "a_time" falls between the old and new coarse times
    Real alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);
    Real dtCoarse = (a_TCoarseNew - a_TCoarseOld);

    // Truncate the fraction to the range [0,1] to remove floating-point
    // subtraction roundoff effects
    Real eps = 0.04 * dtCoarse / m_refineCoarse;

    if (Abs(alpha) < eps)
    {
      alpha = 0.0;
    }

    if (Abs(1.0-alpha) < eps)
    {
      alpha = 1.0;
    }

    // Current time before old coarse time
    if (alpha < 0.0)
    {
      MayDay::Error( "LevelAdvectOperator::evalRHS: alpha < 0.0");
    }

    // Current time after new coarse time
    if (alpha > 1.0)
    {
      MayDay::Error( "LevelAdvectOperator::evalRHS: alpha > 1.0");
    }

    // Interpolate ghost cells from next coarser level using both space
    // and time interpolation.
    {
      // first, interpolate coarse data to the current time
      LevelData<FArrayBox> phiCoarse; phiCoarse.define(a_phiCoarseOld);
      timeInterp(phiCoarse       ,a_time,
                 a_phiCoarseOld  ,a_TCoarseOld,
                 a_phiCoarseNew  ,a_TCoarseNew,
                 a_phi.interval());
      // use current-time coarse data to fill fine ghost cells
      m_patcher.coarseFineInterp(a_phi,phiCoarse);
    }
  }

  // do domain BC's
  m_basicIBCPtr->ghostCellBC(a_phi, m_domain, *m_coordSysPtr, m_dx[0], a_time);

  // Exchange all the ghost cell data between grids at this level
  a_phi.exchange(a_phi.interval());

  // compute computational-space fluxes -
  //   need all three fluxes on each face, so compCoordFlux has
  //   dimension (SpaceDim * nComp)
  LevelData<FluxBox> compCoordFlux( m_grids,
                                    SpaceDim * a_phi.nComp(),
                                    a_phi.ghostVect() );

  LevelData<FluxBox> faceVel( m_grids, SpaceDim, a_phi.ghostVect() );

  if (useFourthOrder())
    {
      computeMappedFourthOrderFlux(compCoordFlux,
                                   faceVel,
                                   a_phi,
                                   m_advVelPtr,
                                   m_FOCS,
                                   m_dx,
                                   limitFaceValues());
    }
  else
    {
      MayDay::Error("2nd-order mapped advection no longer implemented. Ask Dan");
    }

  // set physical boundary conditions on fluxes
  m_basicIBCPtr->fluxBC(compCoordFlux, m_domain, *m_coordSysPtr,
                        m_dx[0], a_time);

  // compute physical-space flux divergence -
  //   We restrict ourselves to nComp=1, since the
  //   mapped grid functionality currently assumes a scalar equation
  CH_assert( a_phi.nComp()==1 );
  LevelData<FluxBox> tempFlux(compCoordFlux.getBoxes(), a_phi.nComp(),
                              compCoordFlux.ghostVect());

  // break up mapped-grid divergence into its parts so that we
  // can accumulate fluxes
  //m_coordSysPtr->mappedGridDivergence( a_LOfPhi, compCoordFlux );
  m_FOCS->computeMetricTermProductAverage(tempFlux, compCoordFlux);

  m_FOCS->simpleDivergence(a_LOfPhi, tempFlux);

  // accumulate fluxes
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisTotalFlux = m_fluxes[dit];
      FluxBox& thisFlux = tempFlux[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          thisTotalFlux[dir].plus(thisFlux[dir], a_weight);
        }
    }

  // compute N^T*vel and place in tempFlux
  m_FOCS->computeMetricTermProductAverage(tempFlux, faceVel);
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisTotalNTvel = m_NTvel[dit];
      FluxBox& thisNTvel = tempFlux[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          thisTotalNTvel[dir].plus(thisNTvel[dir], a_weight);
        }
    }


  // add hyperviscous flux divergence, if needed
  if (useHyperviscosity())
  {
#if NEWHVFD
     //   note that we compute hyperviscous fluxes using physical-space
     //   cell-averaged phi, which is the computational-space cell-averaged phi
     //   to within second-order
    // (7/16/09 -- DFM) -- commented this out because phi_cg is no
    // longer in this scope. Can fix this if we need it.
    MayDay::Error("hyperviscous stuff is broken -- get Dan to fix it");
    //addHyperviscousFluxDivergence( a_LOfPhi, phi_cg );
#else
     //   note that we compute hyperviscous fluxes using cell-averaged phi*J
     addHyperviscousFluxes( a_LOfPhi, a_phi );
#endif
  } // end if hyperviscosity

  // Beginning of loop through patches/grids.
  //DataIterator dit = m_grids.dataIterator();
  Real hDinv = 1.0 / m_dx.product();
  for (dit.begin(); dit.ok(); ++dit)
  {
     // RHS on this box
     FArrayBox& curLOfPhi = a_LOfPhi[dit()];

     // actually RHS is -div; left hand side is <phiJ>
     curLOfPhi.negate();

     // divide by computational cell volume
     curLOfPhi *= hDinv;

     // No support for AMR (yet)
     CH_assert( m_hasFiner==false && m_hasCoarser==false );
#if 0
    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        // Increment coarse flux register between this level and the next
        // finer level - this level is the next coarser level with respect
        // to the next finer level
        if (m_hasFiner && a_finerFluxRegister.isDefined())
          {
            a_finerFluxRegister.incrementCoarse(curFlux[idir],a_weight,
                                                dit(),a_phi.interval(),
                                                a_phi.interval(),idir);
          }

        // Increment fine flux registers between this level and the next
        // coarser level - this level is the next finer level with respect
        // to the next coarser level
        if (m_hasCoarser && a_coarserFluxRegister.isDefined())
          {
            a_coarserFluxRegister.incrementFine(curFlux[idir],a_weight,
                                                dit(),a_phi.interval(),
                                                a_phi.interval(),
                                                idir,Side::Lo);
            a_coarserFluxRegister.incrementFine(curFlux[idir],a_weight,
                                                dit(),a_phi.interval(),
                                                a_phi.interval(),
                                                idir,Side::Hi);
          }
      }
#endif
  }

}


void
LevelAdvectOperator::updateODE(LevelData<FArrayBox>& a_soln,
                               const LevelData<FArrayBox>& a_rhs,
                               Real a_dt)
{
  CH_assert(a_soln.nComp() == a_rhs.nComp());

  DataIterator dit = a_soln.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    FArrayBox& curSoln = a_soln[dit()];
    FArrayBox curRhs(a_rhs[dit()].box(),a_rhs.nComp());
    curRhs.copy(a_rhs[dit()]);
    curRhs *= a_dt;
    curSoln += curRhs;
  }
}

// reset accumulated storage (fluxes, NTvel) to zero
void
LevelAdvectOperator::resetFluxes()
{
  DataIterator dit = m_fluxes.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_fluxes[dit].setVal(0.0);
      m_NTvel[dit].setVal(0.0);
    }

  m_sumWeights = 0.0;
}

/// returns reference to fluxes computed and accumulated by this operator.
LevelData<FluxBox>&
LevelAdvectOperator::getFluxes()
{
  return m_fluxes;
}


// returns reference to face-centered N^T*vel computed and accumulated
// by this operator.
LevelData<FluxBox>&
LevelAdvectOperator::getNTvel()
{
  return m_NTvel;
}
void
LevelAdvectOperator::defineSolnData(LevelData<FArrayBox>& a_newSoln,
                                    const LevelData<FArrayBox>& a_existingSoln)
{
  int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect& ghostVect = a_existingSoln.ghostVect();

  a_newSoln.define(grids, nComp, ghostVect);
}

void
LevelAdvectOperator::defineRHSData(LevelData<FArrayBox>& a_newRHS,
                                   const LevelData<FArrayBox>& a_existingSoln)
{
  // same as defineSolnData, but w/o ghost cells
   int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect ghostVect = IntVect::Zero;

  a_newRHS.define(grids, nComp, ghostVect);
}


/// copy data from src->dest
void
LevelAdvectOperator::copySolnData(LevelData<FArrayBox>& a_dest,
                                  const LevelData<FArrayBox>& a_src)
{
  a_src.copyTo(a_dest);
}

/// set face-centered advection velocities
void
LevelAdvectOperator::setAdvectionVel(LevelData<FArrayBox>* a_advVelPtr)
{
  m_advVelPtr = a_advVelPtr;
}

/// set boundary condition object
void
LevelAdvectOperator::setBC(BasicIBC* a_basicIBCPtr)
{
  m_basicIBCPtr = a_basicIBCPtr;
}

/// set coordinate system object
void
LevelAdvectOperator::setCoordSys(CoordSys<FArrayBox,FluxBox>* a_coordSysPtr)
{
  m_coordSysPtr = a_coordSysPtr;
  m_FOCS = dynamic_cast<FourthOrderCoordSys*>( m_coordSysPtr );
  CH_assert( m_FOCS );
}

void LevelAdvectOperator::avgdown(LevelData<FArrayBox>& a_phi,
                                  LevelData<FArrayBox>& a_phiCoarse)
{
  // Make sure everything is defined

 CH_assert(m_defined);

  // Make sure there is a coarser level.

 CH_assert(m_hasCoarser);

  // Interpolate fine data from next coarser level.
  // We need to use this only if the number of cells over which we are
  // computing the average of l(phi) in FORT_HOAVGDOWN is equal to
  // the size of the refined region. This is always true if nrefine = 2,
  // but need not be otherwise.

  m_patcher.coarseFineInterp(a_phi,a_phiCoarse);

  // Build the stencils for averaging phi and Laplacian(phi).

  Box unitbox(IntVect::Zero,IntVect::Zero);
  Box avStencil = refine(unitbox,m_refineCoarse);
  int halfRefine = m_refineCoarse/2;
  Box lapStencil = refine(unitbox,2) + IntVect::Unit*(halfRefine - 1);
  int sizeLapStencil = lapStencil.volume();
  int sizeAvStencil = avStencil.volume();

  // Exchange all the ghost cell data between grids at this level
  a_phi.exchange(a_phi.interval());

  // Allocate LevelData for corsened version of fine data.

  DisjointBoxLayout coarsenedGrids;
  coarsen(coarsenedGrids,m_grids,m_refineCoarse);

  LevelData<FArrayBox> coarsenedPhi(coarsenedGrids,a_phi.nComp());

  // Beginning of loop through patches/grids.

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    // The current box
    Box curBox = coarsenedGrids.get(dit());

    // The current grid of conserved variables
    FArrayBox& curPhi = a_phi[dit()];
    FArrayBox& curCoarsenedPhi = coarsenedPhi[dit()];

    FORT_HOAVGDOWN(CHF_CONST_FRA(curPhi),
                   CHF_FRA(curCoarsenedPhi),
                   CHF_CONST_INT(m_refineCoarse),
                   CHF_BOX(curBox),
                   CHF_BOX(avStencil),
                   CHF_CONST_INT(sizeAvStencil),
                   CHF_BOX(lapStencil),
                   CHF_CONST_INT(sizeLapStencil)
                   );

  }
  // copy coarsened data onto coarse data.

    coarsenedPhi.copyTo(a_phi.interval(),a_phiCoarse,a_phi.interval());
}


void
LevelAdvectOperator::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

void
LevelAdvectOperator::useHyperviscosity(bool a_useHyperviscosity)
{
   m_useHyperviscosity = a_useHyperviscosity;
}

void
LevelAdvectOperator::hyperviscosity(Real a_hyperviscosity)
{
   CH_assert(a_hyperviscosity>=0);
   m_hyperviscosity = a_hyperviscosity;
}

void
LevelAdvectOperator::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}



inline
RealVect computeFaceAreas( const RealVect& a_dx )
{
   RealVect face_area( RealVect::Unit );
   for (int ndir=0; ndir<SpaceDim; ndir++)
   {
      for (int tdir=0; tdir<SpaceDim; tdir++)
      {
         if (tdir != ndir)
            face_area[ndir] *= a_dx[tdir];
      }
   }
   return face_area;
}

#if NEWHVFD
inline
void zeroLevelData( LevelData<FluxBox>& a_u )
{
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      a_u[dit].setVal( 0.0 );
   }
}

inline
void LevelAdvectOperator::computeFaceGradU( LevelData<FluxBox>& a_grad_u,
                                            const LevelData<FArrayBox>& a_u )
{
   // compute (grad_xi u) at faces
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      const FArrayBox& thisU( a_u[dit] );
      FluxBox& thisGradU( a_grad_u[dit] );
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
         FArrayBox& thisGradUDir( thisGradU[faceDir] );
         Box thisBox( thisGradUDir.box() );
         FORT_SECONDORDERGRADIENT( CHF_FRA(thisGradUDir),
                                   CHF_CONST_FRA(thisU),
                                   CHF_BOX(thisBox),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_INT(faceDir) );
      }
   }
}

inline
void LevelAdvectOperator::computeNJInvGradUOnBox(
   FluxBox&                   a_NJInvGradU,
   const FluxBox&             a_NJInv,
   const FluxBox&             a_gradU,
   const int                  a_comp )
{
   // Compute the vector NJInverse * grad(u)
   const int sdcomp = SpaceDim * a_comp;
   for (int faceDir=0; faceDir<SpaceDim; faceDir++)
   {

      const FArrayBox& thisNJInv( a_NJInv[faceDir] );
      const FArrayBox& thisGradU( a_gradU[faceDir] );
      FArrayBox& thisNJInvGradU( a_NJInvGradU[faceDir] );
      FArrayBox temp( thisNJInvGradU.box(), 1 );

      for (int d=0; d<SpaceDim; d++)
      {
         for (int s=0; s<SpaceDim; s++)
         {
            int NJInv_comp = m_FOCS->getMetricTermComponent( d, s );
            temp.copy( thisNJInv, NJInv_comp, 0, 1 );
            temp.mult( thisGradU, sdcomp + s, 0, 1 );
            thisNJInvGradU.plus( temp, 0, d, 1 );
         }
      }
   }
}

inline
void LevelAdvectOperator::computeNTNJInvGradUOnBox(
   FluxBox&                   a_fluxNormal,
   const FluxBox&             a_N,
   const FluxBox&             a_NJInvGradU,
   const int                  a_comp )
{
   // Compute the vector NT * NJInverse * GradU
   for (int faceDir=0; faceDir<SpaceDim; faceDir++)
   {

      const FArrayBox& thisN( a_N[faceDir] );
      const FArrayBox& thisNJInvGradU( a_NJInvGradU[faceDir] );
      FArrayBox& thisFluxNormal( a_fluxNormal[faceDir] );
      FArrayBox temp( thisFluxNormal.box(), 1 );

      for (int dir=0; dir<SpaceDim; dir++)
      {
         // Get NT component (indices are reversed to get transpose)
         int NTcomp = m_FOCS->getMetricTermComponent( dir, faceDir );
         temp.copy( thisN, NTcomp, 0, 1 );
         temp.mult( thisNJInvGradU, dir, 0, 1 );
         thisFluxNormal.plus( temp, 0, a_comp, 1 );
      }
   }
}

void
LevelAdvectOperator::computeNormalDiffusiveFluxVector(
   LevelData<FluxBox>& a_flux,
   const LevelData<FArrayBox>& a_u )
{
   // initialize normal fluxes to zero
   zeroLevelData( a_flux );

   // get storage parameters from state vector
   const DisjointBoxLayout& grids = a_u.disjointBoxLayout();
   const int ncomp = a_u.nComp();
   IntVect ghostVect( a_flux.ghostVect() );

   // compute grad(u)
   LevelData<FluxBox> grad_u( grids, SpaceDim * ncomp, ghostVect );
   computeFaceGradU( grad_u, a_u );

   // get metric data
   const LevelData<FluxBox>& NJInv( m_FOCS->getNJinverse() );
   const LevelData<FluxBox>& N( m_FOCS->getFaceMetricTerms() );

   // create intermediate product LevelData
   LevelData<FluxBox> NJInvGradU( grids, SpaceDim, ghostVect );

   // Loop over each component
   DataIterator dit = a_flux.dataIterator();
   for (int comp=0; comp<ncomp; comp++)
   {

      // Initialize for accumulation
      zeroLevelData( NJInvGradU );

      // compute average diffusive flux
      for (dit.begin(); dit.ok(); ++dit)
      {

         computeNJInvGradUOnBox( NJInvGradU[dit],
                                 NJInv[dit],
                                 grad_u[dit],
                                 comp );

         computeNTNJInvGradUOnBox( a_flux[dit],
                                   N[dit],
                                   NJInvGradU[dit],
                                   comp );
      }
   }

   /*
     Now remove the extra area factor resulting from the
     fact that the N and NJinverse returned by the
     coordinate system are face integrals rather than
     averages.
   */
   RealVect face_area( computeFaceAreas( m_dx ) );
   for (int dir=0; dir<SpaceDim; dir++)
   {
      double face_area_reciprocal = 1.0 / face_area[dir];
      for (dit.begin(); dit.ok(); ++dit)
      {
         a_flux[dit][dir].mult( face_area_reciprocal );
      }
   }
}

void
LevelAdvectOperator::secondOrderMappedGridDivergence(
   LevelData<FArrayBox>&      a_divF,
   const LevelData<FluxBox>&  a_F )
{
   DataIterator dit = a_divF.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      const FluxBox& thisFlux = a_F[dit];
      FArrayBox& thisDiv = a_divF[dit];

      // first, set divF to 0
      thisDiv.setVal(0.0);

      // since we're computing volume*div, and the fluxes
      // are multiplied by the appropriate areas, don't need
      // a dx here.
      Real fakeDx = 1.0;

      // now loop over directions and increment with directional
      // derivative
      for (int dir=0; dir<SpaceDim; dir++)
      {
         const Box& thisBox = thisDiv.box();
         // use fortran from divergence operator here
         FORT_DIVERGENCE( CHF_CONST_FRA(thisFlux[dir]),
                          CHF_FRA(thisDiv),
                          CHF_BOX(thisBox),
                          CHF_CONST_REAL(fakeDx),
                          CHF_INT(dir) );
      }
   }
}

inline
void divideOutCellVolume( LevelData<FArrayBox> &     a_lap_u,
                          const FourthOrderCoordSys* a_FOCS )
{
   const LevelData<FArrayBox>& cell_volume = a_FOCS->getCellVolumes();
   DataIterator dit = a_lap_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisDiv = a_lap_u[dit];
      const FArrayBox& thisVolume = cell_volume[dit];
      for (int comp=0; comp<thisDiv.nComp(); comp++)
      {
         thisDiv.divide( thisVolume, thisDiv.box(), 0, comp, 1 );
      }
   }
}

void
LevelAdvectOperator::secondOrderLaplacian(
   LevelData<FArrayBox> &       a_lap_u,
   const LevelData<FArrayBox> & a_u,
   const bool                   a_div_cell_volume )
{
   // create a temporary diffusive flux vector
   const DisjointBoxLayout& grids = a_lap_u.disjointBoxLayout();
   const int fluxComp = a_lap_u.nComp();
   const IntVect& ghostVect = a_lap_u.ghostVect();
   LevelData<FluxBox> normalFlux( grids, fluxComp, ghostVect );

   // compute face-averaged normal diffusive fluxes
   computeNormalDiffusiveFluxVector( normalFlux, a_u );
   normalFlux.exchange();

   // compute diffusive flux divergence times physical cell volume
   secondOrderMappedGridDivergence( a_lap_u, normalFlux );

   // Divide by the physcial cell volume to get the average divergence
   if (a_div_cell_volume)
   {
      divideOutCellVolume( a_lap_u, m_FOCS );
   }

   // exchange internal ghost cells
   a_lap_u.exchange();
}

inline
Real effectiveMinimumPhysicalCellLength(
   const DisjointBoxLayout& boxes,
   const FourthOrderCoordSys *a_FOCS )
{
#if 0
   Real vol_min = HUGE_VAL;
   const LevelData<FArrayBox>& cell_volume = a_FOCS->getCellVolumes();
   DataIterator dit = cell_volume.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      vol_min = std::min( vol_min, cell_volume[dit].min() );
   }
   return pow( vol_min, 1.0 / SpaceDim );
#else
   Real h = (a_FOCS->dx())[0];
   Real h_ratio = HUGE_VAL;
   DataIterator dit = boxes.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      BoxIterator bit(boxes[dit]);
      for (bit.begin(); bit.ok(); ++bit)
      {
         IntVect iv = bit();
         RealVect xi( D_DECL( (iv[0]+0.5)*h,
                              (iv[1]+0.5)*h,
                              (iv[2]+0.5)*h) );
         RealVect x( a_FOCS->realCoord( xi ) );
         for (int dir=0; dir<SpaceDim; dir++)
         {
            Real sum = 0.0;
            for (int sdir=0; sdir<SpaceDim; sdir++)
            {
               sum += a_FOCS->dXdXi( x, sdir, dir );
            }
            h_ratio = std::min( h_ratio, fabs(sum) );
         }
      }
   }
   return h_ratio * h;
#endif
}

// this is the one in physical space
void
LevelAdvectOperator::addHyperviscousFluxDivergence(
   LevelData<FArrayBox>&   a_divF,
   const LevelData<FArrayBox>& a_u )
{
   // Add stabilizing hyperviscous fluxes to the SpaceDim-by-nComp
   // face-averaged fluxes in computational space, where a_u is the
   // nComp-dim state vector
   const DisjointBoxLayout& grids( a_u.getBoxes() );
   IntVect nghost( a_u.ghostVect() );
   const int ncomp = a_u.nComp();

   // compute laplacian u
   nghost -= IntVect::Unit;
   LevelData<FArrayBox> lap_u( grids, ncomp, nghost );
   secondOrderLaplacian( lap_u, a_u );
//   secondOrderLaplacian( lap_u, a_u, false );

   // compute laplacian laplacian u
   nghost -= IntVect::Unit;
   LevelData<FArrayBox> lap_lap_u( grids, ncomp, nghost );
   secondOrderLaplacian( lap_lap_u, lap_u );
//   secondOrderLaplacian( lap_lap_u, lap_u, false );

   // compute laplacian laplacian laplacian u
   nghost -= IntVect::Unit;
   LevelData<FArrayBox> hvFluxDiv( grids, ncomp, nghost );
   secondOrderLaplacian( hvFluxDiv, lap_lap_u, false );

   // add to total flux divergence
//   Real h = effectiveMinimumPhysicalCellLength( grids, m_FOCS );
//   cout << "h/m_dx = " << h / m_dx[0] << endl;
   Real h = m_dx[ m_dx.maxDir(false) ];
   Real coeff = m_hyperviscosity * ipow( h, 5 );
//   Real coeff = -m_hyperviscosity * ipow( h, 3 );
//   Real coeff = m_hyperviscosity * h;
   DataIterator dit = a_divF.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisDiv( a_divF[dit] );
      FArrayBox& thisHVFluxDiv( hvFluxDiv[dit] );
//      FArrayBox& thisHVFluxDiv( lap_lap_u[dit] );
//      FArrayBox& thisHVFluxDiv( lap_u[dit] );
      thisHVFluxDiv *= coeff;
      thisDiv.minus( thisHVFluxDiv, thisDiv.box(), 0, 0, thisDiv.nComp() );
   }
}
#else // this is the one in computational space
void
LevelAdvectOperator::addHyperviscousFluxes( LevelData<FArrayBox>& a_divF,
                                            const LevelData<FArrayBox>& a_u )
{
   // Add stabilizing hyperviscous fluxes to the SpaceDim-by-nComp
   // face-averaged fluxes in computational space, where a_u is the
   // nComp-dim state vector

   // compute hyperviscous fluxes
   const RealVect face_area( computeFaceAreas( m_dx ) );
   LevelData<FluxBox> tempFlux( a_u.getBoxes(), a_u.nComp(), IntVect::Zero );
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FluxBox& thisFlux = tempFlux[dit];
      const FArrayBox& thisU = a_u[dit];

      // loop over directions
      for (int dir=0; dir<SpaceDim; dir++)
      {
         Real fakeDx = 1.0;
         FArrayBox& thisFluxDir = thisFlux[dir];
         thisFluxDir.setVal( 0.0 );
         FORT_ADDHYPERVISCOUSFLUX(CHF_FRA(thisFluxDir),
                                  CHF_CONST_FRA(thisU),
                                  CHF_BOX(thisFluxDir.box()),
                                  CHF_REAL(m_hyperviscosity),
                                  CHF_REAL(fakeDx),
                                  CHF_CONST_INT(dir));
         thisFluxDir *= face_area[dir];
      } // end loop over faceDir
   } // end loop over boxes

   // now compute divergence in the usual way
   dit = a_divF.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisDiv = a_divF[dit];
      const Box& thisBox = thisDiv.box();
      FArrayBox tempDiv( thisBox, thisDiv.nComp() );
      tempDiv.setVal(0.0);

      // since we're computing volume*div, and the fluxes
      // are multiplied by the appropriate areas, don't need
      // a dx here.
      Real fakeDx = 1.0;

      // now loop over directions and increment with directional
      // derivative
      FluxBox& thisFlux = tempFlux[dit];
      for (int dir=0; dir<SpaceDim; dir++)
      {
         // use fortran from divergence operator here
         FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                         CHF_FRA(tempDiv),
                         CHF_BOX(thisBox),
                         CHF_CONST_REAL(fakeDx),
                         CHF_INT(dir));

      }
      thisDiv += tempDiv;
   }
}
#endif

inline
void serializeLevelData( const LevelData<FArrayBox>& a_lOfPhi,
                         Real input_vector[] )
{
   const int ncomp = a_lOfPhi.nComp();
   const DisjointBoxLayout& dbl = a_lOfPhi.disjointBoxLayout();
   DataIterator dit = dbl.dataIterator();
   long count = 0;
   for (dit.begin(); dit.ok(); ++dit)
   {
      const FArrayBox& thisLOfPhi = a_lOfPhi[dit];
      BoxIterator bit( dbl[dit] );
      for (bit.begin(); bit.ok(); ++bit)
      {
         for (int n=0; n<ncomp; n++)
         {
            input_vector[ count++ ] = thisLOfPhi( bit(), n);
         }
      }
   }
}

inline
void deserializeLevelData( const Real input_vector[],
                           LevelData<FArrayBox>& a_Phi )

{
   const int ncomp = a_Phi.nComp();
   LevelData<FArrayBox> phiTmp( a_Phi.getBoxes(), ncomp );
   const DisjointBoxLayout& dbl = phiTmp.disjointBoxLayout();
   DataIterator dit = dbl.dataIterator();
   long count = 0;
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisPhi = phiTmp[dit];
      BoxIterator bit( dbl[dit] );
      for (bit.begin(); bit.ok(); ++bit)
      {
         for (int n=0; n<ncomp; n++)
         {
            thisPhi( bit(), n) = input_vector[ count++ ];
         }
      }
   }
   phiTmp.copyTo( a_Phi );
   a_Phi.exchange();
}

int sizeOfSerialOperandVector( const LevelData<FArrayBox>& a_Phi )
{
   const int ncomp = a_Phi.nComp();
   const DisjointBoxLayout& dbl = a_Phi.getBoxes();
   int size = 0;
   DataIterator dit = dbl.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      const Box& box = dbl[dit];
      size += box.numPts() * ncomp;
   }
   return size;
}

inline
void LevelAdvectOperator::mtrxVectMult( double ax[],
                                        const double v[],
                                        LevelData<FArrayBox>& a_Phi,
                                        LevelData<FArrayBox>& a_LOfPhi,
                                        const Real& a_time,
                                        const int& a_ido )
{
   // copy v to phi
   deserializeLevelData( v, a_Phi );
   if ( a_ido==1 )
   {
      // Performing matrix-vector multiplication.
      // In regular mode, w = Av must be performed whenever
      // GetIdo is equal to 1 or -1. GetVector supplies a pointer
      // to the input vector, v, and PutVector a pointer to the
      // output vector, w.
      LevelFluxRegister dummyFR;
      LevelData<FArrayBox> dummyLD( a_Phi.getBoxes(), a_Phi.nComp() );
      Real dummyR = 0;
      evalRHS( a_LOfPhi, a_Phi, dummyFR, dummyFR,
               dummyLD, dummyR, dummyLD, dummyR, a_time, dummyR );
      // copy LOFPhi to ax
      serializeLevelData( a_LOfPhi, ax );
   }
}

inline
void cdaxpy( int n, Real dval, Real v[], Real ax[] )
{
   int one = 1;
   FORTRAN_NAME(DAXPY,daxpy)( &n, &dval, v, &one, ax, &one );
}

inline
Real cdnrm2( int n, Real ax[] )
{
   int one = 1;
   return FORTRAN_NAME(DNRM2,dnrm2)( &n, ax, &one );
}

inline
Real cdlapy2( Real d1, Real d2 )
{
   return FORTRAN_NAME(DLAPY2,dlapy2)( &d1, &d2 );
}

class Ritz
{
   public:
      std::complex<double> lambda;
      double residual;

      Ritz()
        :
        lambda(0.0,0.0),
        residual(0.0)
      {
      }

      Ritz( const Ritz& rhs )
      {
         this->lambda = rhs.lambda;
         this->residual = rhs.residual;
      }

      Ritz& operator=( const Ritz& rhs )
      {
         if (this != &rhs)
         {
            this->lambda = rhs.lambda;
            this->residual = rhs.residual;
         }
         return *this;
      }

      void print( std::ostream& out )
      {
         if ( lambda.real() >= 0 )
            out << " ";
         out << lambda.real();
         if ( lambda.imag() >= 0 )
            out << " + ";
         else
            out << " - ";
         out << abs(lambda.imag())
             << "     " << residual << endl;
      }
};


bool real_less( Ritz a, Ritz b )
{
   return ( a.lambda.real() < b.lambda.real() );
}

bool real_greater( Ritz a, Ritz b )
{
   return ( a.lambda.real() > b.lambda.real() );
}

bool imag_less( Ritz a, Ritz b )
{
   bool result = ( a.lambda.imag() < b.lambda.imag() );
   return result;
}

bool imag_greater( Ritz a, Ritz b )
{
   return ( a.lambda.imag() > b.lambda.imag() );
}

bool mag_less( Ritz a, Ritz b )
{
   return ( std::abs( a.lambda ) < std::abs( b.lambda ) );
}

bool mag_greater( Ritz a, Ritz b )
{
   return ( std::abs( a.lambda ) > std::abs( b.lambda ) );
}

bool res_less( Ritz a, Ritz b )
{
   return ( a.residual < b.residual );
}

void mysort( const std::string& which, std::vector<Ritz>& ritz )
{
   if ( which== "SR" )
      std::sort( ritz.begin(), ritz.end(), real_less );
   else if ( which== "LR" )
      std::sort( ritz.begin(), ritz.end(), real_greater );
   else if ( which== "SI" )
      std::sort( ritz.begin(), ritz.end(), imag_less );
   else if ( which== "LI" )
      std::sort( ritz.begin(), ritz.end(), imag_greater );
   else if ( which== "SM" )
      std::sort( ritz.begin(), ritz.end(), mag_less );
   else if ( which== "LM" )
      std::sort( ritz.begin(), ritz.end(), mag_greater );
   else
      std::sort( ritz.begin(), ritz.end(), res_less );
}


void
LevelAdvectOperator::estimateEigenvalues( const LevelData<FArrayBox>& a_Phi,
                                          const Real&                 a_time,
                                          const Real&                 a_dt )
{

#if 0
//   const int MAXN   = 256*256;
//   const int MAXNEV = 256;
//   const int MAXNCV = 1000;
   const int MAXN   = 400*400;
   const int MAXNEV = 400;
   const int MAXNCV = 1600;
   int LDV    = MAXN;

   // compute the number of entries in work array
   int n   = sizeOfSerialOperandVector( a_Phi );
   int nev = n-2;
   int ncv = n;
   if ( n > MAXN )     MayDay::Abort( "ERROR: N is greater than MAXN" );
   if ( nev > MAXNEV ) MayDay::Abort( "ERROR: NEV is greater than MAXNEV" );
   if ( ncv > MAXNCV ) MayDay::Abort( "ERROR: NCV is greater than MAXNCV" );

   // ARPACK parameters
   std::string bmat( "I" );
   std::vector<std::string> which(2);
   which[0] = "SI";
   which[1] = "LI";
   int lworkl     = ncv * ( 3 * ncv + 6 );
   int info       = 0;
   int maxitr     = 3000;
   int ishfts     = 1;
   int mode       = 1;
   double tol     = 0.0;
   std::vector<int> ipntr(14,0);
   std::vector<int> iparam(11,0);
   iparam[0] = ishfts;
   iparam[2] = maxitr;
   iparam[6] = mode;

   // ARPACK work arrays
   std::vector<double> resid(MAXN);
   std::vector<double> workd(3*MAXN);
   std::vector<double> workl(MAXNCV*(3*MAXNCV+6));
   std::vector<double> v(MAXNCV*LDV);

   // make temporary leveldata
   LevelData<FArrayBox> lOfPhi( a_Phi.getBoxes(), a_Phi.nComp() );
   LevelData<FArrayBox> phiTmp( a_Phi.getBoxes(), a_Phi.nComp(), a_Phi.ghostVect() );
   a_Phi.copyTo( phiTmp );

   // copy a_Phi to workd(ipntr(1))
   {
//      LevelData<FArrayBox> phiTmpNG( a_Phi.getBoxes(), a_Phi.nComp() );
//      a_Phi.copyTo( phiTmpNG );
//      serializeLevelData( phiTmpNG, &(workd[ ipntr[0]-1 ]) );
   }

// two iterations: first find SI, then find LI, writing out each time
   std::ofstream ritzFile( "ritzvalues.dat", ios::out );
   if (!ritzFile)
   {
      MayDay::Abort( "File could not be opened" );
   }

   for ( int factorization=0; factorization<2; factorization++ )
   {

      // Find Arnoldi basis
      int ido = 0;
      for (int i=0; i<14; i++) ipntr[i] = 0;
      int count = 0;
      while (abs(ido)<2)
      {
         FORTRAN_NAME( DNAUPD , dnaupd )( CHF_INT(ido),
                                          bmat.c_str(),
                                          CHF_INT(n),
                                          which[factorization].c_str(),
                                          CHF_INT(nev),
                                          CHF_REAL(tol),
                                          &(resid[0]),
                                          CHF_INT(ncv),
                                          &(v[0]),
                                          CHF_INT(LDV),
                                          &(iparam[0]),
                                          &(ipntr[0]),
                                          &(workd[0]),
                                          &(workl[0]),
                                          CHF_INT(lworkl),
                                          CHF_INT(info) );

         double* x  = &(workd[ipntr[0]-1]);
         double* ax = &(workd[ipntr[1]-1]);
         mtrxVectMult( ax, x, phiTmp, lOfPhi, a_time, abs(ido) );
         count++;
//      cout << count << "\t" << resid[0] << "\t" << resid[1] << endl;
      }

      // Either we have convergence or there is an error.
      if ( info < 0 )
      {
         // Error message, check the documentation in DNAUPD.
         ostringstream msg;
         msg << endl << " Error with _naupd, info = " << info
             << " Check the documentation of _naupd" << endl;
         MayDay::Abort( msg.str().c_str() );
      }
      else
      {
        // No fatal errors occurred.

         /*
          * Post-Process using DNEUPD:
          *   Computed eigenvalues may be extracted.
          *   Eigenvectors may also be computed now if desired. (rvec = true)
          */
         int rvec = 1;
         int ierr = 0;
         char howmany[2] = "A";
         double sigmar = 0.0;
         double sigmai = 0.0;
         std::vector<int> select(MAXNCV,0);
         std::vector<double> d(MAXNCV*3);
         std::vector<double> workev(3*MAXNCV);

         FORTRAN_NAME( DNEUPD , dneupd ) ( &rvec,
                                           howmany,
                                           &(select[0]),
                                           &(d[0]),
                                           &(d[MAXNCV]),
                                           &(v[0]),
                                           CHF_INT(LDV),
                                           CHF_REAL(sigmar),
                                           CHF_REAL(sigmai),
                                           &(workev[0]),
                                           bmat.c_str(),
                                           CHF_INT(n),
                                           which[factorization].c_str(),
                                           CHF_INT(nev),
                                           CHF_REAL(tol),
                                           &(resid[0]),
                                           CHF_INT(ncv),
                                           &(v[0]),
                                           CHF_INT(LDV),
                                           &(iparam[0]),
                                           &(ipntr[0]),
                                           &(workd[0]),
                                           &(workl[0]),
                                           CHF_INT(lworkl),
                                           CHF_INT(ierr) );

         /*
          * The real part of the eigenvalue is returned
          * in the first column of the two dimensional
          * array D, and the imaginary part is returned
          * in the second column of D.  The corresponding
          * eigenvectors are returned in the first NEV
          * columns of the two dimensional array V if
          * requested.  Otherwise, an orthogonal basis
          * for the invariant subspace corresponding to
          * the eigenvalues in D is returned in V.
          */
         if ( ierr != 0 )
         {
            ostringstream msg;
            msg << endl << " Error with _neupd, info = " << ierr
                << " Check the documentation of _neupd" << endl;
            MayDay::Abort( msg.str().c_str() );
         }
         else
         {
           // No fatal errors occurred.

            bool first = true;
            int nconv = iparam[4];
            std::vector<double> axv(MAXN,0.0);
            Real* ax = &(axv[0]);
            for (int j=0; j<nconv; j++)
            {

               /*
                * Compute the residual norm
                *
                * |  A*x - lambda*x ||
                *
                * for the NCONV accurately computed eigenvalues and
                * eigenvectors.  (iparam(4) indicates how many are
                * accurate to the requested tolerance)
                */
               if (d[j+MAXNCV] == 0.0)
               {
                 // Ritz value is real

                  mtrxVectMult( ax, &(v[j*LDV]), phiTmp, lOfPhi, a_time, 1 );
                  cdaxpy( n, -d[j], &(v[j*LDV]), ax );
                  d[j+2*MAXNCV] = cdnrm2( n, ax );
                  d[j+2*MAXNCV] /= fabs( d[j] );

               }
               else if (first)
               {

                  mtrxVectMult( ax, &(v[j*LDV]), phiTmp, lOfPhi, a_time, 1 );
                  cdaxpy( n,       -d[j], &(v[j*LDV]), ax );
                  cdaxpy( n, d[j+MAXNCV],     &(v[(j+1)*LDV]), ax );
                  d[j+2*MAXNCV] = cdnrm2( n, ax );

                  mtrxVectMult( ax, &(v[(j+1)*LDV]), phiTmp, lOfPhi, a_time, 1 );
                  cdaxpy( n, -d[j+MAXNCV], &(v[j*LDV]), ax );
                  cdaxpy( n,        -d[j],     &(v[(j+1)*LDV]), ax );

                  double norm = cdnrm2( n, ax );
                  d[j+2*MAXNCV] = cdlapy2( d[j+2*MAXNCV], norm );
                  d[j+2*MAXNCV] /= cdlapy2( d[j], d[j+MAXNCV] );
                  d[j+1+2*MAXNCV] = d[j+2*MAXNCV];
                  first = false;

               }
               else
               {
                  first = true;
               }
            }

            std::vector<Ritz> ritz( nconv );
            for (int j=0; j<nconv; j++)
            {
               ritz[j].lambda = std::complex<double>( d[j], d[j+MAXNCV] );
               ritz[j].residual = d[j+2*MAXNCV];
            }

            // sort the ritz values
            mysort( which[factorization], ritz );

            // display
            cout.precision(4);
            cout.setf( ios::scientific );
//       cout << "Eigenvalues:" << endl;
            cout << "Sorted:" << endl;
         for (int j=0; j<nconv; j++)
         {
            cout << j+1 << "\t";
            ritz[j].print( cout );
            ritzFile << (ritz[j].lambda.real()*a_dt) << " "
                     << (ritz[j].lambda.imag()*a_dt) << endl;
         }
         cout << endl;

         }

         // Print additional convergence information.
         if ( info == 1)
            cout << endl << " Maximum number of iterations reached." << endl;
         else if ( info == 3)
            cout << endl << " No shifts could be applied during implicit "
                 << "Arnoldi update, try increasing NCV." << endl;

         cout << endl;
         cout << " Size of the matrix is " << n << endl;
         cout << " The number of Ritz values requested is " << nev << endl;
         cout << " The number of Arnoldi vectors generated (NCV) is "
              << ncv << endl;
         cout << " What portion of the spectrum: " <<  which[factorization] << endl;
         cout << " The number of converged Ritz values is " << iparam[4] << endl;
         cout << " The number of Implicit Arnoldi update iterations taken is "
              << iparam[2] << endl;
         cout << " The number of OP*x is " << iparam[8] << endl;
         cout << " The convergence criterion is " << tol << endl;
         cout << endl;
      }
   }

   ritzFile.close();
#endif

}


