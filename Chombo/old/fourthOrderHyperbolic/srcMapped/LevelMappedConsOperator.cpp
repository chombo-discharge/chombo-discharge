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

#include "LevelMappedConsOperator.H"
#include "FourthOrderUtil.H" // functions with names beginning "fourthOrder"
// #include "AMRIO.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "AdvectOpF_F.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
LevelMappedConsOperator::LevelMappedConsOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
LevelMappedConsOperator::~LevelMappedConsOperator()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void LevelMappedConsOperator::define(const DisjointBoxLayout&  a_thisLayout,
                               const DisjointBoxLayout&  a_coarserLayout,
                               const ProblemDomain&      a_domain,
                               const int&                a_refineCoarse,
                               const Real&               a_dx,
                               const GodunovPhysics* const a_gdnvPhysics,
                               const int&                a_numFields,
                               const bool&               a_hasCoarser,
                               const bool&               a_hasFiner)
{
  // Sanity checks
  CH_assert(a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisLayout;

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;

  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;
  m_numFields = a_numFields;

  // petermc added 6 July 2008; changed from 4 to 5 on 25 Aug.
  // changed from 4 to 5 on 20 Aug 2009
  m_numGhost = 5;

  // allocate storage for fluxes
  // LevelData<FluxBox> m_fluxes;
  m_fluxes.define(m_grids, a_numFields);

  // Set up the interpolator if there is a coarser level
  // WAS QuadCFInterp m_patcher;
  if (m_hasCoarser)
    {
      // m_patcher for interpolation from next coarser level to this level
      //      m_patcher.define(a_thisLayout,
      //                       &a_coarserLayout,
      //                       m_dx,
      //                       a_refineCoarse,
      //                       m_numFields,
      //                       a_domain);
      ProblemDomain domainCoarse = coarsen(m_domain, a_refineCoarse);
      // FourthOrderFillPatch m_patcher;
      m_patcher.define(a_thisLayout,
                       a_coarserLayout,
                       m_numFields,
                       domainCoarse,
                       m_refineCoarse,
                       m_numGhost);
    }

  defineFlattening();

  m_doDeconvolution = true;
  m_noPPM = false;
  m_useArtificialViscosity = false;

  // PatchMappedConsOperator m_patchConsOperator;
  m_patchConsOperator.define(m_domain, m_dx, a_gdnvPhysics, m_numFields);
  m_patchConsOperator.spaceOrder(m_spaceOrder);
  m_patchConsOperator.limitFaceValues(m_limitFaceValues);
  m_patchConsOperator.useFlattening(m_useFlattening);
  m_patchConsOperator.noPPM(m_noPPM);
  m_patchConsOperator.doDeconvolution(m_doDeconvolution);
  m_patchConsOperator.doFaceDeconvolution(m_doFaceDeconvolution);
  m_patchConsOperator.useArtificialViscosity(m_useArtificialViscosity);

  // (DFM, 9/29/2005) This is really silly, but the simplest
  // work-around. The problem is that the numConserved and
  // numFluxes functions are not const, but a_gdnvPhysics is.
  // So. the simplest thing is to cast away the const-ness
  // just for this limited context.
  // We may eventually want to make the functions non-const,
  // but that will break all classes derived from GodunovPhysics.
  {
    GodunovPhysics* nonConstPhysicsPtr = (GodunovPhysics*) a_gdnvPhysics;
    m_numCons   = nonConstPhysicsPtr->numConserved();
    m_numFluxes = nonConstPhysicsPtr->numFluxes();
    // Actually we will use SpaceDim * m_numCons as number of fluxes.
  }

  m_exchangeCopier.exchangeDefine(a_thisLayout, m_numGhost*IntVect::Unit);

  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
// Evaluate the operator (-div(F) ) at time a_time,
// on conserved variables a_U:  at coarser level a_UcoarseOld, a_UcoarseNew.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// The flux registers are incremented with the normal derivatives of a_U
// at the boundary, times a fraction of the time step corresponding to
// which stage of an explicit time-stepping scheme (e.g. Runge-Kutta)
// is being invoked.  PATCH BY PATCH.
void LevelMappedConsOperator::evalRHS(
                                LevelData<FArrayBox>&       a_LofU,
                                LevelData<FArrayBox>&       a_U,
                                LevelFluxRegister&          a_finerFluxRegister,
                                LevelFluxRegister&          a_coarserFluxRegister,
                                const LevelData<FArrayBox>& a_UcoarseOld,
                                const Real&                 a_timeCoarseOld,
                                const LevelData<FArrayBox>& a_UcoarseNew,
                                const Real&                 a_timeCoarseNew,
                                Real                        a_time,
                                Real                        a_weight)
{
  // Make sure everything is defined
  CH_assert(m_defined);
  m_evalCount++;

  // Create temporary storage with a layer of "m_numGhost" ghost cells
  IntVect ivGhost = m_numGhost*IntVect::Unit;

  fillGhosts(a_U, a_time, a_timeCoarseOld, a_timeCoarseNew);

  bool setFlattening = (m_useFlattening && (m_evalCount == 1));

  m_patchConsOperator.setCurrentTime(a_time);

  // Obtain <J>
  const LevelData<FArrayBox>& J = m_coordSysPtr->getJ();

  // petermc, 15 Apr 2009
  const LevelData<FluxBox>& faceMetricTerms = m_coordSysPtr->getFaceMetricTerms();

  // Define this here because of need to coordinate with FourthOrderCoordSys.
  // petermc, 7 Apr 2009:  multiplied m_numFluxes by SpaceDim.
  LevelData<FluxBox> FfaceAvg(m_grids, SpaceDim*m_numFluxes, ivGhost);
  // do this fab-by-fab to preserve ghost cells
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& curBox = m_grids[dit];
      m_patchConsOperator.setCurrentBox(curBox);

      m_patchConsOperator.evalRHS(a_LofU[dit],
                                  a_U[dit], // this is <J U>
                                  J[dit], // this is <J>
                                  faceMetricTerms[dit],
                                  m_metricTermComponents,
                                  FfaceAvg[dit],
                                  a_weight,
                                  curBox,
                                  setFlattening,
                                  m_flattening[dit]);
    }

  /*
    Update total fluxes m_fluxes from FfaceAvg:
    m_fluxes[dit] += a_weight * FfaceAvg.

    Update a_finerFluxRegister and a_coarserFluxRegister from FfaceAvg.
  */
  LevelData<FluxBox> normalFfaceAvg(m_grids, m_numFluxes, ivGhost);

  for (int comp = 0; comp < m_numFluxes; comp++)
    {
      // Find one component of a_LofU at a time.
      Interval compInt(comp, comp);
      LevelData<FArrayBox> LofUComp;
      aliasLevelData(LofUComp, &a_LofU, compInt);

      LevelData<FluxBox> FfaceAvgComp(m_grids, SpaceDim, ivGhost);
      // Copy components of FfaceAvg into FfaceAvgComp.
      // These are the components we need for
      // computeMetricTermProductAverage and mappedGridDivergence.

      // Components of FfaceAvg are grouped by dimension, m_numFluxes at a time.
      // FfaceAvgComp has SpaceDim components for a single flux variable.
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FluxBox& FfaceAvgFlub = FfaceAvg[dit];
          FluxBox& FfaceAvgCompFlub = FfaceAvgComp[dit];
          for (int dirFlux = 0; dirFlux < SpaceDim; dirFlux++)
            {
              int srcComp = dirFlux*m_numFluxes + comp;
              int dstComp = dirFlux;
              FfaceAvgCompFlub.copy(FfaceAvgFlub, srcComp, dstComp, 1);
            }
        }

      LevelData<FluxBox> normalFfaceAvgComp(m_grids, 1, ivGhost);
      // You can't alias LevelData<FluxBox>.
      //      aliasLevelData(normalFfaceAvgComp, &normalFfaceAvg, compInt);
      m_coordSysPtr->computeMetricTermProductAverage(normalFfaceAvgComp,
                                                     FfaceAvgComp);

      // Need to call this component by component.
      // LofUComp has 1 component
      // FfaceAvgComp has SpaceDim components
      // m_coordSysPtr->mappedGridDivergence(LofUComp, FfaceAvgComp);
      m_coordSysPtr->simpleDivergence(LofUComp, normalFfaceAvgComp);
      // mappedGridDivergence calls computeMetricTermProductAverage,
      // to find N^T*F.

      // Set normalFfaceAvg[comp] = normalFfaceAvgComp[0] / m_dx .
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FluxBox& normalFfaceAvgCompFlub = normalFfaceAvgComp[dit];
          FluxBox& normalFfaceAvgFlub = normalFfaceAvg[dit];
          normalFfaceAvgFlub.copy(normalFfaceAvgCompFlub, 0, comp, 1);
        }
    }

  // Should this be BEFORE divergence computation?
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& normalFfaceAvgFlub = normalFfaceAvg[dit];
      // kludge that seems to help
      normalFfaceAvgFlub *= 1. / m_dx;
      updateFluxTotalsAndRegisters(normalFfaceAvgFlub,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   dit(), a_weight);
    }

  // actually want -div
  // a_LofU.negate();

  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& LofUFab = a_LofU[dit];
      // Seemed to me that we needed this because I didn't see dx used in
      // mappedGridDivergence.  Actually not sure why dx is squared!
      LofUFab *= -1. / (m_dx * m_dx);
    }

  // added 22 Oct 2008:  these change nothing
  // a_LofU.exchange();
  // a_U.exchange();
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::updateODE(LevelData<FArrayBox>& a_soln,
                             const LevelData<FArrayBox>& a_rhs,
                             Real a_dt)
{
  DataIterator dit = a_soln.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& solnFab = a_soln[dit];
      const FArrayBox& rhsFab = a_rhs[dit];
      m_patchConsOperator.updateODE(solnFab, rhsFab, a_dt);
    }
}


//////////////////////////////////////////////////////////////////////////////
/// reset fluxes contained in this object to zero
void
LevelMappedConsOperator::resetFluxes()
{
  DataIterator dit = m_fluxes.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_fluxes[dit].setVal(0.0);
    }
}


//////////////////////////////////////////////////////////////////////////////
/// returns reference to fluxes computed and accumulated by this operator.
LevelData<FluxBox>&
LevelMappedConsOperator::getFluxes()
{
  return m_fluxes;
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::defineSolnData(LevelData<FArrayBox>& a_newSoln,
                                  const LevelData<FArrayBox>& a_existingSoln)
{
  int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect& ghostVect = a_existingSoln.ghostVect();

  a_newSoln.define(grids, nComp, ghostVect);
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::defineRHSData(LevelData<FArrayBox>& a_newRHS,
                                 const LevelData<FArrayBox>& a_existingSoln)
{
  // same as defineSolnData, but w/o ghost cells
  int nComp = a_existingSoln.nComp();
  const DisjointBoxLayout& grids = a_existingSoln.getBoxes();
  const IntVect ghostVect = IntVect::Zero;

  a_newRHS.define(grids, nComp, ghostVect);
}


//////////////////////////////////////////////////////////////////////////////
/// copy data src->dest
void
LevelMappedConsOperator::copySolnData(LevelData<FArrayBox>& a_dest,
                                const LevelData<FArrayBox>& a_src)
{
  if (a_dest.disjointBoxLayout() == a_src.disjointBoxLayout())
    { // Do this to copy all the ghost cells too.
      for (DataIterator dit = a_src.dataIterator(); dit.ok(); ++dit)
        {
          const FArrayBox& srcFab = a_src[dit];
          FArrayBox& destFab = a_dest[dit];
          destFab.copy(srcFab);
        }
    }
  else
    {
      a_src.copyTo(a_dest);
    }
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
  if (m_defined) m_patchConsOperator.spaceOrder(m_spaceOrder);
}

//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::setCoordSys(FourthOrderCoordSys* a_coordSysPtr)
{
  m_coordSysPtr = a_coordSysPtr;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVect& metricTermComponentsDir = m_metricTermComponents[idir];
      for (int comp = 0; comp < SpaceDim; comp++)
        {
          metricTermComponentsDir[comp] =
            m_coordSysPtr->getMetricTermComponent(idir, comp);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
  if (m_defined) m_patchConsOperator.limitFaceValues(m_limitFaceValues);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
  if (m_defined) m_patchConsOperator.useFlattening(m_useFlattening);
  defineFlattening();
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
  if (m_defined) m_patchConsOperator.noPPM(m_noPPM);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
  if (m_defined) m_patchConsOperator.doDeconvolution(m_doDeconvolution);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
  if (m_defined) m_patchConsOperator.doFaceDeconvolution(m_doFaceDeconvolution);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
  if (m_defined) m_patchConsOperator.useArtificialViscosity(m_useArtificialViscosity);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
  if (m_defined) m_patchConsOperator.artificialViscosity(m_artificialViscosity);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::defineFlattening()
{
  // Need to define m_flattening because we'll take m_flattening[dit].
  //  if (m_useFlattening)
  BoxLayout slopeLayout;
  slopeLayout.deepCopy(m_grids);
  LayoutIterator lit = m_grids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      Box& slopeBox = slopeLayout.ref(lit());
      slopeBox.grow(1);
      slopeBox &= m_domain;
    }
  slopeLayout.close();
  m_flattening.define(slopeLayout, 1);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::resetEvalCount()
{
  m_evalCount = 0;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::evalCountMax(int a_evalCountMax)
{
  m_evalCountMax = a_evalCountMax;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::addArtificialViscosity(LevelData<FArrayBox>&         a_Unew,
                                          const LevelData<FArrayBox>&   a_Uold,
                                          LevelFluxRegister&            a_finerFluxRegister,
                                          LevelFluxRegister&            a_coarserFluxRegister,
                                          Real  a_oldTime,
                                          Real  a_weight)
{
  CH_assert(m_useArtificialViscosity);
  m_patchConsOperator.setCurrentTime(a_oldTime);

  DataIterator dit = a_Uold.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& UnewFab = a_Unew[dit];
      const FArrayBox& UoldFab = a_Uold[dit];

      // The current CELL-centered box
      const Box& curBox = m_grids[dit];
      m_patchConsOperator.setCurrentBox(curBox);

      FluxBox curFlux;
      m_patchConsOperator.addArtificialViscosity(UnewFab,
                                                 UoldFab,
                                                 curFlux,
                                                 curBox,
                                                 a_weight);

      // petermc, 7 Apr 2009:  IGNORE THIS FOR NOW, WITH ONLY ONE LEVEL.
      updateFluxTotalsAndRegisters(curFlux,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   dit(), a_weight);
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::updateFluxTotalsAndRegisters(// used to update flux registers:  it is shifted and shifted back, but data remain unchanged
                                                FluxBox&  a_Fface,
                                                LevelFluxRegister&   a_finerFluxRegister,
                                                LevelFluxRegister&   a_coarserFluxRegister,
                                                const DataIndex&     a_dataIndex,
                                                Real                 a_weight)
{
  // petermc, 7 Apr 2009:  IGNORE THIS FOR NOW, WITH ONLY ONE LEVEL.

  // thisTotalFluxDir += a_weight * curFluxDir
  FluxBox& thisTotalFlux = m_fluxes[a_dataIndex];

  Interval intvlF(0, m_numFluxes-1);
  // Do flux register updates
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // This gets shifted and shifted back, with no data change.
      FArrayBox& curFluxDir = a_Fface[idir];

      // First update m_fluxes:  thisTotalFluxDir += a_weight * curFluxDir
      FArrayBox& thisTotalFluxDir = thisTotalFlux[idir];
      thisTotalFluxDir.plus(curFluxDir, a_weight);

      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner && a_finerFluxRegister.isDefined())
        {
          a_finerFluxRegister.incrementCoarse(curFluxDir, a_weight,
                                              a_dataIndex,
                                              intvlF, intvlF, idir);
        }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser && a_coarserFluxRegister.isDefined())
        {
          a_coarserFluxRegister.incrementFine(curFluxDir, a_weight,
                                              a_dataIndex,
                                              intvlF, intvlF, idir);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

TimeInterpolatorRK4&
LevelMappedConsOperator::getTimeInterpolator()
{
  CH_assert(m_defined);
  return m_patcher.getTimeInterpolator();
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::fillGhosts(LevelData<FArrayBox>&   a_U,
                              const Real&             a_time,
                              const Real&             a_timeCoarseOld,
                              const Real&             a_timeCoarseNew)
{
  if (m_hasCoarser)
    {
      // Check that current fine-level time "a_time" falls between
      // the old and new coarse times.
      Real dtCoarse = a_timeCoarseNew - a_timeCoarseOld;
      Real alpha = (a_time - a_timeCoarseOld) / dtCoarse;

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
          MayDay::Error( "LevelMappedConsOperator::evalRHS: alpha < 0.0");
        }
      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "LevelMappedConsOperator::evalRHS: alpha > 1.0");
        }

      // OK, now we have alpha in the interval [0:1].

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation.

      // We need 4 ghost cell layers of fine grids.
      // To get those, we need 3 ghost cell layers of coarsened fine grids.

      // first, interpolate coarse data to the current time
      //        LevelData<FArrayBox> Ucoarse;
      //        Ucoarse.define(a_UcoarseOld);
      // timeInterp() is a function in BoxTools
      //        timeInterp(Ucoarse       ,a_time,
      //                   a_UcoarseOld  ,a_timeCoarseOld,
      //                   a_UcoarseNew  ,a_timeCoarseNew,
      //                   intvlU);
      // use current-time coarse data to fill fine ghost cells
      // m_patcher.coarseFineInterp(a_U, Ucoarse);
      m_patcher.fillInterp(a_U, alpha, 0, 0, m_numCons);
    }

  // do domain BC's
  // FIX:  In periodic case, this is a no-op.
  // PhysIBC* physIBCPtr = m_gdnvPhysics->getPhysIBC();
  // physIBCPtr->ghostCellBC(a_U, a_time);

  // Exchange all the ghost cell data between grids at this level
  a_U.exchange(m_exchangeCopier);
}

//////////////////////////////////////////////////////////////////////////////

inline void
LevelMappedConsOperator::computeFaceAverages(
                                       LevelData<FluxBox>& a_face_data,
                                       const LevelData<FArrayBox>& a_cell_data ) const
{
  if (m_spaceOrder == 4)
    {
      fourthOrderCellToFace( a_face_data, a_cell_data );
    }
  else if (m_spaceOrder == 2)
    {
      CellToEdge( a_cell_data, a_face_data );
    }
  else
    {
      MayDay::Error("Bad Space Order in LevelMappedConsOperator");
    }
}

//////////////////////////////////////////////////////////////////////////////

inline Real ipow( const Real x, const int p )
{
  Real result = 1.0;
  for (int i=0; i<p; i++) result *= x;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

inline void
LevelMappedConsOperator::getPhysicalCellVolumes(
                                          LevelData<FArrayBox>& a_volumes ) const
{
  a_volumes.define( m_coordSysPtr->getCellVolumes() );
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::computeCompFaceFluxes( LevelData<FluxBox>& a_uTimesV,
                                          const LevelData<FluxBox>& a_u,
                                          const LevelData<FluxBox>& a_v) const
{
   // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
   // space, where a_v is the SpaceDim-dimensional velocity vector and
   // a_u is the nComp-dim state vector
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == SpaceDim);
   CH_assert(a_uTimesV.nComp() == SpaceDim * ncomp);

   // loop over boxes
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
     {
       FluxBox& thisUV = a_uTimesV[dit];
       const FluxBox& thisU = a_u[dit];
       const FluxBox& thisV = a_v[dit];
       m_patchConsOperator.computeCompFaceFluxes(thisUV, thisU, thisV);
     } // end loop over boxes
}

//////////////////////////////////////////////////////////////////////////////
void LevelMappedConsOperator::cellUJToCellU(
                                      LevelData<FArrayBox>& a_Uavg,
                                      const LevelData<FArrayBox>& a_UJavg,
                                      const LevelData<FArrayBox>& a_J) const
{
   // loop over boxes
   DataIterator dit = a_UJavg.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
     {
       FArrayBox& thisU = a_Uavg[dit];
       const FArrayBox& thisUJ = a_UJavg[dit];
       const FArrayBox& thisJ = a_J[dit];
       m_patchConsOperator.cellUJToCellU(thisU, thisUJ, thisJ);
     } // end loop over boxes
}
