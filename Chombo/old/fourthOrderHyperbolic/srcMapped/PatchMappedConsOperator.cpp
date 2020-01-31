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

#include "PatchMappedConsOperator.H"
#include "FourthOrderUtil.H" // functions with names beginning "fourthOrder"
#include "AMRIO.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "AdvectOpF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "UnitNormalsF_F.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
PatchMappedConsOperator::PatchMappedConsOperator()
{
  m_gdnvPhysics = NULL;
  m_isDefined = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet  = false;
  m_dx = 0.0;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
PatchMappedConsOperator::~PatchMappedConsOperator()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void PatchMappedConsOperator::define(
                               const ProblemDomain&      a_domain,
                               const Real&               a_dx,
                               const GodunovPhysics* const a_gdnvPhysics,
                               const int&                a_numFields)
{
  // Sanity checks
  CH_assert(a_dx > 0.0);

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;

  m_numFields = a_numFields;

  // petermc added 6 July 2008; changed from 4 to 5 on 25 Aug.
  // changed from 4 to 5 on 20 Aug 2009
  m_numGhost = 5;

  // GodunovPhysics* m_gdnvPhysics;
  m_gdnvPhysics = a_gdnvPhysics->new_godunovPhysics();
  m_gdnvPhysics->define(m_domain, m_dx);

  m_numFluxes = m_gdnvPhysics->numFluxes();
  // Actually we will use SpaceDim * m_numCons as number of fluxes.

  // GodunovUtilities m_util;
  m_util.define(m_domain, m_dx);
  m_util.highOrderLimiter(true);

  m_doDeconvolution = true;
  m_noPPM = false;
  m_useArtificialViscosity = false;

  // Everything is defined now.
  m_isDefined = true;
}

//////////////////////////////////////////////////////////////////////////////
void PatchMappedConsOperator::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime      = a_currentTime;
  m_isCurrentTimeSet = true;
}

//////////////////////////////////////////////////////////////////////////////
void PatchMappedConsOperator::setCurrentBox(const Box& a_currentBox)
{
  m_currentBox      = a_currentBox;
  m_isCurrentBoxSet = true;

  m_gdnvPhysics->setCurrentBox(a_currentBox);
}

//////////////////////////////////////////////////////////////////////////////
// Evaluate the operator (-div(F) ) at time m_currentTime.
// on conserved variables a_U:  at coarser level a_UcoarseOld, a_UcoarseNew.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// The flux registers are incremented with the normal derivatives of a_U
// at the boundary, times a fraction of the time step corresponding to
// which stage of an explicit time-stepping scheme (e.g. Runge-Kutta)
// is being invoked.
void PatchMappedConsOperator::evalRHS(
                                FArrayBox&          a_LofU,
                                const FArrayBox&    a_JUavgFab,
                                const FArrayBox&    a_JFab,
                                const FluxBox&      a_faceMetricTerms,
                                const Tuple<IntVect, SpaceDim>& a_metricTermComponents,
                                FluxBox&            a_FfaceAvg,
                                Real                a_weight,
                                const Box&          a_gridBox,
                                bool                a_setFlattening,
                                FArrayBox&          a_flatteningFab)
{
  CH_assert(isDefined());
  CH_assert(a_gridBox == m_currentBox);
  /*
    Get UavgFab:  cell-averaged conserved variables
    (from cell-averaged J * U in a_JUavgFab)
  */
  int numU = a_JUavgFab.nComp();
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  Box grownBox = grow(a_gridBox, ivGhost);
  FArrayBox UavgFab(grownBox, numU);
  cellUJToCellU(UavgFab, a_JUavgFab, a_JFab);

  /*
    Get UcellCenFab:  cell-centered conserved variables
    (from cell-averaged conserved variables UavgFab)
  */
  FArrayBox UcellCenFab(grownBox, numU);
  UcellCenFab.copy(UavgFab);
  if (m_doDeconvolution)
    { // UcellCenFab -= laplacian(UavgFab) * m_dx^2 / 24
      // fourthOrderAverageCell(UcellCenFab, -1);
      m_util.deconvolve(UcellCenFab, UavgFab, -1);
    }

  /*
    Get WcellCenFab:  cell-centered primitive variables
    (from cell-centered conserved variables UcellCenFab)
  */
  int numW = m_gdnvPhysics->numConserved();
  FArrayBox WcellCenFab(grownBox, numW);
  m_gdnvPhysics->consToPrim(WcellCenFab, UcellCenFab, grownBox);

  /*
    Get WofUavgFab:  W(UavgFab), where UavgFab is cell-averaged
  */
  FArrayBox WofUavgFab(grownBox, numW);
  m_gdnvPhysics->consToPrim(WofUavgFab, UavgFab, grownBox);

  /*
    Get WcellAvgFab:  cell-averaged primitive variables
    (from cell-centered primitive variables WcellCen)
  */
  FArrayBox WcellAvgFab(grownBox, numW);
  WcellAvgFab.copy(WcellCenFab);
  if (m_doDeconvolution)
    { // WcellAvgFab += laplacian(WcellAvgFab) * m_dx^2 / 24
      // deconvolve(WcellAvgFab, WcellCenFab);
      // 21 Oct 2008:  change this to
      // WcellAvgFab += laplacian(WofUavgFab) * m_dx^2 / 24
      m_util.deconvolve(WcellAvgFab, WofUavgFab);
    }

  // petermc, 15 Apr 2009:
  // we need to find some way to use a_faceMetricTerms
  // to obtain a new fluxbox containing transformations coefficients.
  // it will have SpaceDim*SpaceDim components.
  // i think i might still be able to use the same PolytropicPhysics::riemann
  // function!  A major question is whether primBC needs to be changed.
  // If I use the existing one, it calculates depending on the normal direction.
  FluxBox unitNormals(grownBox, SpaceDim*SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& unitNormalFab = unitNormals[dir];
      const FArrayBox& faceMetricFab = a_faceMetricTerms[dir];
      const Box& bx = unitNormalFab.box();
      const IntVect& normalComponents = a_metricTermComponents[dir];
      // Use normalComponents as component indices into faceMetricFab.
      FORT_GETUNITNORMALS(CHF_FRA(unitNormalFab),
                          CHF_CONST_FRA(faceMetricFab),
                          CHF_CONST_INTVECT(normalComponents),
                          CHF_CONST_INT(dir),
                          CHF_BOX(bx));
    }

  FluxBox WfaceCen(grownBox, numW);
  FluxBox WfaceAvg(grownBox, numW);
  if (m_noPPM)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // WfaceCenDir lives on faces in dimension dir
          FArrayBox& WfaceCenDir = WfaceCen[dir];

          Box faceBox = surroundingNodes(grownBox, dir);
          // leftW and rightW can't be on grownBox:  too small
          FArrayBox leftW(faceBox, numW);
          FArrayBox rightW(faceBox, numW);
          leftW.shiftHalf(dir, -1); // like WPlus
          rightW.shiftHalf(dir, +1); // like WMinus
          leftW.copy(WcellAvgFab);
          rightW.copy(WcellAvgFab);

          const FArrayBox& unitNormalFab = unitNormals[dir];
          // Use normalFab for coordinate transformations.
          // Transform velocity components of leftW and rightW.
          // Get the centerings to be the same as those of unitNormalFab.
          rightW.shiftHalf(dir, -1);
          leftW.shiftHalf(dir, +1);
          forwardBasisTransform(leftW, unitNormalFab);
          forwardBasisTransform(rightW, unitNormalFab);
          // Now shift the centerings back for use in Riemann solver.
          rightW.shiftHalf(dir, +1);
          leftW.shiftHalf(dir, -1);

          // Solve Riemann problem from WcellAvg to get WfaceCen.
          m_gdnvPhysics->riemann(WfaceCenDir, // on FACEs in dir
                                 leftW, // on CELLs
                                 rightW, // on CELLs
                                 WcellAvgFab, // not used in PolytropicPhysics
                                 m_currentTime,
                                 dir,
                                 faceBox); // FACE-centered box

          // Transform back velocity components in WfaceCenDir.
          reverseBasisTransform(WfaceCenDir, unitNormalFab);
        } // end loop over directions
      WfaceAvg.copy(WfaceCen);
    }
  else // use PPM
    {
      /*
        Get WfaceAvg:  face-averaged primitive variables
        (from cell-averaged primitive variables WcellAvg)
      */

      // If we're limiting the face-centered values, do it here.
      // Includes Riemann solver.
      getFaceAvg(WfaceAvg, WcellAvgFab, WofUavgFab, unitNormals,
                 a_gridBox, a_flatteningFab, a_setFlattening);

      /*
        Get WfaceCen:  face-centered primitive variables
        (from face-averaged primitive variables WfaceAvg)
      */
      WfaceCen.copy(WfaceAvg);
      // fourthOrderAverageFace(WfaceCenFlux, -1);
      // 26 Sep 2008
      m_util.deconvolveFace(WfaceCen, WfaceAvg, -1);
    } // whether or not to use PPM

  /*
    Get a_FfaceAvg:  4th-order flux, face-averaged F(W)
    (from primitive variables:  face-averaged WfaceAvg, face-centered WfaceCen)
  */
  // This is defined in LevelConsOperator because of problems with
  // coordinating with FourthOrderCoordSys.
  //  a_FfaceAvg.resize(grownBox, m_numFluxes);
  FluxBox FfaceCen(grownBox, SpaceDim*m_numFluxes);
  // Box nullBox;

  for (int dirFace=0; dirFace<SpaceDim; dirFace++)
    { // dirFace specifies which faces:  either x or y or z
      // box of valid edges for this grid
      // Box faceBox(a_gridBox);
      Box faceBox(grownBox);
      faceBox.surroundingNodes(dirFace);

      const FArrayBox& WfaceAvgDir = WfaceAvg[dirFace];
      const FArrayBox& WfaceCenDir = WfaceCen[dirFace];

      FArrayBox& FfaceAvgDir = a_FfaceAvg[dirFace];
      FArrayBox& FfaceCenDir = FfaceCen[dirFace];

      // Fill in FfaceCenDir and FfaceAvgDir.
      // Each has SpaceDim*m_numFluxes components.
      // We fill in m_numFluxes components at a time.
      // So they will contain:
      // 0:m_numFluxes-1 for dirFlux=0
      // ...
      // (SpaceDim-1)*m_numFluxes:SpaceDim*m_numFluxes-1 for dirFlux=SpaceDim-1
      for (int dirFlux=0; dirFlux<SpaceDim; dirFlux++)
        {
          int fluxIntLo = dirFlux * m_numFluxes;
          int fluxIntHi = fluxIntLo + m_numFluxes-1;
          Interval fluxInt(fluxIntLo, fluxIntHi);
          // FfaceCenDirD = FfaceCenDir[fluxInt]
          // FfaceAvgDirD = FfaceAvgDir[fluxInt]
          FArrayBox FfaceCenDirD(fluxInt, FfaceCenDir);
          FArrayBox FfaceAvgDirD(fluxInt, FfaceAvgDir);
          m_gdnvPhysics->getFlux(FfaceCenDirD, WfaceCenDir, dirFlux, faceBox);
          m_gdnvPhysics->getFlux(FfaceAvgDirD, WfaceAvgDir, dirFlux, faceBox);
        }
    }

  if (m_doFaceDeconvolution)
    { // Set FfaceAvg = FfaceCen + h^2/24 * laplacian(FfaceAvg).
      // That is, set <F> = F(W) + h^2/24 * laplacian(F(<W>)).
      // The subroutine m_gdnvPhysics->getFlux() performs the function F.
      // petermc, 3 Oct 2008, new:
      FluxBox FfaceAvgCopy(grownBox, SpaceDim*m_numFluxes);
      FfaceAvgCopy.copy(a_FfaceAvg);
      a_FfaceAvg.copy(FfaceCen);
      // This should still work with SpaceDim*m_numFluxes components.
      m_util.deconvolveFace(a_FfaceAvg, FfaceAvgCopy);
    }

  /*
    Set a_LofU = -div(FfaceAvg).
  */

  // We don't have this for FABs.  So find it in the calling function.
  // mappedGridDivergence(a_LofU, a_FfaceAvg);

  // actually want -div
  // a_LofU.negate();
  // Seems to me that we need this because I don't see dx used in
  // mappedGridDivergence.
  //  a_LofU *= -1. / m_dx;
}


//////////////////////////////////////////////////////////////////////////////
void
PatchMappedConsOperator::updateODE(FArrayBox&        a_solnFab,
                             const FArrayBox&  a_rhsFab,
                             Real a_dt)
{
  int rhsComp = a_rhsFab.nComp();
  CH_assert(a_solnFab.nComp() == rhsComp);

  const Box& bx = a_rhsFab.box();
  a_solnFab.plus(a_rhsFab, bx, bx, a_dt, 0, 0, rhsComp);
}


//////////////////////////////////////////////////////////////////////////////
void
PatchMappedConsOperator::getFaceAvg(// we'll apply limiter to a_faceW
                                    FluxBox& a_faceW,
                                    const FArrayBox& a_cellW,
                                    const FArrayBox& a_WofUavg,
                                    const FluxBox& a_unitNormals,
                                    // CELL-centered base box
                                    const Box& a_gridBox,
                                    FArrayBox& a_flatteningFab,
                                    bool a_setFlattening)
{
  if (m_spaceOrder == 4)
    { // This uses the coefficients 7/12, -1/12, -1/12, 7/12.
      fourthOrderCellToFace(a_faceW, a_cellW);
    }
  else if (m_spaceOrder == 2)
    {
      // regular 2nd-order averaging from WcellAvg to WfaceAvg
      CellToEdge(a_cellW, a_faceW);
    }
  else
    {
      MayDay::Error("Bad Space Order in PatchMappedConsOperator");
    }

  CH_assert(isDefined());
  CH_assert(a_gridBox == m_currentBox);

  int nComp = a_faceW.nComp();

  // need two laplacian ghost cells for this limiter
  // 23 Sep 2008:  changed from 2 to 3
  Box lapBox = grow(a_gridBox, 3);

  // check to be sure that we have enough W to compute
  // all of these Laplacians
  Box Wbox = a_cellW.box();
  {
    Box LapWBox = grow(lapBox, 1);
    CH_assert(Wbox.contains(LapWBox));
  }
  Wbox &= m_domain;

  // Define the box of cells where slopes will be needed
  // (one larger than the final update box)
  Box slopeBox = grow(a_gridBox, 1);
  slopeBox &= m_domain;

  // added by petermc, 22 Oct 2008, box to send to PPMLimiter
  //  Box slope2Box = grow(a_gridBox, 2);
  //  slope2Box &= m_domain;

  // flattening stencil has size 3
  // CH_assert(Wbox.contains(grow(slopeBox, 3))); // 27 Aug 2008 removed
  // Compute flattening once for all slopes if needed.
  // Each CELL has a flattening coefficient, eta.
  // Flattening means replacing Wface on left and right of each cell
  // by eta*Wface + (1-eta)*Wcell.
  if (m_useFlattening && a_setFlattening)
    {
      Interval velInt    = m_gdnvPhysics->velocityInterval();
      int pressureInd    = m_gdnvPhysics->pressureIndex();
      Real smallPressure = m_gdnvPhysics->smallPressure();
      int bulkInd        = m_gdnvPhysics->bulkModulusIndex();
      // 21 Oct 2008:  changed from a_cellW to a_WofUavg.
      // We calculate a_flatteningFab on cells of base box grown by 1,
      // and we need a_WofUavg on cells of base box grown by 4.
      m_util.computeFlattening(a_flatteningFab, a_WofUavg,
                               velInt, pressureInd, smallPressure, bulkInd,
                               slopeBox);
    }

  for (int dir = 0; dir < SpaceDim; dir++)
    {
      Box faceBox = surroundingNodes(a_gridBox, dir);
      // added 12 Sep 2008
      // 22 Oct 2008:  changed from 1 to 2
      faceBox.grow(2);
      // 11 Aug 2009:  needed for diffW in COLELLASEKORALIMITERF
      // 20 Aug 2009:  changed from 1 to 2
      faceBox.grow(2 * BASISV(dir));
      // 1 Jul 2009:  expand faceBox by 2 in transverse directions
      // 20 Aug 2009:  changed from 2 to m_numGhost - 2
      faceBox.grow((m_numGhost - 2)*(IntVect::Unit - BASISV(dir)));

      // thisFaceWDir = a_faceW[dit][dir] lives on FACEs in dimension dir
      FArrayBox& thisFaceWDir = a_faceW[dir];

      m_util.PPMFaceValues(thisFaceWDir, a_cellW, nComp, m_limitFaceValues,
                           dir, faceBox, m_currentTime, m_gdnvPhysics);

      if (m_limitFaceValues || m_useFlattening)
        {
          // WAS both of these on slopeBox
          FArrayBox WMinus(faceBox, nComp); // on FACEs
          FArrayBox WPlus (faceBox, nComp); // on FACEs
          // Riemann solver makes the opposite shifts
          WPlus .shiftHalf(dir, -1); // on CELLs of a_gridBox grown by 2, including shift down by 1 in direction dir
          WMinus.shiftHalf(dir, +1); // on CELLs of a_gridBox grown by 2, including shift up by 1 in direction dir

          WMinus.setVal(0.);
          WPlus .setVal(0.);

          WMinus -= a_cellW;
          WPlus  -= a_cellW;

          thisFaceWDir.shiftHalf(dir, +1); // [i-e/2] becomes [i]
          WMinus += thisFaceWDir; // WMinus[i] = thisFaceWDir[i-e/2] - a_cellW[i]
          thisFaceWDir.shift(dir, -1); // [i+e/2] becomes [i]
          WPlus  += thisFaceWDir; // WPlus[i] = thisFaceWDir[i+e/2] - a_cellW[i]
          thisFaceWDir.shiftHalf(dir, +1); // final shift back

          // We modify/need WMinus on cells of base box including shift up by 1 in direction dir, and grown by 1 in the other directions,
          // and we modify/need WPlus on cells of base box including shift down by 1 in direction dir, and grown by 1 in the other directions,
          // and we need a_cellW on cells of base box grown by 3 in direction dir, and grown by 1 in the other directions.
          m_util.PPMLimiter(WMinus, WPlus, a_cellW, nComp, dir, slopeBox);

          if (m_useFlattening)
            { // Multiply WMinus and WPlus by a_flatteningFab on slopeBox.
              // We need a_flatteningFab on cells of base box grown by 1 in all directions.
              m_util.applyFlattening(WMinus, a_flatteningFab, slopeBox);
              m_util.applyFlattening(WPlus,  a_flatteningFab, slopeBox);
            }

          // WMinus[i] = thisFaceWDir[i-e/2] on right after limiting
          WMinus += a_cellW;
          // WPlus[i] = thisFaceWDir[i+e/2] on left after limiting
          WPlus  += a_cellW;

          // petermc, 15 Apr 2009:  this is where we transform WMinus and WPlus
          const FArrayBox& unitNormalFab = a_unitNormals[dir];
          // Transform velocity components of leftW and rightW.
          // Get the centerings to be the same as those of unitNormalFab.
          WMinus.shiftHalf(dir, -1);
          WPlus.shiftHalf(dir, +1);
          forwardBasisTransform(WPlus, unitNormalFab);
          forwardBasisTransform(WMinus, unitNormalFab);
          // Now shift the centerings back for use in Riemann solver.
          WMinus.shiftHalf(dir, +1);
          WPlus.shiftHalf(dir, -1);

          // Now solve Riemann problem.
          // Modify thisFaceWDir[i+e/2] based on
          // WPlus[i] = Wleft[i+e/2] and
          // WMinus[i+e] = Wright[i+e/2].
          // We calculate thisFaceWDir on faces of base box in direction dir grown by 1 in the other directions ONLY,
          // and we need WPlus on cells of base box including shift down by 1 in direction dir, and grown by 1 in the other directions,
          // and we need WMinus on cells of base box including shift up by 1 in direction dir, and grown by 1 in the other directions.
          m_gdnvPhysics->riemann(thisFaceWDir, // on FACEs in dir
                                 WPlus, // left, on CELLs
                                 WMinus, // right, on CELLs
                                 a_cellW, // not used in PolytropicPhysics
                                 m_currentTime,
                                 dir,
                                 faceBox); // on FACEs

          // Transform back velocity components in thisFaceWDir.
          // The centering is OK!
          reverseBasisTransform(thisFaceWDir, unitNormalFab);
        }
      else
        { // neither limiting nor flattening,
          // therefore call Riemann solver only on physical boundary.
          PhysIBC* bc = m_gdnvPhysics->getPhysIBC();
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();

              Box boundaryBox;
              bc->getBoundaryFaces(boundaryBox, faceBox, dir, side);
              if (! boundaryBox.isEmpty() )
                {
                  FArrayBox WPlusOrMinus(boundaryBox, nComp); // on FACEs
                  WPlusOrMinus.copy(thisFaceWDir, boundaryBox);

                  int signSide = sign(side);
                  const FArrayBox& unitNormalFab = a_unitNormals[dir];
                  // Transform velocity components.
                  // Get the centerings to be the same as those of unitNormalFab.
                  WPlusOrMinus.shiftHalf(dir, signSide);
                  forwardBasisTransform(WPlusOrMinus, unitNormalFab);
                  // Now shift the centerings back for use in Riemann solver.
                  WPlusOrMinus.shiftHalf(dir, -signSide);
                  // This will fill in thisFaceWDir on boundaryBox only:
                  // it recalculates boundaryBox from thisFaceWDir.box().
                  bc->primBC(thisFaceWDir,
                             WPlusOrMinus, a_cellW, dir, side, m_currentTime);
                  // Transform back velocity components in thisFaceWDir
                  // (actually only on boundaryBox:  take WPlusOrMinus).
                  // The centering is OK!
                  WPlusOrMinus.copy(thisFaceWDir);
                  reverseBasisTransform(WPlusOrMinus, unitNormalFab);
                  thisFaceWDir.copy(WPlusOrMinus);
                }
            }
        }
    } // end loop over directions
}


//////////////////////////////////////////////////////////////////////////////
void
PatchMappedConsOperator::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}


//////////////////////////////////////////////////////////////////////////////
void
PatchMappedConsOperator::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::addArtificialViscosity(FArrayBox&         a_UnewFab,
                                          const FArrayBox&   a_UoldFab,
                                          FluxBox&           a_flux,
                                          const Box&         a_gridBox,
                                          Real  a_weight)
{
  CH_assert(isDefined());
  CH_assert(m_useArtificialViscosity);
  CH_assert(a_gridBox == m_currentBox);
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  Box grownBox = grow(a_gridBox, ivGhost);
  a_flux.resize(grownBox, m_numFluxes);
  FArrayBox divFlux(a_UnewFab.box(), a_UnewFab.nComp());
  divFlux.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& curFluxDir = a_flux[idir];
      curFluxDir.setVal(0.);
      m_gdnvPhysics->artVisc(curFluxDir,
                             a_UoldFab,
                             m_artificialViscosity,
                             m_currentTime,
                             idir,
                             a_gridBox);

      // increment divergence (divFlux) by derivative in direction idir
      FORT_DIVERGENCE(CHF_CONST_FRA(curFluxDir),
                      CHF_FRA(divFlux),
                      CHF_BOX(a_gridBox),
                      CHF_CONST_REAL(m_dx),
                      CHF_INT(idir));
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
  a_UnewFab.plus(divFlux, -a_weight);

  // LevelConsOperator needs to increment m_fluxes and flux registers.
}

//////////////////////////////////////////////////////////////////////////////
bool PatchMappedConsOperator::isDefined() const
{
  return m_isDefined        &&
         m_isCurrentTimeSet &&
         m_isCurrentBoxSet  ;
}


//////////////////////////////////////////////////////////////////////////////
void PatchMappedConsOperator::cellUJToCellU(
                                      FArrayBox& a_UavgFab,
                                      const FArrayBox& a_UJavgFab,
                                      const FArrayBox& a_JFab) const
{
  // return <U> = <UJ>/<J> - h^2/24 * (grad (<UJ>/<J>)) dot (grad <J>) / <J>
  // as described in ESL paper:  Colella, Dorr, Hittinger, and Martin,
  // "High-order, finite-volume methods in mapped coordinates"
  int ncomp = a_UavgFab.nComp();
  // Set JinvFab = 1.0 / a_Jfab .
  FArrayBox JinvFab(a_JFab.box(), 1);
  JinvFab.copy(a_JFab);
  JinvFab.invert(1.0);

  Box intersectBox(a_UJavgFab.box());
  intersectBox &= JinvFab.box();
  intersectBox &= a_UavgFab.box();

  // Initially set <U> = <UJ>/<J> on intersectBox.
  a_UavgFab.copy(a_UJavgFab, intersectBox);
  for (int comp = 0; comp < ncomp; comp++)
    {
      // Multiply a_UavgFab[comp] by JinvFab[0].
      a_UavgFab.mult(JinvFab, intersectBox, 0, comp);
    }

  // compute box over which we can do this
  Box gradIntersectBox(a_UavgFab.box());
  gradIntersectBox.grow(-1);
  gradIntersectBox &= a_JFab.box();
  // Set gradProduct = factor * grad(a_UavgFAb) dot grad(a_JFab).
  FArrayBox gradProduct(gradIntersectBox, ncomp);
  gradProduct.setVal(0.);
  // since dx's cancel, use 1 here
  Real fakeDx = 1.0;
  Real factor = -1.0/12.0; // (dfm - 3/22/09)--WAS 12 instead of 24 -- No, change back!
  for (int dir=0; dir<SpaceDim; dir++)
    {
      // This function is in fourthOrderMappedGrids/util/FourthOrderUtil.cpp
      // sets gradProduct += factor * d(a_UavgFab)/dx[dir] * d(a_JFab)/dx[dir]
      incrementGradProduct(gradProduct,
                           a_UavgFab,
                           a_JFab,
                           gradIntersectBox,
                           fakeDx,
                           factor,
                           dir);
    }

  for (int comp = 0; comp < ncomp; comp++)
    {
      // Multiply gradProduct[comp] by JinvFab[0].
      gradProduct.mult(JinvFab, gradIntersectBox, 0, comp);
    }
  // Now gradProduct = factor * grad(a_UavgFAb) dot grad(a_JFab) / JFab.

  a_UavgFab.plus(gradProduct, gradIntersectBox, 0, 0, ncomp);
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::computeCompFaceFluxes( FluxBox& a_uTimesV,
                                          const FluxBox& a_u,
                                          const FluxBox& a_v) const
{
   // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
   // space, where a_v is the SpaceDim-dimensional velocity vector and
   // a_u is the nComp-dim state vector
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == SpaceDim);
   CH_assert(a_uTimesV.nComp() == SpaceDim * ncomp);

   // loop over faces (index "d" in notes)
   for (int faceDir=0; faceDir<SpaceDim; faceDir++)
     {
       FArrayBox& thisUVdir = a_uTimesV[faceDir];
       const FArrayBox& thisUdir = a_u[faceDir];
       const FArrayBox& thisVdir = a_v[faceDir];

       // compute <u_p><v_d> tensor
       Box intersectBox(thisUdir.box());
       intersectBox &= thisVdir.box();
       intersectBox &= thisUVdir.box();

       thisUVdir.setVal(0.0);
       FORT_INCREMENTFACEPROD(CHF_FRA(thisUVdir),
                              CHF_CONST_FRA(thisUdir),
                              CHF_CONST_FRA(thisVdir),
                              CHF_BOX(intersectBox));

       if (m_spaceOrder == 4)
         {
           // now increment with product of tangential gradients
           // (sum over index "d'" in notes)
           for (int tanDir=0; tanDir<SpaceDim; tanDir++)
             {
               if (tanDir!=faceDir)
                 {
                   // since dx's cancel, use 1 here
                   Real fakeDx = 1.0;
                   Real factor = 1.0/12.0;

                   // compute box over which we can do this
                   Box gradIntersectBox(thisUdir.box());
                   gradIntersectBox &= thisVdir.box();
                   gradIntersectBox.grow(tanDir,-1);
                   gradIntersectBox &= thisUVdir.box();

                   FORT_INCREMENTFACEPRODGRAD(CHF_FRA(thisUVdir),
                                              CHF_CONST_FRA(thisUdir),
                                              CHF_CONST_FRA(thisVdir),
                                              CHF_BOX(gradIntersectBox),
                                              CHF_REAL(fakeDx),
                                              CHF_REAL(factor),
                                              CHF_INT(tanDir));

                 } // end if tangential direction
             } // end loop over tanDir
         } // end if fourth order
     } // end loop over faceDir
}

//////////////////////////////////////////////////////////////////////////////
void
PatchMappedConsOperator::forwardBasisTransform(FArrayBox& a_W,
                                         const FArrayBox& a_unitNormalFab) const
{
  const Box& bx = a_W.box();
  Interval velInt = m_gdnvPhysics->velocityInterval();
  FArrayBox velFab(velInt, a_W);
  FORT_FORWARDTRANSFORMF(CHF_FRA(velFab),
                         CHF_CONST_FRA(a_unitNormalFab),
                         CHF_BOX(bx));
}

//////////////////////////////////////////////////////////////////////////////
void
PatchMappedConsOperator::reverseBasisTransform(FArrayBox& a_W,
                                         const FArrayBox& a_unitNormalFab) const
{
  const Box& bx = a_W.box();
  Interval velInt = m_gdnvPhysics->velocityInterval();
  FArrayBox velFab(velInt, a_W);
  FORT_REVERSETRANSFORMF(CHF_FRA(velFab),
                         CHF_CONST_FRA(a_unitNormalFab),
                         CHF_BOX(bx));
}

//////////////////////////////////////////////////////////////////////////////
// Write these functions:
// computeCompFaceFluxes
// computeFaceAreas computeFaceGradU computeNJInvGradUOnBox
// computeNTNJInvGradUOnBox computeNormalDiffusiveFluxVector
// secondOrderMappedGridDivergence divideOutCellVolume secondOrderLaplacian
