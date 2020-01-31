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
#include "LayoutData.H"
#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "UnitNormalsF_F.H"
#include "LoHiCenter.H"
#include "PhysMappedIBC.H"
#include "HighDerivsF_F.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////

// Constructor - set up some defaults
PatchMappedConsOperator::PatchMappedConsOperator(const int a_numFluxesPerField)
{
  m_molMappedPhysics = NULL;
  m_isDefined = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet  = false;
  m_isMaxSpeedSet = false;
  m_maxSpeed = 0.;
  m_dx = 0.0;
  m_numFluxesPerField = a_numFluxesPerField;
}

//////////////////////////////////////////////////////////////////////////////

// Destructor - free up storage
PatchMappedConsOperator::~PatchMappedConsOperator()
{
  if (m_molMappedPhysics != NULL)
    {
      delete m_molMappedPhysics;
    }
}

//////////////////////////////////////////////////////////////////////////////

// Define the object so that time stepping can begin
void
PatchMappedConsOperator::define(const ProblemDomain&    a_domain,
                                const Real&             a_dx,
                                LevelGridMetrics*       a_levelGridMetrics,
                                const MOLPhysics* const a_molPhysics, // MOLMappedPhysics
                                const int&              a_numFields)
{
  CH_TIME("PatchMappedConsOperator::define");
  // Sanity checks
  CH_assert(a_dx > 0.0);

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;

  m_numFields = a_numFields;

  // Free up old m_molMappedPhysics, if any.
  if (m_molMappedPhysics != NULL)
    {
      delete m_molMappedPhysics;
    }
  // We should have been passed in a MOLMappedPhysics-derived instance
  CH_assert(a_molPhysics != NULL);
  m_molMappedPhysics = 
    dynamic_cast<MOLMappedPhysics*>(a_molPhysics->new_molPhysics());
  CH_assert(m_molMappedPhysics != NULL);

  m_molMappedPhysics->define(m_domain, m_dx);

  m_numFluxes = m_molMappedPhysics->numFluxes();

  m_levelGridMetrics = a_levelGridMetrics;
  // NOTE: We look at a_coordSys, but we do not assume control of it.
  m_coordSys = &(m_levelGridMetrics->getCoordSys());

  // GodunovUtilities m_util;
  m_util.define(m_domain, m_dx);
  m_highOrderLimiter = true;

  m_doDeconvolution = true;
  m_noPPM = false;
  m_useArtificialViscosity = false;
  m_useArtificialDissipation = false;
  m_artificialDissipation = 0.;

  // Everything is defined now.
  m_isDefined = true;
}

//////////////////////////////////////////////////////////////////////////////

void PatchMappedConsOperator::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime      = a_currentTime;
  m_isCurrentTimeSet = true;

  m_molMappedPhysics->setCurrentTime(a_currentTime);
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::setCurrentBox(const DataIndex& a_currentDataIdx,
                                 const Box&       a_currentBox)
{
  CH_TIME("PatchMappedConsOperator::setCurrentBox");
  m_currentDataIdx  = a_currentDataIdx;
  m_currentBox      = a_currentBox;
  m_isCurrentBoxSet = true;

  // Construct a box that covers the current block, plus any needed ghost cells.
  int blockNum = m_coordSys->whichBlock(a_currentBox);
  const Box& blockBox = m_coordSys->mappingBlocks()[blockNum];
  // m_blockDomain = m_coordSys->mappingBlocks()[blockNum];
  m_blockDomain = m_levelGridMetrics->blockDomain(blockBox, 5);
  m_molMappedPhysics->define(m_blockDomain, m_dx);

  // Set the current Box and coordinate system on the physics object.
  m_molMappedPhysics->setCurrentBox(m_currentBox);
  m_molMappedPhysics->setCurrentCoordSys(m_coordSys->getCoordSys(blockNum));
  m_molMappedPhysics->setUnitNormalBasisPtr(&((*m_unitNormalLayPtr)[m_currentDataIdx]));

  m_util.define(m_blockDomain, m_dx);
  m_util.highOrderLimiter(m_highOrderLimiter);
}

//////////////////////////////////////////////////////////////////////////////

// Evaluate the operator F at time m_currentTime.
void
PatchMappedConsOperator::evalFlux(const FArrayBox& a_UavgFab,
                            FluxBox&         a_FfaceAvg,
                            FluxBox&         a_FfromWfaceAvg,
                            Real             a_weight,
                            bool             a_setFlattening,
                            FArrayBox&       a_flatteningFab)
{
  CH_TIME("PatchMappedConsOperator::evalFlux");
  CH_assert(isDefined());
  CH_assert(m_isCurrentBoxSet);
  int numU = a_UavgFab.nComp();
  CH_assert(m_numGhost >= 5);

  /*
    Set UcellCenFab = a_UavgFab - m_dx^2/24 * Laplacian(a_UavgFab)
    on all cells of bx4inDomain = grow(m_currentBox, 4) & m_domain
    using a_UavgFab
    on all cells of bx5inDomain = grow(m_currentBox, 5) & m_domain.
  */
  Box bx5inDomain = grow(m_currentBox, 5);
  bx5inDomain &= m_blockDomain;
  Box bx4inDomain = grow(m_currentBox, 4);
  bx4inDomain &= m_blockDomain;
  FArrayBox UcellCenFab(bx4inDomain, numU);
  UcellCenFab.copy(a_UavgFab);
  if (m_doDeconvolution)
    { // UcellCenFab -= Laplacian(a_UavgFab) * m_dx^2 / 24
      // fourthOrderAverageCell(UcellCenFab, -1);
      m_util.deconvolve(UcellCenFab, a_UavgFab, bx4inDomain, -1);
    }

  /*
    Set WcellCenFab = W(UcellCenFab)
    on all cells of bx4inDomain = grow(m_currentBox, 4) & m_blockDomain
    using UcellCenFab
    on all cells of bx4inDomain.
  */
  int numW = m_molMappedPhysics->numPrimitives();
  FArrayBox WcellCenFab(bx4inDomain, numW);
  m_molMappedPhysics->consToPrim(WcellCenFab, UcellCenFab, bx4inDomain);

  /*
    Set WofUavgFab = W(a_UavgFab)
    on all cells of bx5inDomain = grow(m_currentBox, 5) & m_blockDomain
    using a_UavgFab
    on all cells of bx5inDomain.
  */
  FArrayBox WofUavgFab(bx5inDomain, numW);
  m_molMappedPhysics->consToPrim(WofUavgFab, a_UavgFab, bx5inDomain);

  /*
    Set WcellAvgFab = WofUavgFab + m_dx^2/24 * Laplacian(WofUavgFab)
    on all cells of bx4inDomain = grow(m_currentBox, 4) & m_blockDomain
    using WofUavgFab
    on all cells of bx5inDomain = grow(m_currentBox, 5) & m_blockDomain.
  */
  FArrayBox WcellAvgFab(bx4inDomain, numW);
  WcellAvgFab.copy(WcellCenFab);
  if (m_doDeconvolution)
    {
      // deconvolve(WcellAvgFab, WcellCenFab);
      m_util.deconvolve(WcellAvgFab, WofUavgFab, bx4inDomain);
    }

  Box bx1 = grow(m_currentBox, 1);
  Box bx1inDomain(bx1);
  bx1inDomain &= m_blockDomain;

  FluxBox WfaceAvg(bx1inDomain, numW);
  /*
    Set WfaceAvg
    on all faces of bx1inDomain = grow(m_currentBox, 1) & m_blockDomain,
    using WcellAvgFab
    on bx4inDomain = grow(m_currentBox, 4) & m_blockDomain,
    and WofUavgFab
    on bx4inDomain = grow(m_currentBox, 4) & m_blockDomain,
    and a_flatteningFab
    on bx1inDomain.

    Fourth-order interpolation:
    WfaceAvg[i+e/2] = 7/12 * (WcellAvgFab[i] + WcellAvgFab[i+e])
                    - 1/12 * (WcellAvgFab[i-e] + WcellAvgFab[i+2e]).

    Includes calls to limiter (if used) and flattening (if used),
    and Riemann solver.
  */
  getFaceAvg(WfaceAvg, WcellAvgFab, WofUavgFab,
             a_flatteningFab, a_setFlattening);

  /*
    Set WfaceCen = WfaceAvg - m_dx^2/24 * Laplacian(WfaceAvg)
    on all faces of m_currentBox
    using WfaceAvg
    on all faces of bx1inDomain = grow(m_currentBox, 1) & m_blockDomain.
  */
  FluxBox WfaceCen(m_currentBox, numW);
  // First set WfaceCen = WfaceAvg.
  WfaceCen.copy(WfaceAvg, WfaceCen.box());
  // Then set WfaceCen -= Laplacian(WfaceAvg) * m_dx^2 / 24
  // on all faces of m_currentBox
  // using WfaceAvg
  // on all faces of m_currentBox.
  // Special stencils next to boundary of m_blockDomain.
  m_util.deconvolveFace(WfaceCen, WfaceAvg, m_currentBox, -1);

  /*
    Set a_FfromWfaceAvg = F(WfaceAvg)
    on all faces of bx1inDomain = grow(m_currentBox, 1) & m_blockDomain
    using WfaceAvg
    on all faces of bx1inDomain.

    Set FfaceCen = F(WfaceCen)
    on all faces of m_currentBox
    using WfaceCen
    on all faces of m_currentBox.
  */
  FluxBox FfaceCen(m_currentBox, m_numFluxesPerField * m_numFluxes);
  getAllFluxes(a_FfromWfaceAvg, FfaceCen, WfaceAvg, WfaceCen);
  if (m_doFaceDeconvolution)
    {
      /*
        Set a_FfaceAvg = FfaceCen + m_dx^2/24 * Laplacian(FfromWfaceAvg)
        on all idir-faces of m_currentBox
        using FfaceCen
        on all idir-faces of m_currentBox
        and original a_FfaceAvg
        on all idir-faces of bx1inDomain = grow(m_currentBox, 1) & m_blockDomain.
      */
      // First set a_FfaceAvg = FfaceCen = flux(WfaceCen)
      // on all faces of m_currentBox.
      a_FfaceAvg.copy(FfaceCen, FfaceCen.box());
      // Then set a_FfaceAvg += Laplacian(FfromWfaceAvg) * m_dx^2 / 24
      // on all faces of m_currentBox,
      // using FfromWfaceAvg
      // on all faces of bx1inDomain = grow(m_currentBox, 1) & m_blockDomain.
      // Special stencils next to boundary of m_blockDomain.
      m_util.deconvolveFace(a_FfaceAvg, a_FfromWfaceAvg, m_currentBox);
    }
  else
    {
      a_FfaceAvg.copy(a_FfromWfaceAvg);
    }

  if (m_useArtificialDissipation)
    {
      // Add fifth-derivative correction to flux
      fluxCorrection(a_FfaceAvg, a_UavgFab);
    }

  /*
    Set a_LofU = -div(a_FfaceAvg)
    on all cells of m_currentBox,
    using a_FfaceAvg on all faces of m_currentBox.
    Derived class may also use FfaceCen.
  */
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////

// Evaluate the operator (-div(F) ) at time m_currentTime.
void PatchMappedConsOperator::evalRHS(FArrayBox&          a_LofU,
                                const FluxBox&      a_FfaceAvg)
{
  CH_TIME("PatchMappedConsOperator::evalRHS");
  // Calculate div(F)
  PatchConsOperator::getFluxDivergence(a_LofU, a_FfaceAvg);
  // Actually want -div
  a_LofU.negate();
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::getAllFluxes(FluxBox&        a_FfaceAvg,  // On box+1
                                FluxBox&        a_FfaceCen,  // On box
                                const FluxBox&  a_WfaceAvg,
                                const FluxBox&  a_WfaceCen)
{
  CH_TIME("PatchMappedConsOperator::getAllFluxes");
  // Need to cast non-const in order to alias
//   FluxBox& WfaceAvg = (FluxBox&) a_WfaceAvg;
//   FluxBox& WfaceCen = (FluxBox&) a_WfaceCen;
  CH_assert(m_isCurrentBoxSet);
  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_blockDomain;
  // Sanity checks
  CH_assert(a_FfaceAvg.box().contains(bx1inDomain));
  CH_assert(a_FfaceCen.box().contains(m_currentBox));

  // T - we need the full flux dyad on each face
  // F - we only need the flux in the direction orthogonal to the face
  const bool allFDir = (a_FfaceAvg.nComp() == SpaceDim*m_numFluxes);

  for (int idir = 0; idir < SpaceDim; idir++)
    { // idir specifies which faces:  either x or y or z
      Box faceBox1(bx1inDomain);
      faceBox1.surroundingNodes(idir);

      const FArrayBox& WfaceAvgDir = a_WfaceAvg[idir];
      FArrayBox& FfaceAvgDir = a_FfaceAvg[idir];
      if (allFDir)
        {
          for (int jdir = 0; jdir < SpaceDim; jdir++)
            {
              // Interval
              int fluxIntLo = jdir*m_numFluxes;
              int fluxIntHi = fluxIntLo + m_numFluxes-1;
              Interval fluxIntv(fluxIntLo, fluxIntHi);
              // Alias to pretend the fluxbox has only m_fluxes components
              FArrayBox FfaceAvgDirD(fluxIntv, FfaceAvgDir);
              m_molMappedPhysics->getFlux(FfaceAvgDirD,
                                          WfaceAvgDir,
                                          jdir,
                                          faceBox1);
            }
        }
      else
        {
          // Set FfaceAvgDir = flux(WfaceAvgDir) on all idir-faces of
          // bx1inDomain = grow(m_currentBox, 1) & m_blockDomain.
          m_molMappedPhysics->getFlux(FfaceAvgDir, WfaceAvgDir, idir, faceBox1);
        }

      if (m_doFaceDeconvolution)
        { // otherwise, a_FfaceCen is UNUSED.
          Box faceBox0(m_currentBox);
          faceBox0.surroundingNodes(idir);

          const FArrayBox& WfaceCenDir = a_WfaceCen[idir];
          FArrayBox& FfaceCenDir = a_FfaceCen[idir];
          if (allFDir)
            {
              for (int jdir = 0; jdir < SpaceDim; jdir++)
                {
                  // Interval
                  int fluxIntLo = jdir*m_numFluxes;
                  int fluxIntHi = fluxIntLo + m_numFluxes-1;
                  Interval fluxIntv(fluxIntLo, fluxIntHi);
                  // Alias to pretend the fluxbox has only m_fluxes components
                  FArrayBox FfaceCenDirD(fluxIntv, FfaceCenDir);
                  m_molMappedPhysics->getFlux(FfaceCenDirD,
                                              WfaceCenDir,
                                              jdir,
                                              faceBox0);
                }
            }
          else
            {
              // Set FfaceCenDir = flux(WfaceCenDir)
              // on all idir-faces of m_currentBox
              // using WfaceCenDir
              // on all idir-faces of m_currentBox.
              m_molMappedPhysics->getFlux(FfaceCenDir, WfaceCenDir, idir, faceBox0);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::getFaceAvg(// we'll apply limiter to a_faceW
                              FluxBox& a_faceW,
                              const FArrayBox& a_cellW,
                              const FArrayBox& a_WofUavg,
                              FArrayBox& a_flatteningFab,
                              bool a_setFlattening)
{
  CH_TIME("PatchMappedConsOperator::getFaceAvg");
  CH_assert(isDefined());
  CH_assert(m_isCurrentBoxSet);

  int nComp = a_faceW.nComp();

  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_blockDomain;

  // Flattening stencil has size 3.
  // Compute flattening once for all slopes if needed.
  // Each CELL has a flattening coefficient, eta.
  // Flattening means replacing Wface on left and right of each cell
  // by eta*Wface + (1-eta)*Wcell.
  // = Wcell + eta * (Wface - Wcell).
  // That is:
  // WfaceL[+] = Wcell + eta * (WfaceL[+] - Wcell);
  // WfaceR[-] = Wcell + eta * (WfaceR[-] - Wcell).
  if (m_useFlattening && a_setFlattening)
    {
      Interval velInt    = m_molMappedPhysics->velocityInterval();
      int pressureInd    = m_molMappedPhysics->pressureIndex();
      Real smallPressure = m_molMappedPhysics->smallPressure();
      int bulkInd        = m_molMappedPhysics->bulkModulusIndex();
      // Set a_flatteningFab
      // on all cells of bx1inDomain = grow(m_currentBox, 1) & m_blockDomain
      // using a_WofUavg
      // on all cells of bx4inDomain = grow(m_currentBox, 4) & m_blockDomain.
      const Box& flatteningBox = a_flatteningFab.box();
      CH_assert(flatteningBox == bx1inDomain);
      m_util.computeFlattening(a_flatteningFab, a_WofUavg,
                               velInt, pressureInd, smallPressure, bulkInd,
                               flatteningBox);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // thisFaceWDir lives on idir-FACEs.
      FArrayBox& thisFaceWDir = a_faceW[idir];
      const Box& faceBox = thisFaceWDir.box();
      if (!m_noPPM) // if (m_noPPM) then set thisFaceWDir later.
        {
          /*
            Set thisFaceWDir from a_cellW by 4th-order interpolation
            on all idir-faces of bx1inDomain = grow(m_currentBox, 1) & m_blockDomain
            using a_cellW on bx4inDomain = grow(m_currentBox, 4) & m_blockDomain
            with the formula
            thisFaceWDir[i+e/2] = 7/12 * (a_cellW[i] + a_cellW[i+e])
                                - 1/12 * (a_cellW[i-e] + a_cellW[i+2e]).
          */
          m_util.PPMFaceValues(thisFaceWDir, a_cellW, nComp, m_limitFaceValues,
                               idir, faceBox, m_currentTime, m_molMappedPhysics);
        }

      /*
        On domain boundaries, set thisFaceWDir to a_cellW.
       */
      if (!m_domain.isPeriodic(idir))
        {
          PhysIBC* bc = m_molMappedPhysics->getPhysIBC();
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();
              int isign = sign(side);
              Box boundaryBox;
              bc->getBoundaryFaces(boundaryBox, thisFaceWDir.box(), idir, side);
              boundaryBox.shiftHalf(idir, -isign);
              thisFaceWDir.shiftHalf(idir, -isign);
              thisFaceWDir.copy(a_cellW, boundaryBox);
              thisFaceWDir.shiftHalf(idir, isign);
            }
        }

      if (m_limitFaceValues || m_useFlattening || m_noPPM)
        {
          FArrayBox WMinus(bx1inDomain, nComp);
          FArrayBox WPlus( bx1inDomain, nComp);

          // Set WMinus = thisFaceWDir[-idir/2] - a_cellW
          // on all cells of bx1inDomain = grow(m_currentBox, 1) & m_domain
          // using thisFaceWDir
          // on idir-faces of bx1inDomain except highest face
          // and a_cellW
          // on all cells of bx1inDomain.

          // Set WPlus = thisFaceWDir[+idir/2] - a_cellW
          // on all cells of bx1inDomain = grow(m_currentBox, 1) & m_domain
          // using thisFaceWDir
          // on idir-faces of bx1inDomain except lowest face
          // and a_cellW
          // on all cells of bx1inDomain.

          WMinus.setVal(0.);
          WPlus .setVal(0.);
          if (!m_noPPM)
            {
              // WMinus -= a_cellW;
              // WPlus  -= a_cellW;
              // Call minus() instead of the simpler -= because
              // the domains do not match.
              WMinus.minus(a_cellW, bx1inDomain, 0, 0, nComp);
              WPlus.minus( a_cellW, bx1inDomain, 0, 0, nComp);

              // WMinus[i] = thisFaceWDir[i - e/2] - a_cellW[i]
              // WPlus[i]  = thisFaceWDir[i + e/2] - a_cellW[i]
              thisFaceWDir.shiftHalf(idir, +1); // [i-e/2] becomes [i]
              WMinus.plus(thisFaceWDir, bx1inDomain, 0, 0, nComp);
              thisFaceWDir.shift(idir, -1); // [i+e/2] becomes [i]
              WPlus.plus( thisFaceWDir, bx1inDomain, 0, 0, nComp);
              thisFaceWDir.shiftHalf(idir, +1); // final shift back

              if (m_limitFaceValues)
                {
                  // Note:  petermc, 21 Sep 2010, changed from 2 layers to 3,
                  // with the new check-cubic limiter.

                  // Apply limiter to WMinus and WPlus
                  // on all cells of bx1inDomain = grow(m_currentBox, 1) & m_domain
                  // using a_cellW
                  // on bx4inDomain = grow(m_currentBox, 4) & m_domain.
                  // Recall:
                  // WMinus[i] = thisFaceWDir[i - e/2] - a_cellW.
                  // WPlus[i]  = thisFaceWDir[i + e/2] - a_cellW.
                  m_util.PPMLimiter(WMinus, WPlus, a_cellW, nComp, idir, bx1inDomain);
                }
              if (m_useFlattening)
                {
                  // Multiply WMinus and WPlus by a_flatteningFab
                  // on all cells of bx1inDomain = grow(m_currentBox, 1) & m_domain.
                  m_util.applyFlattening(WMinus, a_flatteningFab, bx1inDomain);
                  m_util.applyFlattening(WPlus,  a_flatteningFab, bx1inDomain);
                }
            } // if (!m_noPPM)
          // WMinus[i] = thisFaceWDir[i-e/2] on right after limiting
          //          WMinus += a_cellW;
          // WPlus[i] = thisFaceWDir[i+e/2] on left after limiting
          //          WPlus  += a_cellW;
          WMinus.plus(a_cellW, bx1inDomain, 0, 0, nComp);
          WPlus.plus( a_cellW, bx1inDomain, 0, 0, nComp);

          // If m_noPPM then WMinus == a_cellW == WPlus,
          // and we fill in thisFaceWDir by solving Riemann problem.

          // Now solve Riemann problem.
          // Modify thisFaceWDir[i + e/2] based on
          // WPlus[i] = WleftExtrap[i + e/2] and
          // WMinus[i + e] = WrightExtrap[i + e/2].

          // Solve Riemann problem:
          // Set thisFaceWDir
          // on interior idir-faces of bx1inDomain = grow(m_currentBox, 1) & m_domain
          // using WPlus and WMinus
          // on cells of bx1inDomain.
          Box riemannBox(bx1inDomain);
          riemannBox.surroundingNodes(idir);
          riemannBox.grow(idir, -1);

          // The riemann() function also sets boundary conditions
          // on thisFaceWDir, even though the external boundary faces
          // are NOT in riemannBox.
          m_molMappedPhysics->primBC(thisFaceWDir, // on idir-FACEs
                                     WPlus,        // on CELLs to left of idir-FACEs
                                     WMinus,       // on CELLs to right of idir-FACEs
                                     a_cellW,      // used? on CELLs
                                     m_currentTime,
                                     idir);        // on idir-FACEs
          
          // The riemann() function also transforms velocities.
          preRiemann(WPlus, WMinus, idir, riemannBox);
          m_molMappedPhysics->riemann(thisFaceWDir, // on idir-FACEs
                                      WPlus,        // on CELLs to left of idir-FACEs
                                      WMinus,       // on CELLs to right of idir-FACEs
                                      a_cellW,      // used? on CELLs
                                      m_currentTime,
                                      idir,
                                      riemannBox);  // on idir-FACEs
          postRiemann(thisFaceWDir, idir, riemannBox);
        }
      else
        { // neither limiting nor flattening,
          // therefore call Riemann solver only on physical boundary.
          PhysIBC* bc = m_molMappedPhysics->getPhysIBC();
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();

              Box boundaryBox;
              // bc->getBoundaryFaces(boundaryBox, faceBox, idir, side);
              // changed by petermc, 16 Mar 2010
              bc->getBoundaryFaces(boundaryBox, thisFaceWDir.box(), idir, side);
              if (! boundaryBox.isEmpty() )
                {
                  FArrayBox WPlusOrMinus(boundaryBox, nComp); // on FACEs
                  WPlusOrMinus.copy(thisFaceWDir, boundaryBox);
                  // This will fill in thisFaceWDir on boundaryBox only:
                  // it recalculates boundaryBox from thisFaceWDir.box().
                  // Notice:
                  // thisFaceWDir and WPlusOrMinus are FACE-centered;
                  // a_cellW is CELL-centered.  Derived classes know this.
                  // They call Fortran subroutines that also know this.
                  bc->primBC(thisFaceWDir,
                             WPlusOrMinus, a_cellW, idir, side, m_currentTime);
                }
            }
        }
    } // end loop over directions
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::preRiemann(FArrayBox&  a_WLeft,
                              FArrayBox&  a_WRight,
                              int         a_dir,
                              const Box&  a_box)
{
  CH_TIME("PatchMappedConsOperator::preRiemann");
  Interval velInt = m_molMappedPhysics->velocityInterval();
  if (velInt.begin() <= velInt.end())
  {
    FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
    FArrayBox& shiftWRight = (FArrayBox&)a_WRight;
    shiftWLeft .shiftHalf(a_dir, 1);
    shiftWRight.shiftHalf(a_dir,-1);
    CH_assert(shiftWLeft .box().contains(a_box));
    CH_assert(shiftWRight.box().contains(a_box));

    FArrayBox velLeftFab(velInt, shiftWLeft);
    FArrayBox velRightFab(velInt, shiftWRight);
    forwardBasisTransform(velLeftFab, a_dir, a_box);
    forwardBasisTransform(velRightFab, a_dir, a_box);

    shiftWLeft .shiftHalf(a_dir,-1);
    shiftWRight.shiftHalf(a_dir, 1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::postRiemann(FArrayBox&  a_Wface,
                               int         a_dir,
                               const Box&  a_box)
{
  CH_TIME("PatchMappedConsOperator::postRiemann");
 // Transform back velocity components in a_Wface.
  Interval velInt = m_molMappedPhysics->velocityInterval();
  if (velInt.begin() <= velInt.end())
  {
    FArrayBox velNewFab(velInt, a_Wface);
    reverseBasisTransform(velNewFab, a_dir, a_box);
  }
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::unitNormals(
  const LayoutData<FluxBox> *const a_unitNormalLayPtr)
{
  m_unitNormalLayPtr = a_unitNormalLayPtr;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::setMaxSpeed(Real a_maxSpeed)
{
  m_maxSpeed = a_maxSpeed;
  m_isMaxSpeedSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::useArtificialDissipation(bool a_useArtificialDissipation)
{
  m_useArtificialDissipation = a_useArtificialDissipation;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::artificialDissipation(Real a_artificialDissipation)
{
  m_artificialDissipation = a_artificialDissipation;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::addArtificialViscosity(FArrayBox&       a_UnewFab,
                                          const FArrayBox& a_UoldFab,
                                          FluxBox&         a_flux,
                                          Real             a_weight)
{
  CH_TIME("PatchMappedConsOperator::addArtificialViscosity");
  CH_assert(isDefined());
  CH_assert(m_useArtificialViscosity);
  CH_assert(m_isCurrentBoxSet);
  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_blockDomain;
  a_flux.define(m_currentBox, m_numFluxes);
  FArrayBox divFlux(m_currentBox, a_UnewFab.nComp());
  divFlux.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& currentFluxDir = a_flux[idir];
      currentFluxDir.setVal(0.);
      // Increment a_flux[idir]
      // on all idir-faces of m_currentBox
      // using a_UoldFab
      // on all cells of bx1inDomain = grow(a_box, 1) & m_blockDomain.
      m_molMappedPhysics->artVisc(currentFluxDir,
                                  a_UoldFab,
                                  m_artificialViscosity,
                                  m_currentTime,
                                  idir,
                                  m_currentBox);
      // Increment divFlux on all cells of m_currentBox
      // using a_flux[idir] on all idir-faces of m_currentBox.
      FORT_DIVERGENCE(CHF_CONST_FRA(currentFluxDir),
                      CHF_FRA(divFlux),
                      CHF_BOX(m_currentBox),
                      CHF_CONST_REAL(m_dx),
                      CHF_INT(idir));
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
  a_UnewFab.plus(divFlux, -a_weight);
  // Now LevelConsOperator needs to increment m_fluxes and flux registers.
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::addMappedArtificialViscosity(FArrayBox&       a_JUnewFab,
                                                const FArrayBox& a_UoldFab,
                                                FluxBox&         a_flux,
                                                const FluxBox&   a_N,
                                                const FArrayBox& a_J,
                                                const Box&       a_box,
                                                Real             a_weight)
{
  CH_TIME("PatchMappedConsOperator::addMappedArtificialViscosity");
  CH_assert(isDefined());
  CH_assert(m_useArtificialViscosity);
  Box bx1inDomain = grow(a_box, 1);
  bx1inDomain &= m_blockDomain;
  FArrayBox divFlux(a_box, m_molMappedPhysics->numFluxes());
  divFlux.setVal(0.);

  m_molMappedPhysics->mappedArtVisc(a_flux,
                              a_UoldFab,
                              a_N,
                              a_J,
                              m_artificialViscosity,
                              m_currentTime,
                              a_box);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& fluxDir = a_flux[idir];
      // Increment divFlux on all cells of m_currentBox
      // using a_flux[idir] on all idir-faces of m_currentBox.
      FORT_DIVERGENCE(CHF_CONST_FRA(fluxDir),
                      CHF_FRA(divFlux),
                      CHF_BOX(a_box),
                      CHF_CONST_REAL(m_dx),
                      CHF_INT(idir));
    }
  a_JUnewFab.plus(divFlux, -a_weight);
  // Now LevelConsOperator needs to increment m_fluxes and flux registers.
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::forwardBasisTransform(FArrayBox& a_W,
                                         const int a_dir,
                                         const Box& a_box)
{
  CH_TIME("PatchMappedConsOperator::forwardBasisTransform");
  CH_assert(m_unitNormalLayPtr != NULL);
  const FArrayBox& unitNormalFab =
    (*m_unitNormalLayPtr)[m_currentDataIdx][a_dir];
  CH_assert(unitNormalFab.box().contains(a_box));
  FORT_FORWARDTRANSFORMF(CHF_FRA(a_W),
                         CHF_CONST_FRA(unitNormalFab),
                         CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::reverseBasisTransform(FArrayBox& a_W,
                                         const int a_dir,
                                         const Box& a_box)
{
  CH_TIME("PatchMappedConsOperator::reverseBasisTransform");
  CH_assert(m_unitNormalLayPtr != NULL);
  const FArrayBox& unitNormalFab =
    (*m_unitNormalLayPtr)[m_currentDataIdx][a_dir];
  CH_assert(unitNormalFab.box().contains(a_box));
  FORT_REVERSETRANSFORMF(CHF_FRA(a_W),
                         CHF_CONST_FRA(unitNormalFab),
                         CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////

void
PatchMappedConsOperator::fluxCorrection(FluxBox&         a_fluxAvg,
                                        const FArrayBox& a_UavgFab)
{
  CH_TIME("PatchMappedConsOperator::fluxCorrection");
  CH_assert(m_isMaxSpeedSet);
  Real factor = -m_maxSpeed * m_artificialDissipation / (SpaceDim * 32.);
  // Loop over the spatial dimensions
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox& fluxAvgDir = a_fluxAvg[dir];
      const Box& fluxBox = fluxAvgDir.box();

      Box box1cells = grow(fluxBox, BASISV(dir));
      box1cells.enclosedCells();

      Box farLoFaces, midLoFaces, nearLoFaces;
      Box farHiFaces, midHiFaces, nearHiFaces;
      Box centerFaces, midCenterFaces, innerCenterFaces, entireFaces;
      int hasLoFaces, hasHiFaces;
      loHiCenterFace6(farLoFaces, midLoFaces, nearLoFaces, hasLoFaces,
                      farHiFaces, midHiFaces, nearHiFaces, hasHiFaces,
                      centerFaces, midCenterFaces, innerCenterFaces,
                      entireFaces,
                      box1cells, m_blockDomain, dir);

      // Recall that shallow-water physics has m_numFluxesPerField > 1.
      FArrayBox fifthDerivOfUdir(fluxBox, m_numFluxes);
      // Compute \frac{\partial^5 U}{\partial x_{dir}^5} h^5 
      // on each dir-face
      FORT_FIFTHDERIVFACES(CHF_FRA(fifthDerivOfUdir),
                           CHF_CONST_FRA(a_UavgFab),
                           CHF_CONST_INT(m_numFluxes),
                           CHF_CONST_INT(dir),
                           CHF_BOX(farLoFaces),
                           CHF_BOX(midLoFaces),
                           CHF_BOX(nearLoFaces),
                           CHF_CONST_INT(hasLoFaces),
                           CHF_BOX(farHiFaces),
                           CHF_BOX(midHiFaces),
                           CHF_BOX(nearHiFaces),
                           CHF_CONST_INT(hasHiFaces),
                           CHF_BOX(innerCenterFaces));
      // Add fifth-derivative term to flux
      if (m_numFluxesPerField == 1)
        {
          fluxAvgDir.plus(fifthDerivOfUdir, factor,
                          0, 0, m_numFluxes);
        }
      else
        {
          // components of fluxAvgDir are grouped by dimension,
          // m_numFluxes at a time.
          CH_assert(m_numFluxesPerField == SpaceDim);
          int dstComp = 0;
          for (int dirFlux = 0; dirFlux < SpaceDim; dirFlux++)
            {
              fluxAvgDir.plus(fifthDerivOfUdir, factor,
                              0, dstComp, m_numFluxes);
              dstComp += m_numFluxes;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

Real
PatchMappedConsOperator::getMaxWaveSpeed(const FArrayBox&   a_U,
                                         const Box&         a_bx)
{
  CH_TIME("PatchMappedConsOperator::getMaxWaveSpeed");
  return m_molMappedPhysics->getMaxWaveSpeed(a_U, a_bx);
}

//////////////////////////////////////////////////////////////////////////////

#include "NamespaceFooter.H"
