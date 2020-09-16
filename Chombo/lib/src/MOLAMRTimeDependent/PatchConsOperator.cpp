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

#include "PatchConsOperator.H"
#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
PatchConsOperator::PatchConsOperator()
{
  m_molPhysics = NULL;
  m_isDefined = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet  = false;
  m_dx = 0.0;
  m_numFluxesPerField = 1;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
PatchConsOperator::~PatchConsOperator()
{
  if (m_molPhysics != NULL)
    {
      delete m_molPhysics;
    }
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void PatchConsOperator::define(
                               const ProblemDomain&      a_domain,
                               const Real&               a_dx,
                               const MOLPhysics* const a_molPhysics,
                               const int&                a_numFields)
{
  CH_TIME("PatchConsOperator::define");
  // Sanity checks
  CH_assert(a_dx > 0.0);

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;

  m_numFields = a_numFields;

  // MOLPhysics* m_molPhysics;
  m_molPhysics = a_molPhysics->new_molPhysics();
  m_molPhysics->define(m_domain, m_dx);

  m_numFluxes = m_molPhysics->numFluxes();

  // MOLUtilities m_util;
  m_util.define(m_domain, m_dx);
  m_highOrderLimiter = true;

  m_doDeconvolution = true;
  m_noPPM = false;
  m_useArtificialViscosity = false;

  // Everything is defined now.
  m_isDefined = true;
}

//////////////////////////////////////////////////////////////////////////////
void PatchConsOperator::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime      = a_currentTime;
  m_isCurrentTimeSet = true;
}

//////////////////////////////////////////////////////////////////////////////
void PatchConsOperator::setCurrentBox(const Box& a_currentBox)
{
  m_currentBox      = a_currentBox;
  m_isCurrentBoxSet = true;

  m_molPhysics->setCurrentBox(m_currentBox);
}


//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::getNormalFlux(
                                 FluxBox&            a_FfaceAvg,
                                 const FArrayBox&    a_UavgFab,
                                 Real                a_weight,
                                 bool                a_setFlattening,
                                 FArrayBox&          a_flatteningFab)
{
  CH_TIME("PatchConsOperator::getNormalFlux");
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
  bx5inDomain &= m_domain;
  Box bx4inDomain = grow(m_currentBox, 4);
  bx4inDomain &= m_domain;
  FArrayBox UcellCenFab(bx4inDomain, numU);
  UcellCenFab.copy(a_UavgFab);
  if (m_doDeconvolution)
    { // UcellCenFab -= Laplacian(a_UavgFab) * m_dx^2 / 24
      // fourthOrderAverageCell(UcellCenFab, -1);
      m_util.deconvolve(UcellCenFab, a_UavgFab, bx4inDomain, -1);
    }

  /*
    Set WcellCenFab = W(UcellCenFab)
    on all cells of bx4inDomain = grow(m_currentBox, 4) & m_domain
    using UcellCenFab
    on all cells of bx4inDomain.
  */
  int numW = m_molPhysics->numConserved();
  FArrayBox WcellCenFab(bx4inDomain, numW);
  m_molPhysics->consToPrim(WcellCenFab, UcellCenFab, bx4inDomain);

  /*
    Set WofUavgFab = W(a_UavgFab)
    on all cells of bx5inDomain = grow(m_currentBox, 5) & m_domain
    using a_UavgFab
    on all cells of bx5inDomain.
  */
  FArrayBox WofUavgFab(bx5inDomain, numW);
  m_molPhysics->consToPrim(WofUavgFab, a_UavgFab, bx5inDomain);

  /*
    Set WcellAvgFab = WofUavgFab + m_dx^2/24 * Laplacian(WofUavgFab)
    on all cells of bx4inDomain = grow(m_currentBox, 4) & m_domain
    using WofUavgFab
    on all cells of bx5inDomain = grow(m_currentBox, 5) & m_domain.
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
  bx1inDomain &= m_domain;

  FluxBox WfaceAvg(bx1inDomain, numW);
  /*
    Set WfaceAvg
    on all faces of bx1inDomain = grow(m_currentBox, 1) & m_domain,
    using WcellAvgFab
    on bx4inDomain = grow(m_currentBox, 4) & m_domain,
    and WofUavgFab
    on bx4inDomain = grow(m_currentBox, 4) & m_domain,
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
    on all faces of bx1inDomain = grow(m_currentBox, 1) & m_domain.
  */
  FluxBox WfaceCen(m_currentBox, numW);
  // First set WfaceCen = WfaceAvg.
  WfaceCen.copy(WfaceAvg);
  // Then set WfaceCen -= Laplacian(WfaceAvg) * m_dx^2 / 24
  // on all faces of m_currentBox
  // using WfaceAvg
  // on all faces of m_currentBox.
  // Special stencils next to boundary of m_domain.
  m_util.deconvolveFace(WfaceCen, WfaceAvg, m_currentBox, -1);

  /*
    Set a_FfromWfaceAvg = F(WfaceAvg)
    on all faces of bx1inDomain = grow(m_currentBox, 1) & m_domain
    using WfaceAvg
    on all faces of bx1inDomain.

    Set FfaceCen = F(WfaceCen)
    on all faces of m_currentBox
    using WfaceCen
    on all faces of m_currentBox.
  */
  FluxBox FfaceCen(m_currentBox, m_numFluxesPerField * m_numFluxes);
  FluxBox FfromWfaceAvg(bx1, m_numFluxesPerField * m_numFluxes);
  getAllFluxes(FfromWfaceAvg, FfaceCen, WfaceAvg, WfaceCen);
  a_FfaceAvg.define(bx1inDomain, m_numFluxesPerField * m_numFluxes);
  if (m_doFaceDeconvolution)
    {
      /*
        Set a_FfaceAvg = FfaceCen + m_dx^2/24 * Laplacian(FfromWfaceAvg)
        on all idir-faces of m_currentBox
        using FfaceCen
        on all idir-faces of m_currentBox
        and original a_FfaceAvg
        on all idir-faces of bx1inDomain = grow(m_currentBox, 1) & m_domain.
      */
      // First set a_FfaceAvg = FfaceCen = flux(WfaceCen)
      // on all faces of m_currentBox.
      a_FfaceAvg.copy(FfaceCen);
      // Then set a_FfaceAvg += Laplacian(FfromWfaceAvg) * m_dx^2 / 24
      // on all faces of m_currentBox,
      // using FfromWfaceAvg
      // on all faces of bx1inDomain = grow(m_currentBox, 1) & m_domain.
      // Special stencils next to boundary of m_domain.
      m_util.deconvolveFace(a_FfaceAvg, FfromWfaceAvg, m_currentBox);
    }
  else
    {
      a_FfaceAvg.copy(FfromWfaceAvg);
    }

  reduceFlux(a_FfaceAvg, FfromWfaceAvg);
}


//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::reduceFlux(FluxBox&     a_FfaceAvg,
                              FluxBox&     a_FfaceForGrad)
{
  // This function intentionally left blank.
  // Derived classes may use it.
}

//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::getFluxDivergence(FArrayBox&   a_LofU,
                                     const FluxBox&     a_FfaceAvg)
{
  CH_TIME("PatchConsOperator::getFluxDivergence");
  a_LofU.setVal(0.0);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const FArrayBox& currentFluxDir = a_FfaceAvg[idir];
      // Increment a_LofU on all cells of m_currentBox
      // using a_FfaceAvg[idir] on all idir-faces of m_currentBox.
      FORT_DIVERGENCE(CHF_CONST_FRA(currentFluxDir),
                      CHF_FRA(a_LofU),
                      CHF_BOX(m_currentBox),
                      CHF_CONST_REAL(m_dx),
                      CHF_INT(idir));
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
}


//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::getAllFluxes(FluxBox&        a_FfaceAvg,
                                FluxBox&        a_FfaceCen,
                                const FluxBox&  a_WfaceAvg,
                                const FluxBox&  a_WfaceCen)
{
  CH_TIME("PatchConsOperator::getAllFluxes");
  CH_assert(m_isCurrentBoxSet);
  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_domain;
  for (int idir = 0; idir < SpaceDim; idir++)
    { // idir specifies which faces:  either x or y or z
      Box faceBox1(bx1inDomain);
      faceBox1.surroundingNodes(idir);

      const FArrayBox& WfaceAvgDir = a_WfaceAvg[idir];
      FArrayBox& FfaceAvgDir = a_FfaceAvg[idir];
      // Set FfaceAvgDir = flux(WfaceAvgDir)
      // on all idir-faces of bx1inDomain = grow(m_currentBox, 1) & m_domain.
      m_molPhysics->getFlux(FfaceAvgDir, WfaceAvgDir, idir, faceBox1);

      if (m_doFaceDeconvolution)
        { // otherwise, a_FfaceCen is UNUSED.
          Box faceBox0(m_currentBox);
          faceBox0.surroundingNodes(idir);

          const FArrayBox& WfaceCenDir = a_WfaceCen[idir];
          FArrayBox& FfaceCenDir = a_FfaceCen[idir];
          // Set FfaceCenDir = flux(WfaceCenDir)
          // on all idir-faces of m_currentBox
          // using WfaceCenDir
          // on all idir-faces of m_currentBox.
          m_molPhysics->getFlux(FfaceCenDir, WfaceCenDir, idir, faceBox0);
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::updateODE(FArrayBox&        a_solnFab,
                             const FArrayBox&  a_rhsFab,
                             Real a_dt)
{
  CH_TIME("PatchConsOperator::updateODE");
  int rhsComp = a_rhsFab.nComp();
  CH_assert(a_solnFab.nComp() == rhsComp);

  const Box& bx = a_rhsFab.box();
  a_solnFab.plus(a_rhsFab, bx, bx, a_dt, 0, 0, rhsComp);
}


//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::getFaceAvg(// we'll apply limiter to a_faceW
                              FluxBox& a_faceW,
                              const FArrayBox& a_cellW,
                              const FArrayBox& a_WofUavg,
                              FArrayBox& a_flatteningFab,
                              bool a_setFlattening)
{
  CH_TIME("PatchConsOperator::getFaceAvg");
  CH_assert(isDefined());
  CH_assert(m_isCurrentBoxSet);

  int nComp = a_faceW.nComp();

  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_domain;

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
      Interval velInt    = m_molPhysics->velocityInterval();
      int pressureInd    = m_molPhysics->pressureIndex();
      Real smallPressure = m_molPhysics->smallPressure();
      int bulkInd        = m_molPhysics->bulkModulusIndex();
      // Set a_flatteningFab
      // on all cells of bx1inDomain = grow(m_currentBox, 1) & m_domain
      // using a_WofUavg
      // on all cells of bx4inDomain = grow(m_currentBox, 4) & m_domain.
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
            on all idir-faces of bx1inDomain = grow(m_currentBox, 1) & m_domain
            using a_cellW on bx4inDomain = grow(m_currentBox, 4) & m_domain
            with the formula
            thisFaceWDir[i+e/2] = 7/12 * (a_cellW[i] + a_cellW[i+e])
                                - 1/12 * (a_cellW[i-e] + a_cellW[i+2e]).
          */
          m_util.PPMFaceValues(thisFaceWDir, a_cellW, nComp, m_limitFaceValues,
                               idir, faceBox, m_currentTime, m_molPhysics);
        }

      /*
        On domain boundaries, set thisFaceWDir to a_cellW.
       */
      if (!m_domain.isPeriodic(idir))
        {
          PhysIBC* bc = m_molPhysics->getPhysIBC();
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
          // The riemann() function also transforms velocities.
          preRiemann(WPlus, WMinus, idir, riemannBox);
          m_molPhysics->riemann(thisFaceWDir, // on idir-FACEs
                                WPlus, // on CELLs to left of idir-FACEs
                                WMinus, // on CELLs to right of idir-FACEs
                                a_cellW, // used? on CELLs
                                m_currentTime,
                                idir,
                                riemannBox); // on idir-FACEs
          postRiemann(thisFaceWDir, idir, riemannBox);
          // dummy statement in order to get around gdb bug
          int dummy_unused = 0; dummy_unused = 0;
        }
      else
        { // neither limiting nor flattening,
          // therefore call Riemann solver only on physical boundary.
          PhysIBC* bc = m_molPhysics->getPhysIBC();
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
PatchConsOperator::preRiemann(FArrayBox&  a_WLeft,
                              FArrayBox&  a_WRight,
                              int         a_dir,
                              const Box&  a_box)
{
}


//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::postRiemann(FArrayBox&  a_Wface,
                               int         a_dir,
                               const Box&  a_box)
{
}


//////////////////////////////////////////////////////////////////////////////
void
PatchConsOperator::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}


void
PatchConsOperator::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

void
PatchConsOperator::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
  m_util.highOrderLimiter(m_highOrderLimiter);
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::numGhost(int a_numGhost)
{
  m_numGhost = a_numGhost;
}

//////////////////////////////////////////////////////////////////////////////

void
PatchConsOperator::addArtificialViscosity(FArrayBox&         a_UnewFab,
                                          const FArrayBox&   a_UoldFab,
                                          FluxBox&           a_flux,
                                          Real  a_weight)
{
  CH_TIME("PatchConsOperator::addArtificialViscosity");
  CH_assert(isDefined());
  CH_assert(m_useArtificialViscosity);
  CH_assert(m_isCurrentBoxSet);
  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_domain;
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
      // on all cells of bx1inDomain = grow(a_box, 1) & m_domain.
      m_molPhysics->artVisc(currentFluxDir,
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
bool PatchConsOperator::isDefined() const
{
  return m_isDefined        &&
         m_isCurrentTimeSet &&
         m_isCurrentBoxSet  ;
}

#include "NamespaceFooter.H"
