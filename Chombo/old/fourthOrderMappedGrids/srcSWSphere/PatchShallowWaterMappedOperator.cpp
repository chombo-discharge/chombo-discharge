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

#include "PatchShallowWaterMappedOperator.H"
#include "MOLShallowWaterPhysics.H"
#include "SWintegrator.H"
#include "MatrixVectorTransformF_F.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
PatchShallowWaterMappedOperator::PatchShallowWaterMappedOperator()
  :
  PatchMappedConsOperator()
{
  m_contravariantMetricFacePtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
PatchShallowWaterMappedOperator::~PatchShallowWaterMappedOperator()
{
}


//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
void PatchShallowWaterMappedOperator::setCurrentBox(const Box& a_currentBox)
{
  PatchMappedConsOperator::setCurrentBox(a_currentBox);
  MOLShallowWaterPhysics* swPhysics = (MOLShallowWaterPhysics*) m_molPhysics;
  swPhysics->setCurrentCoordSys(m_coordSysCurrentPtr);
}


//////////////////////////////////////////////////////////////////////////////
void
PatchShallowWaterMappedOperator::getNormalFlux(
                                       FluxBox&            a_FfaceAvg,
                                       FluxBox&            a_WfaceOrthoCen,
                                       const FArrayBox&    a_UavgFab,
                                       Real                a_weight,
                                       bool                a_setFlattening,
                                       FArrayBox&          a_flatteningFab)
{
  CH_TIME("PatchShallowWaterMappedOperator::getNormalFlux");
  CH_assert(isDefined());
  CH_assert(m_isCurrentBoxSet);

  /*
    Get face-centered N.
  */
  Box bx1 = grow(m_currentBox, 1);
  int nFaceMetricTerms = SpaceDim*SpaceDim;
  // m_coordSysCurrentPtr and m_faceMetricTerms will be used in reduceFlux().
  m_faceMetricTerms = new FluxBox(bx1, nFaceMetricTerms);
  m_coordSysCurrentPtr->getN(*m_faceMetricTerms, bx1);
  //  PatchConsOperator::getNormalFlux(a_FfaceAvg, a_UavgFab, a_weight,
  //                                   a_setFlattening, a_flatteningFab);
  // Copy this whole function only because you need a_WfaceCen!
  //  CH_assert(isDefined());
  //  CH_assert(m_isCurrentBoxSet);
  int numU = a_UavgFab.nComp();
  CH_assert(m_numGhost >= 5);
  Box bx5inDomain = grow(m_currentBox, 5);
  bx5inDomain &= m_domain;
  Box bx4inDomain = grow(m_currentBox, 4);
  bx4inDomain &= m_domain;
  FArrayBox UcellCenFab(bx4inDomain, numU);
  UcellCenFab.copy(a_UavgFab);
  m_util.deconvolve(UcellCenFab, a_UavgFab, bx4inDomain, -1);
  int numW = m_molPhysics->numConserved();
  FArrayBox WcellCenFab(bx4inDomain, numW);
  m_molPhysics->consToPrim(WcellCenFab, UcellCenFab, bx4inDomain);
  FArrayBox WofUavgFab(bx5inDomain, numW);
  m_molPhysics->consToPrim(WofUavgFab, a_UavgFab, bx5inDomain);
  FArrayBox WcellAvgFab(bx4inDomain, numW);
  WcellAvgFab.copy(WcellCenFab);
  m_util.deconvolve(WcellAvgFab, WofUavgFab, bx4inDomain);
  // Box bx1 = grow(m_currentBox, 1);
  Box bx1inDomain(bx1);
  bx1inDomain &= m_domain;
  Box bx2inDomain = grow(m_currentBox, 2);
  bx2inDomain &= m_domain;
  FluxBox WfaceAvg(bx2inDomain, numW); // WAS on bx1inDomain
  getFaceAvg(WfaceAvg, WcellAvgFab, WofUavgFab,
             a_flatteningFab, a_setFlattening);
  FluxBox WfaceCen(bx1inDomain, numW); // WAS on m_currentBox
  WfaceCen.copy(WfaceAvg);
  m_util.deconvolveFace(WfaceCen, WfaceAvg, bx1inDomain, -1);
  FluxBox FfaceCen(bx1, m_numFluxesPerField * m_numFluxes); // WAS on m_currentBox
  FluxBox FfromWfaceAvg(bx1, m_numFluxesPerField * m_numFluxes);
  getAllFluxes(FfromWfaceAvg, FfaceCen, WfaceAvg, WfaceCen);
  a_FfaceAvg.resize(bx1inDomain, m_numFluxesPerField * m_numFluxes);
  a_FfaceAvg.copy(FfaceCen);
  m_util.deconvolveFace(a_FfaceAvg, FfromWfaceAvg, m_currentBox);
  reduceFlux(a_FfaceAvg, FfromWfaceAvg);

  // Need to return a_WfaceOrthoCen.
  a_WfaceOrthoCen.resize(bx1, numW); // WAS m_currentBox
  Interval velInt = m_molPhysics->velocityInterval();
  Vector<int> velocityVector;
  for (int velComp = velInt.begin(); velComp <= velInt.end(); velComp++)
    {
      velocityVector.push_back(velComp);
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const FArrayBox& WfaceCenDir = WfaceCen[idir];
      const Box& faceBox = WfaceCenDir.box();

      const FArrayBox& orthoMatrix = (*m_orthoMatrixFacePtr)[idir];
      FArrayBox& WfaceOrthoCenDir = a_WfaceOrthoCen[idir];
      // csWcomps: component indices of velocity in coordinate-system basis
      // IntVect csWcomps = IntVect(D_DECL(WVELX, WVELY, WVELZ));
      // orthoWcomps: component indices of velocity in orthonormal basis
      // IntVect orthoWcomps = csWcomps;
      // a_WfaceOrthoCen[idir][WHGT] := WfaceCen[idir][WHGT] == h
      WfaceOrthoCenDir.copy(WfaceCenDir, WHGT, WHGT);
      // a_WfaceOrthoCen[idir][WVELX:WVELY] :=
      // orthonormalize(WfaceCen[idir][WVELX:WVELY]) == velocity
      // panelCS->orthonormalize(WfaceCen[idir], a_WfaceOrthoCen[idir], faceBox,
      //                         idir, csWcomps, orthoWcomps);
      FORT_MATVECMULTONBOX(CHF_CONST_FRA(orthoMatrix),
                           CHF_CONST_FRA(WfaceCenDir),
                           CHF_FRA(WfaceOrthoCenDir),
                           CHF_BOX(faceBox),
                           CHF_CONST_VI(velocityVector),
                           CHF_CONST_VI(velocityVector));
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
  delete m_faceMetricTerms;
}


//////////////////////////////////////////////////////////////////////////////
void
PatchShallowWaterMappedOperator::getAllFluxes(FluxBox&        a_FfaceAvg,
                                              FluxBox&        a_FfaceCen,
                                              const FluxBox&  a_WfaceAvg,
                                              const FluxBox&  a_WfaceCen)
{
  CH_TIME("PatchShallowWaterMappedOperator::getAllFluxes");
  CH_assert(m_isCurrentBoxSet);
  MOLShallowWaterPhysics* swPhysics = (MOLShallowWaterPhysics*) m_molPhysics;
  swPhysics->setContravariantMetricFace(m_contravariantMetricFacePtr);

  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_domain;
  for (int idir = 0; idir < SpaceDim; idir++)
    { // idir specifies which faces:  either x or y or z
      Box faceBox1(bx1inDomain);
      faceBox1.surroundingNodes(idir);
      getDirFluxes(a_FfaceAvg[idir], a_WfaceAvg[idir], faceBox1);

      // Box faceBox0(m_currentBox);
      // faceBox0.surroundingNodes(idir);
      getDirFluxes(a_FfaceCen[idir], a_WfaceCen[idir], faceBox1);
    }
}

//////////////////////////////////////////////////////////////////////////////
void
PatchShallowWaterMappedOperator::preRiemann(FArrayBox&  a_WLeft,
                                      FArrayBox&  a_WRight,
                                      int         a_dir,
                                      const Box&  a_box)
{
/*
  // No need to transform the advection velocity.
  FArrayBox& advVelFab = (*m_advVelFacePtr)[a_dir];
  CH_assert(advVelFab.box().contains(a_box));
  MOLAdvectionPhysics* advectionPhysics =
    (MOLAdvectionPhysics*) m_molPhysics;
  advectionPhysics->setVelocityFab(&advVelFab);
*/
}

//////////////////////////////////////////////////////////////////////////////
void
PatchShallowWaterMappedOperator::getFaceAvg(// we'll apply limiter to a_faceW
                                            FluxBox& a_faceW,
                                            const FArrayBox& a_cellW,
                                            const FArrayBox& a_WofUavg,
                                            FArrayBox& a_flatteningFab,
                                            bool a_setFlattening)
{
  CH_TIME("PatchShallowWaterMappedOperator::getFaceAvg");
  CH_assert(isDefined());
  CH_assert(m_isCurrentBoxSet);

  int nComp = a_faceW.nComp();

  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_domain;

  Box bx2inDomain = grow(m_currentBox, 2);
  bx2inDomain &= m_domain;

  // No slope flattening for shallow-water equations.

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // thisFaceWDir lives on idir-FACEs.
      FArrayBox& thisFaceWDir = a_faceW[idir];
      // faceBox for shallow-water equations
      // const Box& faceBox = thisFaceWDir.box();
      Box faceBox(thisFaceWDir.box());
      faceBox.grow(idir, -1);
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
              bc->getBoundaryFaces(boundaryBox, faceBox, idir, side);
              boundaryBox.shiftHalf(idir, -isign);
              thisFaceWDir.shiftHalf(idir, -isign);
              thisFaceWDir.copy(a_cellW, boundaryBox);
              thisFaceWDir.shiftHalf(idir, isign);
            }
        }

      if (m_limitFaceValues || m_noPPM)
        {
          FArrayBox WMinus(bx2inDomain, nComp);
          FArrayBox WPlus( bx2inDomain, nComp);

          WMinus.setVal(0.);
          WPlus .setVal(0.);
          if (!m_noPPM)
            {
              WMinus.minus(a_cellW, bx2inDomain, 0, 0, nComp);
              WPlus.minus( a_cellW, bx2inDomain, 0, 0, nComp);

              thisFaceWDir.shiftHalf(idir, +1); // [i-e/2] becomes [i]
              WMinus.plus(thisFaceWDir, bx2inDomain, 0, 0, nComp);
              thisFaceWDir.shift(idir, -1); // [i+e/2] becomes [i]
              WPlus.plus( thisFaceWDir, bx2inDomain, 0, 0, nComp);
              thisFaceWDir.shiftHalf(idir, +1); // final shift back

              if (m_limitFaceValues)
                {
                  m_util.PPMLimiter(WMinus, WPlus, a_cellW, nComp, idir, bx2inDomain);
                }
            } // if (!m_noPPM)
          WMinus.plus(a_cellW, bx2inDomain, 0, 0, nComp);
          WPlus.plus( a_cellW, bx2inDomain, 0, 0, nComp);

          Box riemannBox(bx2inDomain);
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
              bc->getBoundaryFaces(boundaryBox, faceBox, idir, side);
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
PatchShallowWaterMappedOperator::postRiemann(FArrayBox&  a_Wface,
                                       int         a_dir,
                                       const Box&  a_box)
{
}

//////////////////////////////////////////////////////////////////////////////
void
PatchShallowWaterMappedOperator::setContravariantMetricFace(const FluxBox* a_contravariantMetricFacePtr)
{
  m_contravariantMetricFacePtr = (FluxBox*) a_contravariantMetricFacePtr;
}

//////////////////////////////////////////////////////////////////////////////
void
PatchShallowWaterMappedOperator::setOrthoMatrix(const FluxBox* a_orthoMatrixFacePtr)
{
  m_orthoMatrixFacePtr = (FluxBox*) a_orthoMatrixFacePtr;
}

#include "NamespaceFooter.H"
