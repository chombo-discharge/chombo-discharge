#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
#include <sstream>
using std::cout;
using std::endl;
using std::ofstream;
#include <limits> // contains numeric_limits<Real>

#include "parstream.H"
#include "CubedSphereShellCS.H"
#include "CubedSphereShellPanelCS.H"
#include "ClimateIBCF_F.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "LoadBalance.H"
#include "BlockRegister.H"
#include "MBMiscUtil.H"
#include "computeNorm.H"
#include "DebugOut.H"

#include "FABView.H"
#include "DebugDump.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testCubedSphereBlockRegister";
static const char* indent = "   ";
static const char* indent2 = "      ";
static bool verbose = false ;
// static bool verbose = true ;

static bool showSame = false;
static bool showFlux = false;
static bool showAux = false;
static bool showScalar = false;
static bool showVector = false;

static char letts[] = {'x', 'y', 'z'};

static Real precision = 2 * numeric_limits<Real>::epsilon();

//-----------------------------------------------------------------------
string
nameFace(int a_dir, Side::LoHiSide a_side )
{
  stringstream ss;
  ss << letts[a_dir];
  if (a_side == Side::Lo)
    ss << "Lo";
  else if (a_side == Side::Hi)
    ss << "Hi";
  else
    {
      MayDay::Error("Calling nameFace with undefined side");
      return "";
    }
  return ss.str();
}

//-----------------------------------------------------------------------
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if ( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------

void evalScalarFunction(FArrayBox& a_funFab,
                        const FArrayBox& a_xiFab,
                        const NewCoordSys* a_coordSysPtr)
{
  const Box& bx = a_xiFab.box();
  FArrayBox rllFab(bx, SpaceDim);
  const CubedSphereShellPanelCS* cssPanelPtr =
    dynamic_cast<const CubedSphereShellPanelCS*>(a_coordSysPtr);
  cssPanelPtr->fabTransformEquiangularToLonLat(a_xiFab, rllFab);
  Real magnitude = 2.0 * M_PI;
  Real shift = 0.;
  Real angle = M_PI / 3.0;
  FArrayBox vecRLLFab(bx, SpaceDim);
  FORT_CSPVECLONLATSOLIDBODYFLOW(CHF_FRA(vecRLLFab),
                                 CHF_CONST_FRA(rllFab),
                                 CHF_CONST_REAL(magnitude),
                                 CHF_CONST_REAL(angle),
                                 CHF_CONST_REAL(shift));
#if CH_SPACEDIM >= 3
  // vecRLLFab.setVal(1., 2);
  // Set vecRLLFab[2] = vecRLLFab[0]^2 + vecRLLFab[1]^2.
  FArrayBox velLon2Fab(bx, 1);
  velLon2Fab.copy(vecRLLFab, 0, 0);
  velLon2Fab.mult(vecRLLFab, 0, 0);
  FArrayBox velLat2Fab(bx, 1);
  velLat2Fab.copy(vecRLLFab, 1, 0);
  velLat2Fab.mult(vecRLLFab, 1, 0);
  vecRLLFab.copy(velLon2Fab, 0, 2);
  vecRLLFab.plus(velLat2Fab, 0, 2);
#endif
  FArrayBox vecEquFab(bx, SpaceDim);
  cssPanelPtr->fabVectorTransformLatLonToEquiangular(a_xiFab, vecRLLFab, vecEquFab);
  int dirNormal = faceDimension(bx);
  // Return normal component of velocity.
  a_funFab.copy(vecEquFab, dirNormal, 0);
}

//-----------------------------------------------------------------------

void evalVectorFunction(FArrayBox& a_funFab,
                        const FArrayBox& a_xiFab,
                        const NewCoordSys* a_coordSysPtr)
{
  const Box& bx = a_xiFab.box();
  FArrayBox rllFab(bx, SpaceDim);
  const CubedSphereShellPanelCS* cssPanelPtr =
    dynamic_cast<const CubedSphereShellPanelCS*>(a_coordSysPtr);
  cssPanelPtr->fabTransformEquiangularToLonLat(a_xiFab, rllFab);
  Real magnitude = 2.0 * M_PI;
  Real shift = 0.;
  Real angle = M_PI / 3.0;
  FArrayBox vecRLLFab(bx, SpaceDim);
  FORT_CSPVECLONLATSOLIDBODYFLOW(CHF_FRA(vecRLLFab),
                                 CHF_CONST_FRA(rllFab),
                                 CHF_CONST_REAL(magnitude),
                                 CHF_CONST_REAL(angle),
                                 CHF_CONST_REAL(shift));
#if CH_SPACEDIM >= 3
  // vecRLLFab.setVal(1., 2);
  // Set vecRLLFab[2] = vecRLLFab[0]^2 + vecRLLFab[1]^2.
  FArrayBox velLon2Fab(bx, 1);
  velLon2Fab.copy(vecRLLFab, 0, 0);
  velLon2Fab.mult(vecRLLFab, 0, 0);
  FArrayBox velLat2Fab(bx, 1);
  velLat2Fab.copy(vecRLLFab, 1, 0);
  velLat2Fab.mult(vecRLLFab, 1, 0);
  vecRLLFab.copy(velLon2Fab, 0, 2);
  vecRLLFab.plus(velLat2Fab, 0, 2);
#endif
  FArrayBox vecEquFab(bx, SpaceDim);
  cssPanelPtr->fabVectorTransformLatLonToEquiangular(a_xiFab, vecRLLFab, vecEquFab);
  int dirNormal = faceDimension(bx);
  // Return velocity vector multiplied by normal component of velocity.
  a_funFab.copy(vecEquFab);
  for (int comp = 0; comp < SpaceDim; comp++)
    {
      a_funFab.mult(vecEquFab, dirNormal, comp);
    }
}

//-----------------------------------------------------------------------

int
cubedSphereTest()
{
  // Create a cubed sphere coordinate system with six N-by-N blocks.
#if CH_SPACEDIM >= 3
  // Avoid irritating unused variable warning if dimension < 3.
  int nlayers = 16;
#endif
  int N = 32;

  //  int layer = 0; // which layer in 3D
  // CubedSphereShellCS* coordSys = new CubedSphereShellCS(N, nlayers);

  IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
  ProblemDomain pd(IntVect::Zero, domainLen-IntVect::Unit);

  CubedSphereShellCSFactory* cssFactory = new CubedSphereShellCSFactory();
  cssFactory->setHeight(1.);
  cssFactory->setRadius(1.);
  RealVect probLen = RealVect(D_DECL(11*(M_PI/2.), M_PI/2., 1.));
  RealVect dx = probLen / RealVect(domainLen);
  MultiBlockCoordSys* coordSys = cssFactory->getCoordSys(pd, dx);

  const Vector<Box>& blockBoxes =
    coordSys->mappingBlocks();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSys->boundaries();

  int nblocks = blockBoxes.size();

  int refRatio = 2;
  IntVect totalLen = IntVect(D_DECL(N, N, nlayers));
  IntVect partialLen = totalLen / refRatio;
  Box partialBaseBox(IntVect::Zero, partialLen-IntVect::Unit);

  Box refBox(IntVect::Zero, (refRatio-1) * IntVect::Unit);

  Vector<Box> boxes;
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      const Box& blockBox = blockBoxes[iblock];

      const IntVect& blockLo = blockBox.smallEnd();
      for (BoxIterator bitRef(refBox); bitRef.ok(); ++bitRef)
        {
          IntVect bu = bitRef();
          IntVect thisShift = partialLen * bu + blockLo;
          Box thisBox(partialBaseBox);
          thisBox.shift(thisShift);
          boxes.push_back(thisBox);
        }
    }

  mortonOrdering(boxes);
  Vector<int> procs(boxes.size());
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, pd);
  DataIterator dit = grids.dataIterator();

  Real worstRatio = 0.;
  bool bad;

  int status = 0;
  //  int verbosity = 2;
  
  Interval baseVectorIntv(0, SpaceDim-1);
  Interval vectorIntv(baseVectorIntv);
  int nscalar = 1;
  int iscalar = SpaceDim;
  Interval scalarIntv(iscalar, iscalar+nscalar-1);
  int ncomp = SpaceDim + nscalar;

  int numGhost = 1;
  IntVect ghostVect = numGhost * IntVect::Unit;

  LevelData<FluxBox> flux(grids, ncomp);
  LevelData<FluxBox> aux(grids, ncomp, ghostVect);

  BlockRegister br;
  br.define(coordSys, grids, numGhost);

  LayoutData<int> blockNumber(grids);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = grids[dit];
      int thisBlock = coordSys->whichBlock(baseBox);
      blockNumber[dit] = thisBlock;
    }

  for (dit.begin(); dit.ok(); ++dit)
    {
      int blockID = blockNumber[dit];
      const NewCoordSys* thisCoordSysPtr = coordSys->getCoordSys(blockID);

      FluxBox& fluxFlub = flux[dit];
      FluxBox& auxFlub = aux[dit];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox& fluxFab = fluxFlub[idir];
          FArrayBox& auxFab = auxFlub[idir];

          const Box& fluxBox = fluxFab.box();
          FArrayBox xiFab(fluxBox, SpaceDim);
          // returns face-centered xiFab
          thisCoordSysPtr->getCenterMappedCoordinates(xiFab, fluxBox);

          FArrayBox fluxScalarFab(scalarIntv, fluxFab);
          FArrayBox fluxVectorFab(vectorIntv, fluxFab);
          evalScalarFunction(fluxScalarFab, xiFab, thisCoordSysPtr);
          evalVectorFunction(fluxVectorFab, xiFab, thisCoordSysPtr);

          // Note that we store aux only on valid cells.
          // The BlockRegister should fill in within-block ghost cells.
          FArrayBox auxScalarFab(scalarIntv, auxFab);
          FArrayBox auxVectorFab(vectorIntv, auxFab);
          evalScalarFunction(auxScalarFab, xiFab, thisCoordSysPtr);
          evalVectorFunction(auxVectorFab, xiFab, thisCoordSysPtr);

          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();
              if (br.hasInterface(dit(), idir, side))
                {
                  br.storeFlux(fluxFab, dit(), idir, side);
                  br.storeAux(auxFab, dit(), idir, side);
                }
            }
        }
    }
  br.close();
  aux.exchange();

  const int MAXNORM = 0;

  Real scalarMaxNorm = norm(flux, scalarIntv, MAXNORM);
  Real acceptableScalar = 12. * precision * scalarMaxNorm;

  RealVect vectorMaxNorm;
  RealVect acceptableVector;
  for (int comp = vectorIntv.begin(), compBase = 0; comp <= vectorIntv.end(); comp++, compBase++)
    {
      Interval compIntvl(comp, comp);
      vectorMaxNorm[compBase] = norm(flux, compIntvl, MAXNORM);
      acceptableVector[compBase] = 12. * precision * vectorMaxNorm[compBase];
    }

  if (verbose)
    {
      pout() << "Scalar maxnorm=" << scalarMaxNorm
             << " so acceptable error is " << acceptableScalar << endl;
      pout() << "Vector maxnorm=" << vectorMaxNorm
             << " so acceptable error is " << acceptableVector << endl;
    }
  for (dit.begin(); dit.ok(); ++dit)
    {
      int blockID = blockNumber[dit];
      const NewCoordSys* thisCoordSysPtr = coordSys->getCoordSys(blockID);
      const CubedSphereShellPanelCS* thisPanelCS =
        dynamic_cast<const CubedSphereShellPanelCS*>(thisCoordSysPtr);

      const Tuple<BlockBoundary, 2*SpaceDim>& thisBoundaries =
        boundaries[blockID];
      const Box& blockBox = blockBoxes[blockID];

      const Box& baseBox = grids[dit];
        
      const FluxBox& fluxFlub = flux[dit];
      const FluxBox& auxFlub = aux[dit];
      int faceID = 0;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const FArrayBox& fluxFab = fluxFlub[idir];
              const FArrayBox& auxFab = auxFlub[idir];

              if (br.hasInterface(dit(), idir, side))
                {
                  // maybe better if this is done inside BlockRegister
                  const BlockBoundary& bb = thisBoundaries[faceID];
                  int reorientFace = bb.reorientFace(idir);
                  int idirOther = bb.dirOther(idir);
                  Side::LoHiSide sideOther = flip(side);
                  int blockOther = bb.neighbor();

                  Side::LoHiSide sideOtherTransformed = sideOther;
                  if (reorientFace == -1) sideOtherTransformed = side;

                  Box faceBlockBox = adjCellBox(blockBox, idir, side, 1);
                  faceBlockBox.shiftHalf(idir, -sign(side));

                  Box faceBox = adjCellBox(baseBox, idir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(idir, -sign(side));

                  if (verbose)
                    {
                      pout() << "Box " << faceBox
                             << " of block " << blockID
                             << " " << nameFace(idir, side)
                             << " to block " << blockOther
                             << " " << nameFace(idirOther, sideOtherTransformed);
                      if (reorientFace == -1) pout() << " reorient";
                      pout() << endl;
                    }

                  IntVect growTan = numGhost * (IntVect::Unit - BASISV(idir));
                  Box faceGrownBox = grow(faceBox, growTan);
                  faceGrownBox &= faceBlockBox;

                  FArrayBox fluxThisFab(faceBox, ncomp);
                  FArrayBox fluxOtherFab(faceBox, ncomp);
                  FArrayBox auxThisFab(faceGrownBox, ncomp);
                  FArrayBox auxOtherFab(faceGrownBox, ncomp);

                  // Note the quirk that
                  // we need to ask for sideOther to get this side,
                  // and to ask for side to get other side.
                  br.getFlux(fluxThisFab, dit(), idir, side, sideOther);
                  br.getFlux(fluxOtherFab, dit(), idir, side, side);
                  br.getAux(auxThisFab, dit(), idir, side, sideOther);
                  br.getAux(auxOtherFab, dit(), idir, side, side);

                  /*
                    Compare correct flux and aux against what is
                    stored on the face.
                   */

                  FArrayBox diffFluxThisFab(faceBox, ncomp);
                  diffFluxThisFab.copy(fluxThisFab);
                  diffFluxThisFab.minus(fluxFab, 0, 0, ncomp);
                  bad = false;
                  Vector<Real> maxnormDiffFluxThis(ncomp);
                  for (int comp = scalarIntv.begin(); comp <= scalarIntv.end(); comp++)
                    {
                      maxnormDiffFluxThis[comp] =
                        diffFluxThisFab.norm(MAXNORM, comp);
                      worstRatio = Max(worstRatio,
                                       maxnormDiffFluxThis[comp]/scalarMaxNorm);
                      // should be 0 because we're just copying.
                      if (maxnormDiffFluxThis[comp] > 0.) bad = true;
                    }
                  for (int comp = vectorIntv.begin(), compBase = 0; comp <= vectorIntv.end(); comp++, compBase++)
                    {
                      maxnormDiffFluxThis[comp] =
                        diffFluxThisFab.norm(MAXNORM, comp);
                      worstRatio = Max(worstRatio,
                                       maxnormDiffFluxThis[comp]/vectorMaxNorm[compBase]);
                      // should be 0 because we're just copying.
                      if (maxnormDiffFluxThis[comp] > 0.) bad = true;
                    }
                  if (bad) status--;
                  if (showSame)
                    {
                      pout() << "maxnorm(diffFluxThisFab) = "
                             << maxnormDiffFluxThis;
                      if (bad) pout() << indent2 << indent2 << "<= N.B.";
                      pout() << endl;
                    }

                  FArrayBox diffAuxThisFab(faceGrownBox, ncomp);
                  diffAuxThisFab.copy(auxThisFab);
                  diffAuxThisFab.minus(auxFab, 0, 0, ncomp);
                  Vector<Real> maxnormDiffAuxThis(ncomp);
                  for (int comp = scalarIntv.begin(); comp <= scalarIntv.end(); comp++)
                    {
                      maxnormDiffAuxThis[comp] =
                        diffAuxThisFab.norm(MAXNORM, comp);
                      worstRatio = Max(worstRatio,
                                       maxnormDiffAuxThis[comp]/scalarMaxNorm);
                      // should be 0 because we're just copying.
                      if (maxnormDiffAuxThis[comp] > 0.) bad = true;
                    }
                  for (int comp = vectorIntv.begin(), compBase = 0; comp <= vectorIntv.end(); comp++, compBase++)
                    {
                      maxnormDiffAuxThis[comp] =
                        diffAuxThisFab.norm(MAXNORM, comp);
                      worstRatio = Max(worstRatio,
                                       maxnormDiffAuxThis[comp]/vectorMaxNorm[compBase]);
                      // should be 0 because we're just copying.
                      if (maxnormDiffAuxThis[comp] > 0.) bad = true;
                    }
                  if (bad) status--;
                  if (showSame)
                    {
                      pout() << "maxnorm(diffAuxThisFab) = "
                             << maxnormDiffAuxThis;
                      if (bad) pout() << indent2 << indent2 << "<= N.B.";
                      pout() << endl;
                    }

                  /*
                    Compare scalar fluxes on the face.
                  */

                  // Reorient SCALAR fluxes (but not vector fluxes) if necessary.
                  if (reorientFace == -1)
                    {
                      for (int iflux = 0; iflux < ncomp; iflux++)
                        {
                          if (!vectorIntv.contains(iflux))
                            {
                              fluxOtherFab.negate(iflux, 1);
                              auxOtherFab.negate(iflux, 1);
                            }
                        }
                    }

                  FArrayBox diffFluxScalarFab(faceBox, nscalar);
                  diffFluxScalarFab.copy(fluxThisFab, iscalar, 0, nscalar);
                  diffFluxScalarFab.minus(fluxOtherFab, iscalar, 0, nscalar);
                  Real maxnormDiffFluxScalar =
                    diffFluxScalarFab.norm(MAXNORM, 0, nscalar);
                  worstRatio =
                    Max(worstRatio, maxnormDiffFluxScalar/acceptableScalar);
                  bad = (maxnormDiffFluxScalar > acceptableScalar);
                  if (bad) status--;
                  if (showFlux && showScalar)
                    {
                      pout() << "maxnorm(diffFluxScalarFab) = "
                             << maxnormDiffFluxScalar;
                      if (bad) pout() << indent2 << indent2 << "<= N.B.";
                      pout() << endl;
                    }
                      
                  FArrayBox diffAuxScalarFab(faceGrownBox, nscalar);
                  diffAuxScalarFab.copy(auxThisFab, iscalar, 0, nscalar);
                  diffAuxScalarFab.minus(auxOtherFab, iscalar, 0, nscalar);
                  Real maxnormDiffAuxScalar =
                    diffAuxScalarFab.norm(MAXNORM, 0, nscalar);
                  worstRatio =
                    Max(worstRatio, maxnormDiffAuxScalar/acceptableScalar);
                  bad = (maxnormDiffAuxScalar > acceptableScalar);
                  if (bad) status--;
                  if (showAux && showScalar)
                    {
                      pout() << "maxnorm(diffAuxScalarFab) = "
                             << maxnormDiffAuxScalar;
                      if (bad) pout() << indent2 << indent2 << "<= N.B.";
                      pout() << endl;
                    }

                  /*
                    Compare vector fluxes on the face.
                  */

                  FArrayBox fluxVectorThisFab(vectorIntv, fluxThisFab);
                  FArrayBox fluxVectorOtherFab(vectorIntv, fluxOtherFab);
                  FArrayBox auxVectorThisFab(vectorIntv, auxThisFab);
                  FArrayBox auxVectorOtherFab(vectorIntv, auxOtherFab);

                  IndicesTransformation it = bb.getTransformation();
                  Box faceOtherBox = it.transform(faceBox);
                  Box faceOtherGrownBox = it.transform(faceGrownBox);

                  /*
                    Orthonormalize vector flux on this flux,
                    just so we can compare it in debugger.
                   */
                  FArrayBox fluxVectorThisOrthoFab(faceBox, SpaceDim);
                  fluxVectorThisOrthoFab.copy(fluxVectorThisFab);
                  thisPanelCS->
                    orthonormalizeVectorFluxes(fluxVectorThisOrthoFab,
                                               faceBox,
                                               baseVectorIntv);

                  /*
                    Transform the vector fluxes from the other block.
                   */
                  FArrayBox fluxVectorOtherTransFab(faceOtherBox, SpaceDim);
                  FArrayBox auxVectorOtherTransFab(faceOtherGrownBox, SpaceDim);
                  it.transformFwd(fluxVectorOtherTransFab,
                                  fluxVectorOtherFab,
                                  faceBox,
                                  baseVectorIntv, baseVectorIntv);
                  it.transformFwd(auxVectorOtherTransFab,
                                  auxVectorOtherFab,
                                  faceGrownBox,
                                  baseVectorIntv, baseVectorIntv);
                  const NewCoordSys* otherCoordSysPtr =
                    coordSys->getCoordSys(blockOther);
                  const CubedSphereShellPanelCS* otherPanelCS =
                    dynamic_cast<const CubedSphereShellPanelCS*>(otherCoordSysPtr);
                  otherPanelCS->
                    orthonormalizeVectorFluxes(fluxVectorOtherTransFab,
                                               faceOtherBox,
                                               baseVectorIntv);
                  otherPanelCS->
                    orthonormalizeVectorFluxes(auxVectorOtherTransFab,
                                               faceOtherGrownBox,
                                               baseVectorIntv);

                  // Reorient orthonormalized flux and aux in each
                  // tangential direction where necessary.
                  // (Never necessary in normal direction.)
                  int orthoDir = 0; // component in orthonormal frame
                  for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
                    if (tanDir != idir)
                      {
                        orthoDir++;
                        int sgnFace = bb.reorientFace(tanDir);
                        if (sgnFace * reorientFace == -1)
                          {
                            fluxVectorOtherTransFab.negate(orthoDir, 1);
                            auxVectorOtherTransFab.negate(orthoDir, 1);
                          }
                      }

                  it.transformBack(fluxVectorOtherFab,
                                   fluxVectorOtherTransFab,
                                   faceOtherBox,
                                   baseVectorIntv, baseVectorIntv);
                  it.transformBack(auxVectorOtherFab,
                                   auxVectorOtherTransFab,
                                   faceOtherGrownBox,
                                   baseVectorIntv, baseVectorIntv);
                  thisPanelCS->
                    deorthonormalizeVectorFluxes(fluxVectorOtherFab,
                                                 faceBox,
                                                 baseVectorIntv);
                  thisPanelCS->
                    deorthonormalizeVectorFluxes(auxVectorOtherFab,
                                                 faceGrownBox,
                                                 baseVectorIntv);

                  FArrayBox diffFluxVectorFab(faceBox, SpaceDim);
                  diffFluxVectorFab.copy(fluxVectorThisFab, 0, 0, SpaceDim);
                  diffFluxVectorFab.minus(fluxVectorOtherFab, 0, 0, SpaceDim);
                  RealVect maxnormDiffFluxVector;
                  bad = false;
                  for (int comp = 0; comp < SpaceDim; comp++)
                    {
                      maxnormDiffFluxVector[comp] =
                        diffFluxVectorFab.norm(MAXNORM, comp, 1);
                      worstRatio =
                        Max(worstRatio, maxnormDiffFluxVector[comp]/acceptableVector[comp]);
                      bad |= (maxnormDiffFluxVector[comp] > acceptableVector[comp]);
                    }
                  if (bad) status--;
                  if (showFlux && showVector)
                    {
                      pout() << "maxnorm(diffFluxVectorFab) = "
                             << maxnormDiffFluxVector;
                      if (bad) pout() << indent2 << indent2 << "<= N.B.";
                      pout() << endl;
                    }

                  FArrayBox diffAuxVectorFab(faceGrownBox, SpaceDim);
                  diffAuxVectorFab.copy(auxVectorThisFab, 0, 0, SpaceDim);
                  diffAuxVectorFab.minus(auxVectorOtherFab, 0, 0, SpaceDim);
                  RealVect maxnormDiffAuxVector;
                  bad = false;
                  for (int comp = 0; comp < SpaceDim; comp++)
                    {
                      maxnormDiffAuxVector[comp] =
                        diffAuxVectorFab.norm(MAXNORM, comp, 1);
                      worstRatio =
                        Max(worstRatio, maxnormDiffAuxVector[comp]/acceptableVector[comp]);
                      bad |= (maxnormDiffAuxVector[comp] > acceptableVector[comp]);
                    }
                  if (bad) status--;
                  if (showAux && showVector)
                    {
                      pout() << "maxnorm(diffAuxVectorFab) = "
                             << maxnormDiffAuxVector;
                      if (bad) pout() << indent2 << indent2 << "<= N.B.";
                      pout() << endl;
                    }
                }
              faceID++;
            }
        }
    }

  if (verbose)
    {
      pout() << "worst ratio is " << worstRatio << endl;
    }
  // Clean up.
  delete coordSys;
  delete cssFactory;

  // Git.
  return status;
}

//-----------------------------------------------------------------------
int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = 0;

#if CH_SPACEDIM >= 2

  status = cubedSphereTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereTest failed with return code " 
      << status << endl ;

#endif

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize ();
#endif
  return status ;
}
//-----------------------------------------------------------------------
