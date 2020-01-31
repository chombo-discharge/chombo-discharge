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
using std::cout;
using std::endl;
using std::ofstream;
#include <limits> // contains numeric_limits<Real>

#include "parstream.H"
#include "MultiBlockFluxRegister.H"
#include "SingleBlockCSAdaptor.H"
#include "CartesianBlockCS.H"
#include "TripleCartesianCS.H"
#include "MultiCartesianCS.H"
#include "CubedSphereShellCS.H"
#include "CubedSphereShellPanelCS.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "LoadBalance.H"
#include "SphericalHarmonic.H"
#include "ClimateIBCF_F.H"
#include "MBMiscUtil.H"
#include "FourthOrderUtil.H"
#include "computeNorm.H"
#include "DebugOut.H"

#include "FABView.H"
#include "DebugDump.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testMultiBlockFluxRegister";
static const char* indent = "   ";
static const char* indent2 = "      ";
static bool verbose = false ;
static int verbosity = 0;

static Real precision = 2 * numeric_limits<Real>::epsilon();

static int MAXNORM = 0;

static int subRatio = 4;
static IntVect subRatioVect = subRatio * IntVect::Unit;

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
              // Check whether next argument is a number.
              // If it is, set verbosity to that.
              if ( i + 1 < argc )
                {
                  int num = atoi(argv[i+1]);
                  if (num > 0) verbosity = num;
                  pout() << "verbosity = " << verbosity << endl;
                }
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
void reduceMax(Real& a_val)
{
  Real recv = a_val;
#ifdef CH_MPI
  int resultMPI = MPI_Allreduce(&a_val, &recv,
                                1, MPI_CH_REAL,
                                MPI_MAX, Chombo_MPI::comm);

  if (resultMPI != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on reduceMax");
    }
#endif
  a_val = recv;
}

//-----------------------------------------------------------------------

Vector< List<int> > getCubedSphereCoarseBlockIndices()
{
  Vector< List<int> > coarseBlockIndices;
  List<int> cblist;
  // Blocks 0, 2, 3, 4, 5 are coarse; block 1 is fully refined.
  cblist.clear();
  cblist.append(0);
  cblist.append(2);
  cblist.append(3);
  cblist.append(4);
  cblist.append(5);
  coarseBlockIndices.push_back(cblist);
  // Blocks 1, 2, 3, 4, 5 are coarse; block 0 is fully refined.
  cblist.clear();
  cblist.append(1);
  cblist.append(2);
  cblist.append(3);
  cblist.append(4);
  cblist.append(5);
  coarseBlockIndices.push_back(cblist);
  // Blocks 0 and 2 are coarse; all other blocks are fully refined.
  cblist.clear();
  cblist.append(0);
  cblist.append(2);
  coarseBlockIndices.push_back(cblist);
  // Blocks 1 and 3 are coarse; all other blocks are fully refined.
  cblist.clear();
  cblist.append(1);
  cblist.append(3);
  coarseBlockIndices.push_back(cblist);
  // Blocks 4 and 5 are coarse; all other blocks are fully refined.
  cblist.clear();
  cblist.append(4);
  cblist.append(5);
  coarseBlockIndices.push_back(cblist);
  // Blocks 1, 3, 4, 5 are coarse; all other blocks are fully refined.
  cblist.clear();
  cblist.append(1);
  cblist.append(3);
  cblist.append(4);
  cblist.append(5);
  coarseBlockIndices.push_back(cblist);
  // Blocks 0, 2, 4, 5 are coarse; all other blocks are fully refined.
  cblist.clear();
  cblist.append(0);
  cblist.append(2);
  cblist.append(4);
  cblist.append(5);
  coarseBlockIndices.push_back(cblist);
  // Blocks 0, 1, 2, 3 are coarse; all other blocks are fully refined.
  cblist.clear();
  cblist.append(0);
  cblist.append(1);
  cblist.append(2);
  cblist.append(3);
  coarseBlockIndices.push_back(cblist);
  return coarseBlockIndices;
}

//-----------------------------------------------------------------------
int
singleBlockTest()
{
  // NOTE: This test demonstrates that the multi-block flux register can
  // NOTE: perform as a single-block flux register. Note that
  // NOTE: the coarse-fine grids here are properly nested in this test.
  // NOTE: This is not so with the multi-block tests, since the boxes
  // NOTE: in multi-block calculations are necessarily disjoint. So the
  // NOTE: logic is different for single- and multi-block cases.

  // Create a single-block Cartesian coordinate system.
  RealVect dxCoarse = 1.0*RealVect::Unit;
  Box domainBox(IntVect::Zero, 16*IntVect::Unit);
  ProblemDomain coarseDomain(domainBox);
  NewCoordSys* cartesianCS = new CartesianBlockCS(0, IntVect::Zero, dxCoarse, domainBox);
  // Note that coordSys takes control of cartesianCS!
  MultiBlockCoordSys* coordSys = new SingleBlockCSAdaptor(cartesianCS, coarseDomain);

  // Create the coarse-level box (two grid cells at (7,7), (8,7)).
  Vector<Box> coarseBoxes;
  IntVect Ic = 7*IntVect::Unit, Icf;
  Icf = Ic, Icf[0] += 1;
  Box coarseBox(Ic, Icf);
  coarseBoxes.push_back(coarseBox);
  Vector<int> procs(1,0);
  DisjointBoxLayout coarseGrids(coarseBoxes, procs, coarseDomain);

  // Create the fine-level box: a 2 x 2 square spanning (14, 14) x (15, 15).
  Vector<Box> fineBoxes;
  IntVect If = 14*IntVect::Unit;
  Box fineBox(If, If + IntVect::Unit);
  fineBoxes.push_back(fineBox);
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(2);
  DisjointBoxLayout fineGrids(fineBoxes, procs, fineDomain);

  Box fineDomainBox(domainBox);
  fineDomainBox.refine(2);
  NewCoordSys* cartesianFineCS = new CartesianBlockCS(0, IntVect::Zero, dxCoarse/2., fineDomainBox);
  // Note that coordSysFine takes control of cartesianFineCS!
  MultiBlockCoordSys* coordSysFine = new SingleBlockCSAdaptor(cartesianFineCS, fineDomain);

  // Create the flux register.
  MultiBlockFluxRegister fluxReg(coordSys, coordSysFine);
  fluxReg.define(fineGrids, coarseGrids, coarseDomain, 2, 1, true);
  fluxReg.setToZero();

  // Create a zero solution U the coarse grids with a single ghost layer.
  LevelData<FArrayBox> Ucoarse(coarseGrids, 1, IntVect::Unit);
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      Ucoarse[dit()].setVal(0.0);
    }

  // Create coarse and fine fluxes that straddle the coarse-fine boundary.
  LevelData<FluxBox> coarseFlux(coarseGrids, 1),
                     fineFlux(fineGrids, 1);
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      coarseFlux[dit()].setVal(1.0);
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      fineFlux[dit()].setVal(0.6);
    }

  // Stick the coarse and fine fluxes into the register to compute the defect.
  Interval i(0, 0);
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      fluxReg.incrementCoarse(coarseFlux[dit()][0], 1.0, dit(), i, i, 0);
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      fluxReg.incrementFine(fineFlux[dit()][0], 1.0, dit(), i, i, 0, Side::Hi);
    }

  // Apply the defect to the solution.
//  fluxReg.poutCoarseRegisters();
//  fluxReg.poutFineRegisters();
  //vector interval is weird
  Interval vectorInt;
  fluxReg.reflux(Ucoarse, i, vectorInt, i, 1.0);

  // Check the answer.
  int status = 0;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox& U = Ucoarse[dit()];
    IntVect Ig = Ic;
    Ig[0] += 1;
    if (Abs(U(Ig, 0) - 0.4) > precision)
      {
        status = -1;
      }
  }

  // Clean up.
  delete coordSys;
  delete coordSysFine;

  // Git.
  return status;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
tripleCartesianTest()
{
  // Create a triple-Cartesian coordinate system.
  // This creates a 3 x 3 set of N/5 x N/5 patches in indexspace,
  // separated by buffers of width N/5.
  RealVect dxCoarse = 1.0*RealVect::Unit;
  int coarseBoxlen = 16;
  Box domainBox(IntVect::Zero, (5*coarseBoxlen - 1)*IntVect::Unit);
  ProblemDomain coarseDomain(domainBox);
  TripleCartesianCS* coordSysCoarse = new TripleCartesianCS();
  coordSysCoarse->define(coarseDomain, dxCoarse);

  int fluxDir = 0;

  // Create the coarse-level box containing a grid cell at (32,40),
  // abutting the boundary between the center-left and center block.
  Vector<Box> coarseBoxes;
  IntVect Icoarse = (2*coarseBoxlen + coarseBoxlen/2)*IntVect::Unit;
  Icoarse[fluxDir] = 2*coarseBoxlen;
  Box coarseBox(Icoarse, Icoarse);
  coarseBoxes.push_back(coarseBox);
  Vector<int> procs(1,0);
  DisjointBoxLayout coarseGrids(coarseBoxes, procs, coarseDomain);

  // Create the fine-level box, which is a 2 x 2 patch spanning (30, 80) x (31 x 81).
  // This puts the patch at (15, 40) in the coarse index space, which is in a
  // different block from the coarse-level box (but abuts it in the multi-
  // block sense).
  int refRatio = 2;
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatio);
  RealVect dxFine = dxCoarse / RealVect(refRatio*IntVect::Unit);

  Vector<Box> fineBoxes;
  IntVect Ifine = refRatio*Icoarse;
  Ifine[fluxDir] = refRatio*(coarseBoxlen-1);
  Box fineBox(Ifine, Ifine + (refRatio-1)*IntVect::Unit);
  fineBoxes.push_back(fineBox);
  DisjointBoxLayout fineGrids(fineBoxes, procs, fineDomain);

  TripleCartesianCS* coordSysFine = new TripleCartesianCS();
  coordSysFine->define(fineDomain, dxFine);

  // Create the flux register.
  MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
  int ncomp = 1;
  bool scaleFlag = true;
  fluxReg.define(fineGrids, coarseGrids, coarseDomain, refRatio, ncomp, scaleFlag);
  fluxReg.setToZero();

  // Create a constant solution U the coarse grids with a single ghost layer.
  Real startVal = 0.1;
  LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      Ucoarse[dit].setVal(startVal);
    }

  // Create coarse and fine fluxes that straddle the coarse-fine boundary.
  LevelData<FluxBox> coarseFlux(coarseGrids, ncomp),
                     fineFlux(fineGrids, ncomp);
  Real coarseVal = 1.0;
  Real fineVal = 0.6;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      coarseFlux[dit].setVal(coarseVal);
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      fineFlux[dit].setVal(fineVal);
    }

  // Stick the coarse and fine fluxes into the register to compute the defect.
  Interval zeroIntvl(0, 0);
  Real scale1 = 1.0;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      fluxReg.incrementCoarse(coarseFlux[dit][fluxDir], scale1, dit(), zeroIntvl, zeroIntvl, fluxDir);
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      fluxReg.incrementFine(fineFlux[dit][fluxDir], scale1, dit(), zeroIntvl, zeroIntvl, fluxDir, Side::Hi);
    }

  // Apply the defect to the solution.
//  fluxReg.poutCoarseRegisters();
//  fluxReg.poutFineRegisters();
  Interval vectorIntvl;
  fluxReg.reflux(Ucoarse, zeroIntvl, vectorIntvl, zeroIntvl, scale1);

  // Check the answer.
  int status = 0;
  Real endVal = startVal + coarseVal - fineVal;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox& U = Ucoarse[dit];
    IntVect Ig = Icoarse;
//cout << "U(" << Ig << ", 0) = " << U(Ig, 0) << endl;
    if (Abs(U(Ig, 0) - endVal) > precision)
      {
        status = -1;
      }
  }

  // Clean up.
  delete coordSysCoarse;
  delete coordSysFine;

  // Git.
  return status;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
multiCartesianTest()
{
  // Create a multi-Cartesian coordinate system 
  // with 2 blocks in other-than-highest dimension, with wraparound,
  // and 1 block in highest dimension, no wraparound.
  // This creates a (2 x) 2 x 1 set of (N x) N x N patches in indexspace,
  // separated by buffers of width N.
  RealVect dxCoarse = 1.0*RealVect::Unit;
  int coarseBoxlen = 16;

  IntVect count = 2*IntVect::Unit;
  if (SpaceDim > 1)
    {
      count[SpaceDim-1] = 1;
    }
  IntVect domainLength = coarseBoxlen * (2*count - IntVect::Unit);
  Box domainBox(IntVect::Zero, domainLength - IntVect::Unit);
  ProblemDomain coarseDomain(domainBox);

  int fluxDir = 0;

  // Wrap around in every dimension except the last one.
  IntVect wraparound = IntVect::Unit - BASISV(SpaceDim-1);
  MultiCartesianCS* coordSysCoarse = new MultiCartesianCS(count, wraparound);
  coordSysCoarse->define(coarseDomain, dxCoarse);

  // Create the coarse-level box containing a grid cell at (32, 8, 8),
  // abutting the boundary between two blocks.
  Vector<Box> coarseBoxes;
  IntVect Icoarse = (coarseBoxlen/2)*IntVect::Unit;
  Icoarse[fluxDir] = 2*coarseBoxlen;
  Box coarseBox(Icoarse, Icoarse);
  coarseBoxes.push_back(coarseBox);
  Vector<int> procs(1,0);
  DisjointBoxLayout coarseGrids(coarseBoxes, procs, coarseDomain);

  int refRatio = 2;
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatio);
  RealVect dxFine = dxCoarse / RealVect(refRatio*IntVect::Unit);

  // Create the fine-level box, which is a 2 x 2 x 2 patch
  // spanning (30:31, 16:17, 16:17) in fine index space.
  // This puts the patch at (15, 8, 8) in coarse index space, which is in a
  // different block from the coarse-level box (but abuts it in the multi-
  // block sense).
  Vector<Box> fineBoxes;
  IntVect Ifine = refRatio*Icoarse;
  Ifine[fluxDir] = refRatio*(coarseBoxlen-1);
  Box fineBox(Ifine, Ifine + (refRatio-1)*IntVect::Unit);
  fineBoxes.push_back(fineBox);
  DisjointBoxLayout fineGrids(fineBoxes, procs, fineDomain);

  MultiCartesianCS* coordSysFine = new MultiCartesianCS(count, wraparound);
  coordSysFine->define(fineDomain, dxFine);

  // Create the flux register.
  MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
  int ncomp = 1;
  bool scaleFlag = true;
  fluxReg.define(fineGrids, coarseGrids, coarseDomain, refRatio, ncomp, scaleFlag);
  fluxReg.setToZero();

  // Create a constant solution U the coarse grids with a single ghost layer.
  Real startVal = 0.1;
  LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      Ucoarse[dit].setVal(startVal);
    }

  // Create coarse and fine fluxes that straddle the coarse-fine boundary.
  LevelData<FluxBox> coarseFlux(coarseGrids, 1),
                     fineFlux(fineGrids, 1);
  Real coarseVal = 1.0;
  Real fineVal = 0.6;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      coarseFlux[dit].setVal(coarseVal);
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      fineFlux[dit].setVal(fineVal);
    }

  // Stick the coarse and fine fluxes into the register to compute the defect.
  Interval zeroIntvl(0, 0);
  Real scale1 = 1.0;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      fluxReg.incrementCoarse(coarseFlux[dit][fluxDir], scale1, dit(), zeroIntvl, zeroIntvl, fluxDir);
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      fluxReg.incrementFine(fineFlux[dit][fluxDir], scale1, dit(), zeroIntvl, zeroIntvl, fluxDir, Side::Hi);
    }

  // Apply the defect to the solution.
//  fluxReg.poutCoarseRegisters();
//  fluxReg.poutFineRegisters();
  Interval vectorIntvl;
  fluxReg.reflux(Ucoarse, zeroIntvl, vectorIntvl, zeroIntvl, scale1);

  // Check the answer.
  int status = 0;
  Real endVal = startVal + coarseVal - fineVal;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& U = Ucoarse[dit];
      IntVect Ig = Icoarse;
      //cout << "U(" << Ig << ", 0) = " << U(Ig, 0) << endl;
      if (Abs(U(Ig, 0) - endVal) > precision)
        {
          status = -1;
        }
    }

  // Clean up.
  delete coordSysCoarse;
  delete coordSysFine;

  return status;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
cubedSphereCellTest()
{
  // Create a cubed sphere coordinate system with six N-by-N blocks.
  int nlayers = 8;
  int N = 10;
  int nRefRatio = 2;

  int interpDims = SpaceDim;
  if (SpaceDim == 3) interpDims = SpaceDim-1;

  int layer = 0; // which layer in 3D

  CubedSphereShellCS* coordSys = new CubedSphereShellCS(N, nlayers);
  IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
  ProblemDomain coarseDomain(IntVect::Zero, domainLen-IntVect::Unit);

  IntVect refRatioVect = IntVect::Unit;
  for (int idir = 0; idir < interpDims; idir++)
    {
      refRatioVect[idir] = nRefRatio;
    }
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatioVect);

  CubedSphereShellCS* coordSysFine = new CubedSphereShellCS(N/nRefRatio, nlayers);

  const Vector<Box>& blockBoxesCoarse =
    coordSys->mappingBlocks();
  //  const Vector<Box>& blockBoxesFine =
  //    coordSysFine->mappingBlocks();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSys->boundaries();

  int status = 0;
  int ncomp = 1;
  Interval allIntv(0, ncomp-1);
  Interval vectorIntv; // empty

  int nblocks = blockBoxesCoarse.size();
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      const Box& blockBox = blockBoxesCoarse[iblock];
      const IntVect& blockBoxLo = blockBox.smallEnd();
      const Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries =
        boundaries[iblock];
      int faceID = 0;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir < interpDims)
                {
                  const BlockBoundary& bb = blockBoundaries[faceID];
                  IndicesTransformation it = bb.getTransformation();
                  int idirOther = bb.dirOther(idir);
                  int reorient = bb.reorientFace(idir);
                  Side::LoHiSide sideOther =
                    (reorient == -1) ? side : flip(side);
                  
                  // Coarse-level box containing a grid cell in middle of
                  // this face of block.
                  Box adjBox = adjCellBox(blockBox, idir, side, 1);
                  for (int idirTan = 0; idirTan < interpDims; idirTan++)
                    if (idirTan != idir)
                      {
                        adjBox.setRange(idirTan, blockBoxLo[idirTan]+N/2);
                      }
                  if (SpaceDim == 3) adjBox.setRange(SpaceDim-1, layer);
                  // Now adjBox consists of 1 cell, just outside block.
                  // Set coarseBox to 1 cell, just inside block.
                  Box coarseBox(adjBox);
                  coarseBox.shift(idir, -sign(side));
                  CH_assert(coarseBox.numPts() == 1);
                  IntVect coarseCell = coarseBox.smallEnd();
                  Vector<Box> coarseBoxes;
                  coarseBoxes.push_back(coarseBox);
                  
                  Box adjOtherBox = it.transform(adjBox);
                  Box fineBox = refine(adjOtherBox, refRatioVect);
                  Vector<Box> fineBoxes;
                  fineBoxes.push_back(fineBox);
                  
                  Vector<int> coarseProcs(coarseBoxes.size(), 0);
                  Vector<int> fineProcs(coarseBoxes.size(), 0);
                  DisjointBoxLayout coarseGrids(coarseBoxes, coarseProcs, coarseDomain);
                  DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
                  
                  // Create the flux register.
                  MultiBlockFluxRegister fluxReg(coordSys, coordSysFine);
                  // If SpaceDim = 3, don't refine in the r (last) direction (anisotropic refinement)
                  int unrefDir = (SpaceDim == 3) ? SpaceDim-1 : -1;
                  bool scaleFlag = true;
                  fluxReg.define(fineGrids, coarseGrids, coarseDomain,
                                 nRefRatio, ncomp, scaleFlag, unrefDir);
                  fluxReg.setToZero();
                  
                  // Create a solution Ucoarse with startVal on
                  // the coarse grids with a single ghost layer.
                  Real startVal = 0.1;
                  LevelData<FArrayBox> Ucoarse(coarseGrids, 1, IntVect::Unit);
                  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
                    {
                      Ucoarse[dit].setVal(startVal);
                    }
                  
                  // Create coarse and fine fluxes that straddle the
                  // coarse-fine boundary.
                  LevelData<FluxBox> coarseFlux(coarseGrids, 1);
                  LevelData<FluxBox> fineFlux(fineGrids, 1);
                  Real coarseVal = 1.0;
                  Real fineVal = 0.6;
                  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
                    {
                      coarseFlux[dit].setVal(coarseVal);
                    }
                  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
                    { // reorient fine value if necessary
                      fineFlux[dit].setVal(reorient * fineVal);
                    }
                  
                  // Stick the coarse and fine fluxes into the register
                  // to compute the defect.
                  Real scale1 = 1.0;
                  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
                    {
                      fluxReg.incrementCoarse(coarseFlux[dit][idir], scale1, dit(),
                                              allIntv, allIntv, idir, flip(side));
                    }
                  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
                    {
                      fluxReg.incrementFine(fineFlux[dit][idirOther], scale1, dit(),
                                            allIntv, allIntv, idirOther, sideOther);
                    }
                  
                  // Apply the defect to the solution.
                  fluxReg.reflux(Ucoarse, vectorIntv, scale1);
                  
                  // Check the answer.
                  Real endVal = startVal + sign(side)*(fineVal - coarseVal);
                  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
                    {
                      const FArrayBox& UFab = Ucoarse[dit];
                      Real UcoarseVal = UFab(coarseCell, 0);
#if 0
                      pout() << "Block face " << iblock
                             << ((idir == 0) ? "x" : "y")
                             << ((side == Side::Lo) ? "Lo" : "Hi")
                             << " to " << bb.neighbor()
                             << ((idirOther == 0) ? "x" : "y")
                             << ((sideOther == Side::Lo) ? "Lo" : "Hi")
                             << " reorient=" << reorient
                             << ": UcoarseVal = " << UcoarseVal << endl;
#endif
                      if (Abs(UcoarseVal - endVal) > precision)
                        {
                          status--;
                        }
                    }
                }
              faceID++;
            }
        }
    }
  
  // Clean up.
  delete coordSys;
  delete coordSysFine;

  // Git.
  return status;
}

//-----------------------------------------------------------------------
// Returns the low and high faces of the given box in the given dimension.
Tuple<Box, 2> adjFaces(const Box& a_bx,
                       int a_dir)
{
  Tuple<Box, 2> faces;
  faces[0] = adjCellLo(a_bx, a_dir);
  faces[0].shiftHalf(a_dir, 1);
  faces[1] = adjCellHi(a_bx, a_dir);
  faces[1].shiftHalf(a_dir, -1);
  return faces;
}

//-----------------------------------------------------------------------

void evalScalarFunction(FArrayBox& a_funFab,
                        const Box& a_bx,
                        int a_dirNormal,
                        const NewCoordSys* a_coordSysPtr,
                        const Box& a_blockBox)
{
  CH_assert(a_funFab.box().contains(a_bx));
  CH_assert(faceDimension(a_bx) == a_dirNormal);
  // blockFaces[0] and blockFaces[1] are
  // end faces of a_blockBox, low and high in direction a_dirNormal
  Tuple<Box, 2> blockFaces = adjFaces(a_blockBox, a_dirNormal);
  for (int iface = 0; iface < 2; iface++)
    {
      Box faceBox = a_bx & blockFaces[iface];
      if ( !faceBox.isEmpty() )
        {
          FArrayBox xiFab(faceBox, SpaceDim);
          a_coordSysPtr->getCenterMappedCoordinates(xiFab, faceBox);
          FArrayBox rllFab(faceBox, SpaceDim);
          const CubedSphereShellPanelCS* cssPanelPtr =
            dynamic_cast<const CubedSphereShellPanelCS*>(a_coordSysPtr);
          cssPanelPtr->fabTransformEquiangularToLonLat(xiFab, rllFab);
          Real magnitude = 2.0 * M_PI;
          Real shift = 0.;
          Real angle = M_PI / 3.0;
          FArrayBox vecRLLFab(faceBox, SpaceDim);
          FORT_CSPVECLONLATSOLIDBODYFLOW(CHF_FRA(vecRLLFab),
                                         CHF_CONST_FRA(rllFab),
                                         CHF_CONST_REAL(magnitude),
                                         CHF_CONST_REAL(angle),
                                         CHF_CONST_REAL(shift));
#if CH_SPACEDIM >= 3
          // vecRLLFab.setVal(1., 2);
          // Set vecRLLFab[2] = vecRLLFab[0]^2 + vecRLLFab[1]^2.
          FArrayBox velLon2Fab(faceBox, 1);
          velLon2Fab.copy(vecRLLFab, 0, 0);
          velLon2Fab.mult(vecRLLFab, 0, 0);
          FArrayBox velLat2Fab(faceBox, 1);
          velLat2Fab.copy(vecRLLFab, 1, 0);
          velLat2Fab.mult(vecRLLFab, 1, 0);
          vecRLLFab.copy(velLon2Fab, 0, 2);
          vecRLLFab.plus(velLat2Fab, 0, 2);
#endif
          FArrayBox vecEquFab(faceBox, SpaceDim);
          cssPanelPtr->fabVectorTransformLatLonToEquiangular(xiFab, vecRLLFab, vecEquFab);
          // Return normal component of velocity.
          a_funFab.copy(vecEquFab, a_dirNormal, 0);
        }
    }
}

//-----------------------------------------------------------------------
void evalVectorFunction(FArrayBox& a_funFab,
                        const Box& a_bx,
                        int a_dirNormal,
                        const NewCoordSys* a_coordSysPtr,
                        const Box& a_blockBox)
{
  CH_assert(a_funFab.box().contains(a_bx));
  // CH_assert(faceDimension(a_bx) == a_dirNormal);
  //  CH_assert(a_bx.ixType() == IndexType::TheNodeType());
  //  pout() << "evalVectorFunction box " << a_bx.ixType() << " normal=" << a_dirNormal << endl;
  if (a_bx.ixType() != IndexType::TheNodeType())
    {
      CH_assert(faceDimension(a_bx) == a_dirNormal);
    }
  // blockFaces[0] and blockFaces[1] are
  // end faces of a_blockBox, low and high in direction a_dirNormal
  Tuple<Box, 2> blockFaces = adjFaces(a_blockBox, a_dirNormal);
  for (int iface = 0; iface < 2; iface++)
    {
      Box faceBox = a_bx;
      if (a_bx.ixType() == IndexType::TheNodeType())
        {
          Box blockNodes = surroundingNodes(blockFaces[iface]);
          IntVect tan1Vect = IntVect::Unit - BASISV(a_dirNormal);
          blockNodes.grow(4*tan1Vect);
          faceBox &= blockNodes; // will be NODE-centered now
        }
      else
        { // dirNormal-face-centered 
          faceBox &= blockFaces[iface];
        }
      if ( !faceBox.isEmpty() )
        {
          FArrayBox xiFab(faceBox, SpaceDim);
          a_coordSysPtr->getCenterMappedCoordinates(xiFab, faceBox);
          FArrayBox rllFab(faceBox, SpaceDim);
          const CubedSphereShellPanelCS* cssPanelPtr =
            dynamic_cast<const CubedSphereShellPanelCS*>(a_coordSysPtr);
          cssPanelPtr->fabTransformEquiangularToLonLat(xiFab, rllFab);
          Real magnitude = 2.0 * M_PI;
          Real shift = 0.;
          Real angle = M_PI / 3.0;
          FArrayBox vecRLLFab(faceBox, SpaceDim);
          FORT_CSPVECLONLATSOLIDBODYFLOW(CHF_FRA(vecRLLFab),
                                         CHF_CONST_FRA(rllFab),
                                         CHF_CONST_REAL(magnitude),
                                         CHF_CONST_REAL(angle),
                                         CHF_CONST_REAL(shift));
#if CH_SPACEDIM >= 3
          // vecRLLFab.setVal(1., 2);
          // Set vecRLLFab[2] = vecRLLFab[0]^2 + vecRLLFab[1]^2.
          FArrayBox velLon2Fab(faceBox, 1); // set this to vecRLLFab[0]^2
          velLon2Fab.copy(vecRLLFab, 0, 0);
          velLon2Fab.mult(vecRLLFab, 0, 0);
          FArrayBox velLat2Fab(faceBox, 1); // set this to vecRLLFab[1]^2
          velLat2Fab.copy(vecRLLFab, 1, 0);
          velLat2Fab.mult(vecRLLFab, 1, 0);
          vecRLLFab.copy(velLon2Fab, 0, 2);
          vecRLLFab.plus(velLat2Fab, 0, 2);
#endif
          FArrayBox vecEquFab(faceBox, SpaceDim);
          cssPanelPtr->fabVectorTransformLatLonToEquiangular(xiFab, vecRLLFab, vecEquFab);
          // Return velocity vector multiplied by normal component of velocity.
          a_funFab.copy(vecEquFab);
          for (int comp = 0; comp < SpaceDim; comp++)
            {
              a_funFab.mult(vecEquFab, a_dirNormal, comp);
            }
        }
    }
}

//-----------------------------------------------------------------------

void evalVectorAvgFunction(FArrayBox& a_funFab,
                           const Box& a_bx,
                           int a_dirNormal,
                           const NewCoordSys* a_coordSysSubPtr,
                           const Box& a_blockSubBox)
{
  CH_assert(a_funFab.box().contains(a_bx));
  CH_assert(faceDimension(a_bx) == a_dirNormal);
  // int dirNormal = faceDimension(a_bx);
  // Use Boole's rule: weights (7, 32, 12, 32, 7)/90 on subnodes.
  Box bxSub = refine(a_bx, subRatioVect);
  Box bxSubNodes = surroundingNodes(bxSub);

  IntVect integrateNodesMax = subRatioVect;
  integrateNodesMax[a_dirNormal] = 0;
  Box integrateNodesBox(IntVect::Zero, integrateNodesMax,
                        IndexType::TheNodeType());
  Real booleCoeffs[] = {7./90., 32./90., 12./90., 32./90., 7./90.};
  FArrayBox integrateNodesFab(integrateNodesBox, 1);
  integrateNodesFab.setVal(1.);
  for (BoxIterator bitN(integrateNodesBox); bitN.ok(); ++bitN)
    {
      IntVect ivN = bitN();
      for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
        if (tanDir != a_dirNormal)
          {
            int ind = ivN[tanDir]; // will be in the range 0:4
            integrateNodesFab(ivN, 0) *= booleCoeffs[ind];
          }
    }

  Tuple<Box, 2> blockFacesSub = adjFaces(a_blockSubBox, a_dirNormal);
  IntVect tan1Vect = IntVect::Unit - BASISV(a_dirNormal);
  for (int iface = 0; iface < 2; iface++)
    {
      Box blockNodesSub = surroundingNodes(blockFacesSub[iface]);
      blockNodesSub.grow(4*tan1Vect);
      Box bxSubFaceNodes = bxSubNodes & blockNodesSub;
      if ( !bxSubFaceNodes.isEmpty() )
        {
          Box bxSubFaces = bxSubFaceNodes;
          for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
            if (tanDir != a_dirNormal)
              bxSubFaces.enclosedCells(tanDir);

          Box faceBox = coarsen(bxSubFaces, subRatioVect);
          // a_funFab.setVal(0., faceBox, 0, SpaceDim);

          FArrayBox nodesSubFunFab(bxSubFaceNodes, SpaceDim);
          evalVectorFunction(nodesSubFunFab, bxSubFaceNodes, a_dirNormal,
                             a_coordSysSubPtr, a_blockSubBox);
          for (BoxIterator bit(faceBox); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              IntVect ivSubLo = subRatioVect * iv;
              IntVect ivSubHi = ivSubLo + integrateNodesMax;
              Box faceSubNodes(ivSubLo, ivSubHi, IndexType::TheNodeType());
              FArrayBox nodesSubThisFaceFab(faceSubNodes, SpaceDim);
              nodesSubThisFaceFab.copy(nodesSubFunFab);
              nodesSubThisFaceFab.shift(-ivSubLo);
              for (int comp = 0; comp < SpaceDim; comp++)
                {
                  nodesSubThisFaceFab.mult(integrateNodesFab, 0, comp);
                  a_funFab(iv, comp) = nodesSubThisFaceFab.sum(comp);
                }
            }
          CH_assert(0 == 0); // dummy for gdb
        }
    }
}

//-----------------------------------------------------------------------
// Returns max norm of the restriction of face-centered a_fab to
// the end faces of a_blockBox.
Real maxnormFaceFab(const FArrayBox& a_fab,
                    const Box& a_blockBox)
{
  const Box& bx = a_fab.box();
  int ncomp = a_fab.nComp();
  int dirNormal = faceDimension(bx);
  // blockFaces[0] and blockFaces[1] are
  // end faces of a_blockBox, low and high in direction dirNormal
  Tuple<Box, 2> blockFaces = adjFaces(a_blockBox, dirNormal);
  Real normFace = 0.;
  for (int iface = 0; iface < 2; iface++)
    {
      Box bxFace = a_fab.box() & blockFaces[iface];
      if ( !bxFace.isEmpty() )
        {
          normFace = a_fab.norm(bxFace, MAXNORM, 0, ncomp);
        }
    }
  return normFace;
}

//-----------------------------------------------------------------------
int
cubedSphereWholeTest()
{
  // Create a cubed sphere coordinate system with six N-by-N blocks.
#if CH_SPACEDIM >= 3
  // Avoid irritating unused variable warning if dimension < 3.
  int nlayers = 8;
#endif
  int N = 32;
  int nRefRatio = 2;

  int interpDims = SpaceDim;
  if (SpaceDim == 3) interpDims = SpaceDim-1;

  //  int layer = 0; // which layer in 3D
  // CubedSphereShellCS* coordSys = new CubedSphereShellCS(N, nlayers);

  IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
  ProblemDomain coarseDomain(IntVect::Zero, domainLen-IntVect::Unit);

  IntVect refRatioVect = IntVect::Unit;
  for (int idir = 0; idir < interpDims; idir++)
    {
      refRatioVect[idir] = nRefRatio;
    }
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatioVect);

  CubedSphereShellCSFactory* cssFactory = new CubedSphereShellCSFactory();
  cssFactory->setHeight(1.);
  cssFactory->setRadius(1.);
  RealVect probLen = RealVect(D_DECL(11*(M_PI/2.), M_PI/2., 1.));
  RealVect dxCoarse = probLen / RealVect(domainLen);
  RealVect dxFine = dxCoarse / RealVect(refRatioVect);
  MultiBlockCoordSys* coordSysCoarse =
    cssFactory->getCoordSys(coarseDomain, dxCoarse);
  MultiBlockCoordSys* coordSysFine =
    cssFactory->getCoordSys(fineDomain, dxFine);

  const Vector<Box>& blockBoxesCoarse =
    coordSysCoarse->mappingBlocks();
  const Vector<Box>& blockBoxesFine =
    coordSysFine->mappingBlocks();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSysCoarse->boundaries();

  int nblocks = blockBoxesCoarse.size();

  Box refinementBox(IntVect::Zero, refRatioVect-IntVect::Unit);
  IntVect baseLen = IntVect(D_DECL(N, N, nlayers));
  Box baseBox(IntVect::Zero, baseLen-IntVect::Unit);

  Vector< List<int> > coarseBlockIndices = getCubedSphereCoarseBlockIndices();

  int ncomp = 1;
  Interval allIntv(0, ncomp-1);
  Interval vectorIntv; // empty

  Real worstRatio = 0.;

  int status = 0;
  for (int icbi = 0; icbi < coarseBlockIndices.size(); icbi++)
    {
      Vector<Box> coarseBoxes;
      Vector<Box> fineBoxes;
      const List<int>& coarseList = coarseBlockIndices[icbi];

      if (verbosity >= 2)
        {
          pout() << "Coarse blocks";
          for (ListIterator<int> li(coarseList); li.ok(); ++li)
            {
              int blockNum = coarseList[li];
              pout() << " " << blockNum;
            }
        }

      for (int iblock = 0; iblock < nblocks; iblock++)
        {
          const Box& blockBox = blockBoxesCoarse[iblock];
          const IntVect& blockBoxLo = blockBox.smallEnd();
          if (coarseList.includes(iblock))
            {
              coarseBoxes.push_back(blockBox);
            }
          else
            {
              IntVect fineBlockBoxLo = refRatioVect * blockBoxLo;
              for (BoxIterator bitRef(refinementBox); bitRef.ok(); ++bitRef)
                {
                  IntVect br = bitRef();
                  IntVect thisFineShift = baseLen * br + fineBlockBoxLo;
                  Box thisFineBox(baseBox);
                  thisFineBox.shift(thisFineShift);
                  fineBoxes.push_back(thisFineBox);
                }
            }
        }

      mortonOrdering(coarseBoxes);
      Vector<int> coarseProcs(coarseBoxes.size());
      LoadBalance(coarseProcs, coarseBoxes);
      DisjointBoxLayout coarseGrids(coarseBoxes, coarseProcs, coarseDomain);
      DataIterator ditC = coarseGrids.dataIterator();

      mortonOrdering(fineBoxes);
      Vector<int> fineProcs(fineBoxes.size());
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
      DataIterator ditF = fineGrids.dataIterator();

      // Create a solution Ucoarse with startVal on
      // the coarse grids with a single ghost layer.
      Real startVal = 0.1;
      LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          Ucoarse[ditC].setVal(startVal);
        }

      // Create coarse and fine fluxes that straddle the
      // coarse-fine boundary.
      LevelData<FluxBox> coarseFlux(coarseGrids, ncomp);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& thisBlockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysCoarsePtr =
            coordSysCoarse->getCoordSys(thisBlock);

          FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& coarseFluxFab = coarseFlub[idir];
              const Box& faceBox = coarseFluxFab.box();
              evalScalarFunction(coarseFluxFab, faceBox, idir,
                                 thisCoordSysCoarsePtr, thisBlockBox);
            }
        }

      LevelData<FluxBox> fineFlux(fineGrids, ncomp);
      // Don't just set to norm(fineFlux), because most of it is E+298.
      Real fluxMaxNorm = 0.; // norm(fineFlux, allIntv, MAXNORM);
      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const Box& baseBox = fineGrids[ditF];
          int thisBlock = coordSysFine->whichBlock(baseBox);
          const Box& thisBlockBox = blockBoxesFine[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& fineFluxFab = fineFlub[idir];
              const Box& faceBox = fineFluxFab.box();
              evalScalarFunction(fineFluxFab, faceBox, idir,
                                 thisCoordSysFinePtr, thisBlockBox);
              Real thisMaxNorm = maxnormFaceFab(fineFluxFab, thisBlockBox);
              CH_assert(0 == 0); // dummy for gdb
              fluxMaxNorm = Max(fluxMaxNorm, thisMaxNorm);
            }
        }
      reduceMax(fluxMaxNorm); // Get global max across processors.
      if (verbosity >= 2)
        {
          pout() << " fluxMaxNorm = " << fluxMaxNorm << " maxnorms";
        }

      // Create the flux register.
      MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
      // If SpaceDim = 3, don't refine in the r (last) direction (anisotropic refinement)
      int unrefDir = (SpaceDim == 3) ? SpaceDim-1 : -1;
      bool scaleFlag = true;
      fluxReg.define(fineGrids, coarseGrids, coarseDomain,
                     nRefRatio, ncomp, scaleFlag, unrefDir);
      fluxReg.setToZero();
                  
      // Stick the coarse and fine fluxes into the register
      // to compute the defect.
      Real scale1 = 1.0;
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementCoarse(coarseFlub[idir], scale1, ditC(),
                                      allIntv, allIntv, idir);
            }
        }

      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementFine(fineFlub[idir], scale1, ditF(),
                                    allIntv, allIntv, idir);
            }
        }
      
      // Apply the defect to the solution Ucoarse.
      fluxReg.reflux(Ucoarse, vectorIntv, scale1);

      // Find what the defect and Ucoarse should be.
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];

          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& blockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          const Box& blockFineBox = blockBoxesFine[thisBlock];
          const Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries =
            boundaries[thisBlock];

          FArrayBox defectFab(baseBox, ncomp);
          defectFab.setVal(0.);

          int faceID = 0;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (idir < interpDims)
                    {
                      // Check that Ucoarse is correct on this face.
                      const BlockBoundary& bb = blockBoundaries[faceID];
                      int otherBlock = bb.neighbor();
                      if ( !coarseList.includes(otherBlock) )
                        {
                          // There is fine patch on other side of this face.

                          IndicesTransformation it = bb.getTransformation();

                          int idirOther = bb.dirOther(idir);
                          int reorient = bb.reorientFace(idir);
                          Side::LoHiSide sideOther =
                            (reorient == -1) ? side : flip(side);

                          Box adjBox = adjCellBox(blockBox, idir, side, 1);
                          Box adjOtherBox = it.transform(adjBox);

                          // coarse flux on fullFaceBox
                          Box fullFaceBox(adjBox);
                          fullFaceBox.shiftHalf(idir, -sign(side));

                          // fine flux on fullFaceOtherBox
                          Box fullFaceOtherBox(adjOtherBox);
                          fullFaceOtherBox.shiftHalf(idirOther, sideOther);

                          // Ucoarse adjusted on fullCellBox
                          Box fullCellBox(adjBox);
                          fullCellBox.shift(idir, -sign(side));

                          const FArrayBox& coarseFluxFab = coarseFlub[idir];
                          defectFab.shiftHalf(idir, sign(side));
                          defectFab.plus(coarseFluxFab, fullFaceBox, fullFaceBox, 
                                         -sign(side)*1., 0, 0, ncomp);

                          // Add in sign(side) * reorient * fine flux.

                          // face-centered box
                          Box fullFaceFineBox =
                            refine(fullFaceBox, refRatioVect);
                          FArrayBox fineFluxFab(fullFaceFineBox, ncomp);
                          evalScalarFunction(fineFluxFab, fullFaceFineBox, idir,
                                             thisCoordSysFinePtr, blockFineBox);
                          Box refinementFaceBox(refinementBox);
                          refinementFaceBox.setRange(idir, 0);
                          Real scaling =
                             1.0 / Real(refinementFaceBox.numPts());
                          FArrayBox otherFluxFab(fullFaceBox, ncomp);
                          otherFluxFab.setVal(0.);
                          for (BoxIterator bit(fullFaceBox); bit.ok(); ++bit)
                            {
                              IntVect iv = bit();
                              for (BoxIterator bitRef(refinementFaceBox); bitRef.ok(); ++bitRef)
                                {
                                  IntVect br = bitRef();
                                  for (int comp = 0; comp < ncomp; comp++)
                                    {
                                      otherFluxFab(iv, comp) +=
                                        scaling *
                                        fineFluxFab(refRatioVect*iv + br, comp);
                                    }
                                }
                            }
                          // otherFluxFab *= reorient;
                          defectFab.plus(otherFluxFab, fullFaceBox, fullFaceBox, 
                                         sign(side)*1., 0, 0, ncomp);
                          defectFab.shiftHalf(idir, -sign(side));
                        }
                    }
                  faceID++;
                }
            }
          const FArrayBox& UFab = Ucoarse[ditC];

          FArrayBox UshouldFab(baseBox, ncomp);
          UshouldFab.setVal(startVal);
          UshouldFab += defectFab;

          FArrayBox diffFab(baseBox, ncomp);
          diffFab.copy(UshouldFab);
          diffFab -= UFab;
          Real diffFabMaxNorm = diffFab.norm(baseBox, MAXNORM, 0, ncomp);
          if (verbosity >= 2)
            {
              pout() << " " << diffFabMaxNorm;
            }
          Real relativeMaxNorm = diffFabMaxNorm / fluxMaxNorm;
          Real thisRatio = relativeMaxNorm / precision;
          worstRatio = Max(worstRatio, thisRatio);
          if (thisRatio > 16)
            {
              status--;
            }
          CH_assert(0 == 0); // dummy for gdb
        }
      if (verbosity >= 2)
        {
          pout() << endl;
        }
      CH_assert(0 == 0); // dummy for gdb
    }
  if (verbosity >= 1)
    {
      pout() << "worst relative difference / precision = " << worstRatio << endl;
    }
  
  // Clean up.
  delete coordSysCoarse;
  delete coordSysFine;
  delete cssFactory;

  // Git.
  return status;
}

//-----------------------------------------------------------------------
int
cubedSphereScalarTest()
{
  // Create a cubed sphere coordinate system with six N-by-N blocks.
#if CH_SPACEDIM >= 3
  // Avoid irritating unused variable warning if dimension < 3.
  int nlayers = 8;
#endif
  int N = 32;
  int nRefRatio = 2;

  int interpDims = SpaceDim;
  if (SpaceDim == 3) interpDims = SpaceDim-1;

  //  int layer = 0; // which layer in 3D
  // CubedSphereShellCS* coordSys = new CubedSphereShellCS(N, nlayers);

  IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
  ProblemDomain coarseDomain(IntVect::Zero, domainLen-IntVect::Unit);

  IntVect refRatioVect = IntVect::Unit;
  for (int idir = 0; idir < interpDims; idir++)
    {
      refRatioVect[idir] = nRefRatio;
    }
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatioVect);

  CubedSphereShellCSFactory* cssFactory = new CubedSphereShellCSFactory();
  cssFactory->setHeight(1.);
  cssFactory->setRadius(1.);
  RealVect probLen = RealVect(D_DECL(11*(M_PI/2.), M_PI/2., 1.));
  RealVect dxCoarse = probLen / RealVect(domainLen);
  RealVect dxFine = dxCoarse / RealVect(refRatioVect);
  MultiBlockCoordSys* coordSysCoarse =
    cssFactory->getCoordSys(coarseDomain, dxCoarse);
  MultiBlockCoordSys* coordSysFine =
    cssFactory->getCoordSys(fineDomain, dxFine);

  const Vector<Box>& blockBoxesCoarse =
    coordSysCoarse->mappingBlocks();
  const Vector<Box>& blockBoxesFine =
    coordSysFine->mappingBlocks();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSysCoarse->boundaries();

  int nblocks = blockBoxesCoarse.size();

  Box refinementBox(IntVect::Zero, refRatioVect-IntVect::Unit);
  IntVect baseLen = IntVect(D_DECL(N, N, nlayers));
  Box baseBox(IntVect::Zero, baseLen-IntVect::Unit);

  Vector< List<int> > coarseBlockIndices = getCubedSphereCoarseBlockIndices();

  int ncomp = 1;
  Interval allIntv(0, ncomp-1);
  Interval vectorIntv; // empty

  Real worstRatio = 0.;

  int status = 0;
  for (int icbi = 0; icbi < coarseBlockIndices.size(); icbi++)
    {
      Vector<Box> coarseBoxes;
      Vector<Box> fineBoxes;
      const List<int>& coarseList = coarseBlockIndices[icbi];

      if (verbosity >= 2)
        {
          pout() << "Coarse blocks";
          for (ListIterator<int> li(coarseList); li.ok(); ++li)
            {
              int blockNum = coarseList[li];
              pout() << " " << blockNum;
            }
          pout() << ": maxnorms";
        }

      for (int iblock = 0; iblock < nblocks; iblock++)
        {
          const Box& blockBox = blockBoxesCoarse[iblock];
          const IntVect& blockBoxLo = blockBox.smallEnd();
          if (coarseList.includes(iblock))
            {
              coarseBoxes.push_back(blockBox);
            }
          else
            {
              IntVect fineBlockBoxLo = refRatioVect * blockBoxLo;
              for (BoxIterator bitRef(refinementBox); bitRef.ok(); ++bitRef)
                {
                  IntVect br = bitRef();
                  IntVect thisFineShift = baseLen * br + fineBlockBoxLo;
                  Box thisFineBox(baseBox);
                  thisFineBox.shift(thisFineShift);
                  fineBoxes.push_back(thisFineBox);
                }
            }
        }

      mortonOrdering(coarseBoxes);
      Vector<int> coarseProcs(coarseBoxes.size());
      LoadBalance(coarseProcs, coarseBoxes);
      DisjointBoxLayout coarseGrids(coarseBoxes, coarseProcs, coarseDomain);
      DataIterator ditC = coarseGrids.dataIterator();

      mortonOrdering(fineBoxes);
      Vector<int> fineProcs(fineBoxes.size());
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
      DataIterator ditF = fineGrids.dataIterator();

      // Create a solution Ucoarse with startVal on
      // the coarse grids with a single ghost layer.
      Real startVal = 0.1;
      LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          Ucoarse[ditC].setVal(startVal);
        }

      // Create coarse and fine fluxes that straddle the
      // coarse-fine boundary.
      LevelData<FluxBox> coarseFlux(coarseGrids, ncomp);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& thisBlockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysCoarsePtr =
            coordSysCoarse->getCoordSys(thisBlock);

          FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& coarseFluxFab = coarseFlub[idir];
              const Box& faceBox = coarseFluxFab.box();
              evalScalarFunction(coarseFluxFab, faceBox, idir,
                                 thisCoordSysCoarsePtr, thisBlockBox);
            }
        }

      LevelData<FluxBox> fineFlux(fineGrids, ncomp);
      // Don't just set to norm(fineFlux), because most of it is E+298.
      Real fluxMaxNorm = 0.; // norm(fineFlux, allIntv, MAXNORM);
      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const Box& baseBox = fineGrids[ditF];
          int thisBlock = coordSysFine->whichBlock(baseBox);
          const Box& thisBlockBox = blockBoxesFine[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& fineFluxFab = fineFlub[idir];
              const Box& faceBox = fineFluxFab.box();
              evalScalarFunction(fineFluxFab, faceBox, idir,
                                 thisCoordSysFinePtr, thisBlockBox);
              fluxMaxNorm =
                Max(fluxMaxNorm, maxnormFaceFab(fineFluxFab, thisBlockBox));
            }
        }
      reduceMax(fluxMaxNorm); // Get global max across processors.

      // Create the flux register.
      MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
      // If SpaceDim = 3, don't refine in the r (last) direction (anisotropic refinement)
      int unrefDir = (SpaceDim == 3) ? SpaceDim-1 : -1;
      bool scaleFlag = true;
      fluxReg.define(fineGrids, coarseGrids, coarseDomain,
                     nRefRatio, ncomp, scaleFlag, unrefDir);
      fluxReg.setToZero();
                  
      // Stick the coarse and fine fluxes into the register
      // to compute the defect.
      Real scale1 = 1.0;
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementCoarse(coarseFlub[idir], scale1, ditC(),
                                      allIntv, allIntv, idir);
            }
        }

      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementFineCentered(fineFlub[idir], scale1, ditF(),
                                            vectorIntv, allIntv, allIntv, idir);
            }
        }
      
      // Apply the defect to the solution Ucoarse.
      fluxReg.reflux(Ucoarse, vectorIntv, scale1);

      // Find what the defect and Ucoarse should be.
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];

          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& blockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          const Box& blockFineBox = blockBoxesFine[thisBlock];
          const Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries =
            boundaries[thisBlock];

          FArrayBox defectFab(baseBox, ncomp);
          defectFab.setVal(0.);

          int faceID = 0;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (idir < interpDims)
                    {
                      // Check that Ucoarse is correct on this face.
                      const BlockBoundary& bb = blockBoundaries[faceID];
                      int otherBlock = bb.neighbor();
                      if ( !coarseList.includes(otherBlock) )
                        {
                          // There is fine patch on other side of this face.

                          IndicesTransformation it = bb.getTransformation();

                          int idirOther = bb.dirOther(idir);
                          int reorient = bb.reorientFace(idir);
                          Side::LoHiSide sideOther =
                            (reorient == -1) ? side : flip(side);

                          Box adjBox = adjCellBox(blockBox, idir, side, 1);
                          Box adjOtherBox = it.transform(adjBox);

                          // coarse flux on fullFaceBox
                          Box fullFaceBox(adjBox);
                          fullFaceBox.shiftHalf(idir, -sign(side));

                          // fine flux on fullFaceOtherBox
                          Box fullFaceOtherBox(adjOtherBox);
                          fullFaceOtherBox.shiftHalf(idirOther, sideOther);

                          // Ucoarse adjusted on fullCellBox
                          Box fullCellBox(adjBox);
                          fullCellBox.shift(idir, -sign(side));

                          const FArrayBox& coarseFluxFab = coarseFlub[idir];
                          defectFab.shiftHalf(idir, sign(side));
                          defectFab.plus(coarseFluxFab, fullFaceBox, fullFaceBox, 
                                         -sign(side)*1., 0, 0, ncomp);

                          // Add in sign(side) * reorient * fine flux.

                          // face-centered box
                          Box fullFaceFineBox =
                            refine(fullFaceBox, refRatioVect);
                          FArrayBox fineFluxFab(fullFaceFineBox, ncomp);
                          evalScalarFunction(fineFluxFab, fullFaceFineBox, idir,
                                             thisCoordSysFinePtr, blockFineBox);
                          Box refinementFaceBox(refinementBox);
                          refinementFaceBox.setRange(idir, 0);
                          Real scaling =
                             1.0 / Real(refinementFaceBox.numPts());
                          FArrayBox otherFluxFab(fullFaceBox, ncomp);
                          otherFluxFab.setVal(0.);
                          for (BoxIterator bit(fullFaceBox); bit.ok(); ++bit)
                            {
                              IntVect iv = bit();
                              for (BoxIterator bitRef(refinementFaceBox); bitRef.ok(); ++bitRef)
                                {
                                  IntVect br = bitRef();
                                  for (int comp = 0; comp < ncomp; comp++)
                                    {
                                      otherFluxFab(iv, comp) +=
                                        scaling *
                                        fineFluxFab(refRatioVect*iv + br, comp);
                                    }
                                }
                            }
                          // otherFluxFab *= reorient;
                          defectFab.plus(otherFluxFab, fullFaceBox, fullFaceBox, 
                                         sign(side)*1., 0, 0, ncomp);
                          defectFab.shiftHalf(idir, -sign(side));
                        }
                    }
                  faceID++;
                }
            }
          const FArrayBox& UFab = Ucoarse[ditC];

          FArrayBox UshouldFab(baseBox, ncomp);
          UshouldFab.setVal(startVal);
          UshouldFab += defectFab;

          FArrayBox diffFab(baseBox, ncomp);
          diffFab.copy(UshouldFab);
          diffFab -= UFab;
          Real diffFabMaxNorm = diffFab.norm(baseBox, MAXNORM, 0, ncomp);
          if (verbosity >= 2)
            {
              pout() << " " << diffFabMaxNorm;
            }
          Real relativeMaxNorm = diffFabMaxNorm / fluxMaxNorm;
          Real thisRatio = relativeMaxNorm / precision;
          worstRatio = Max(worstRatio, thisRatio);
          if (thisRatio > 16)
            {
              status--;
            }
          CH_assert(0 == 0); // dummy for gdb
        }
      if (verbosity >= 2)
        {
          pout() << endl;
        }
      CH_assert(0 == 0); // dummy for gdb
    }
  if (verbosity >= 2)
    {
      pout() << "worst relative difference / precision = " << worstRatio << endl;
    }
  
  // Clean up.
  delete coordSysCoarse;
  delete coordSysFine;
  delete cssFactory;

  // Git.
  return status;
}

//-----------------------------------------------------------------------

int
cubedSphereVectorTest()
{
  // Create a cubed sphere coordinate system with six N-by-N blocks.
#if CH_SPACEDIM >= 3
  // Avoid irritating unused variable warning if dimension < 3.
  int nlayers = 8;
#endif
  int N = 16;
  int nRefRatio = 2;

  int interpDims = SpaceDim;
  if (SpaceDim == 3) interpDims = SpaceDim-1;

  //  int layer = 0; // which layer in 3D
  // CubedSphereShellCS* coordSys = new CubedSphereShellCS(N, nlayers);

  IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
  ProblemDomain coarseDomain(IntVect::Zero, domainLen-IntVect::Unit);

  IntVect refRatioVect = IntVect::Unit;
  for (int idir = 0; idir < interpDims; idir++)
    {
      refRatioVect[idir] = nRefRatio;
    }
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatioVect);

  CubedSphereShellCSFactory* cssFactory = new CubedSphereShellCSFactory();
  cssFactory->setHeight(1.);
  cssFactory->setRadius(1.);
  RealVect probLen = RealVect(D_DECL(11*(M_PI/2.), M_PI/2., 1.));
  RealVect dxCoarse = probLen / RealVect(domainLen);
  RealVect dxFine = dxCoarse / RealVect(refRatioVect);
  MultiBlockCoordSys* coordSysCoarse =
    cssFactory->getCoordSys(coarseDomain, dxCoarse);
  MultiBlockCoordSys* coordSysFine =
    cssFactory->getCoordSys(fineDomain, dxFine);

  const Vector<Box>& blockBoxesCoarse =
    coordSysCoarse->mappingBlocks();
  const Vector<Box>& blockBoxesFine =
    coordSysFine->mappingBlocks();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSysCoarse->boundaries();

  int nblocks = blockBoxesCoarse.size();

  Box refinementBox(IntVect::Zero, refRatioVect-IntVect::Unit);
  IntVect baseLen = IntVect(D_DECL(N, N, nlayers));
  Box baseBox(IntVect::Zero, baseLen-IntVect::Unit);

  Vector< List<int> > coarseBlockIndices = getCubedSphereCoarseBlockIndices();

  int ncomp = SpaceDim;
  Interval allIntv(0, ncomp-1);
  Interval vectorIntv = allIntv;

  Real worstRatio = 0.;

  int status = 0;
  for (int icbi = 0; icbi < coarseBlockIndices.size(); icbi++)
    {
      Vector<Box> coarseBoxes;
      Vector<Box> fineBoxes;
      const List<int>& coarseList = coarseBlockIndices[icbi];

      if (verbosity >= 1)
        {
          pout() << "Coarse blocks";
          for (ListIterator<int> li(coarseList); li.ok(); ++li)
            {
              int blockNum = coarseList[li];
              pout() << " " << blockNum;
            }
          if (verbosity == 1)
            pout() << ": maxnorms";
          else
            pout() << endl;
          CH_assert(0 == 0); // dummy for gdb
        }

      for (int iblock = 0; iblock < nblocks; iblock++)
        {
          const Box& blockBox = blockBoxesCoarse[iblock];
          const IntVect& blockBoxLo = blockBox.smallEnd();
          if (coarseList.includes(iblock))
            {
              coarseBoxes.push_back(blockBox);
            }
          else
            {
              IntVect fineBlockBoxLo = refRatioVect * blockBoxLo;
              for (BoxIterator bitRef(refinementBox); bitRef.ok(); ++bitRef)
                {
                  IntVect br = bitRef();
                  IntVect thisFineShift = baseLen * br + fineBlockBoxLo;
                  Box thisFineBox(baseBox);
                  thisFineBox.shift(thisFineShift);
                  fineBoxes.push_back(thisFineBox);
                }
            }
        }

      mortonOrdering(coarseBoxes);
      Vector<int> coarseProcs(coarseBoxes.size());
      LoadBalance(coarseProcs, coarseBoxes);
      DisjointBoxLayout coarseGrids(coarseBoxes, coarseProcs, coarseDomain);
      DataIterator ditC = coarseGrids.dataIterator();

      mortonOrdering(fineBoxes);
      Vector<int> fineProcs(fineBoxes.size());
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
      DataIterator ditF = fineGrids.dataIterator();

      // Create a solution Ucoarse with startVal on
      // the coarse grids with a single ghost layer.
      Real startVal = 0.1;
      LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          Ucoarse[ditC].setVal(startVal);
        }

      // Create coarse and fine fluxes that straddle the
      // coarse-fine boundary.
      LevelData<FluxBox> coarseFlux(coarseGrids, ncomp);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& blockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysCoarsePtr =
            coordSysCoarse->getCoordSys(thisBlock);

          FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& coarseFluxFab = coarseFlub[idir];
              const Box& faceBox = coarseFluxFab.box();
              evalVectorFunction(coarseFluxFab, faceBox, idir,
                                 thisCoordSysCoarsePtr, blockBox);
            }
        }

      LevelData<FluxBox> fineFlux(fineGrids, ncomp);
      // Don't just set to norm(fineFlux), because most of it is E+298.
      Real fluxMaxNorm = 0.; // norm(fineFlux, allIntv, MAXNORM);
      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const Box& baseBox = fineGrids[ditF];
          Box coarseBaseBox = coarsen(baseBox, refRatioVect);
          int thisBlock = coordSysCoarse->whichBlock(coarseBaseBox);
          const Box& blockBox = blockBoxesFine[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& fineFluxFab = fineFlub[idir];
              const Box& faceBox = fineFluxFab.box();
              evalVectorFunction(fineFluxFab, faceBox, idir,
                                 thisCoordSysFinePtr, blockBox);
              fluxMaxNorm =
                Max(fluxMaxNorm, maxnormFaceFab(fineFluxFab, blockBox));
            }
        }
      reduceMax(fluxMaxNorm); // Get global max across processors.

      // Create the flux register.
      MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
      // If SpaceDim = 3, don't refine in the r (last) direction (anisotropic refinement)
      int unrefDir = (SpaceDim == 3) ? SpaceDim-1 : -1;
      bool scaleFlag = true;
      fluxReg.define(fineGrids, coarseGrids, coarseDomain,
                     nRefRatio, ncomp, scaleFlag, unrefDir);
      fluxReg.setToZero();
                  
      // Stick the coarse and fine fluxes into the register
      // to compute the defect.
      Real scale1 = 1.0;
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementCoarse(coarseFlub[idir], scale1, ditC(),
                                      allIntv, allIntv, idir);
            }
        }

      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementFineCentered(fineFlub[idir], scale1, ditF(),
                                            vectorIntv, allIntv, allIntv, idir);
            }
        }
      
      // Apply the defect to the solution Ucoarse.
      fluxReg.reflux(Ucoarse, vectorIntv, scale1);

      // Find what the defect and Ucoarse should be.
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];

          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& blockBox = blockBoxesCoarse[thisBlock];
          const Box& blockFineBox = blockBoxesFine[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          const Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries =
            boundaries[thisBlock];

          FArrayBox defectFab(baseBox, ncomp);
          defectFab.setVal(0.);

          int faceID = 0;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (idir < interpDims)
                    {
                      // Check that Ucoarse is correct on this face.
                      const BlockBoundary& bb = blockBoundaries[faceID];
                      int otherBlock = bb.neighbor();
                      if ( !coarseList.includes(otherBlock) )
                        {
                          // There is fine patch on other side of this face.

                          IndicesTransformation it = bb.getTransformation();

                          int idirOther = bb.dirOther(idir);
                          int reorient = bb.reorientFace(idir);
                          Side::LoHiSide sideOther =
                            (reorient == -1) ? side : flip(side);

                          Box adjBox = adjCellBox(blockBox, idir, side, 1);
                          Box adjOtherBox = it.transform(adjBox);

                          // coarse flux on fullFaceBox
                          Box fullFaceBox(adjBox);
                          fullFaceBox.shiftHalf(idir, -sign(side));

                          // fine flux on fullFaceOtherBox
                          Box fullFaceOtherBox(adjOtherBox);
                          fullFaceOtherBox.shiftHalf(idirOther, sideOther);

                          // Ucoarse adjusted on fullCellBox
                          Box fullCellBox(adjBox);
                          fullCellBox.shift(idir, -sign(side));

                          const FArrayBox& coarseFluxFab = coarseFlub[idir];
                          defectFab.shiftHalf(idir, sign(side));
                          defectFab.plus(coarseFluxFab, fullFaceBox, fullFaceBox, 
                                         -sign(side)*1., 0, 0, ncomp);

                          // Add in sign(side) * reorient * fine flux.

                          // face-centered box
                          Box fullFaceFineBox =
                            refine(fullFaceBox, refRatioVect);
                          FArrayBox fineFluxFab(fullFaceFineBox, ncomp);
                          evalVectorFunction(fineFluxFab, fullFaceFineBox, idir,
                                             thisCoordSysFinePtr, blockFineBox);
                          Box refinementFaceBox(refinementBox);
                          refinementFaceBox.setRange(idir, 0);
                          Real scaling =
                             1.0 / Real(refinementFaceBox.numPts());
                          FArrayBox otherFluxFab(fullFaceBox, ncomp);
                          otherFluxFab.setVal(0.);
                          for (BoxIterator bit(fullFaceBox); bit.ok(); ++bit)
                            {
                              IntVect iv = bit();
                              for (BoxIterator bitRef(refinementFaceBox); bitRef.ok(); ++bitRef)
                                {
                                  IntVect br = bitRef();
                                  for (int comp = 0; comp < ncomp; comp++)
                                    {
                                      otherFluxFab(iv, comp) +=
                                        scaling *
                                        fineFluxFab(refRatioVect*iv + br, comp);
                                    }
                                }
                            }
                          // otherFluxFab *= reorient;
                          defectFab.plus(otherFluxFab, fullFaceBox, fullFaceBox, 
                                         sign(side)*1., 0, 0, ncomp);
                          defectFab.shiftHalf(idir, -sign(side));
                        }
                    }
                  faceID++;
                }
            }
          const FArrayBox& UFab = Ucoarse[ditC];

          FArrayBox UshouldFab(baseBox, ncomp);
          UshouldFab.setVal(startVal);
          UshouldFab += defectFab;

          FArrayBox diffFab(baseBox, ncomp);
          diffFab.copy(UshouldFab);
          diffFab -= UFab;
          RealVect diffFabMaxNorm;
          Real diffFabMaxNormAll = 0.;
          for (int comp = 0; comp < SpaceDim; comp++)
            {
              diffFabMaxNorm[comp] = diffFab.norm(baseBox, MAXNORM, comp, 1);
              if (diffFabMaxNorm[comp] > diffFabMaxNormAll)
                diffFabMaxNormAll = diffFabMaxNorm[comp];
            }
          if (verbosity >= 1)
            {
              if (verbosity >= 2) pout() << "block=" << thisBlock
                                         << " maxnorm";
              pout() << " " << diffFabMaxNorm;
            }
          if (verbosity >= 3)
            {
              // Look at difference on each face.
              int faceIDdebug = 0;
              for (SideIterator sitdebug; sitdebug.ok(); ++sitdebug)
                {
                  Side::LoHiSide sidedebug = sitdebug();
                  for (int idirdebug = 0; idirdebug < SpaceDim; idirdebug++)
                    {
                      Box adjWithinBox =
                        adjCellBox(baseBox, idirdebug, sidedebug, -1);
                      // Omit corners.
                      adjWithinBox.grow(BASISV(idirdebug) - IntVect::Unit);
                      RealVect diffFabFaceMaxNorm;
                      for (int comp = 0; comp < SpaceDim; comp++)
                        {
                          diffFabFaceMaxNorm[comp] =
                            diffFab.norm(adjWithinBox, 0, comp, 1);
                        }
                      char letts[] = {'x', 'y', 'z'};
                      pout() << " " << letts[idirdebug]
                             << ((sidedebug == Side::Lo) ? "Lo" : "Hi")
                             << "=" << diffFabFaceMaxNorm;
                      faceIDdebug++;
                    }
                }
            }
          if (verbosity >= 3) pout() << endl;

          Real relativeMaxNorm = diffFabMaxNormAll / fluxMaxNorm;
          Real thisRatio = relativeMaxNorm / precision;
          worstRatio = Max(worstRatio, thisRatio);
          if (thisRatio > 16)
            {
              status--;
            }
          CH_assert(0 == 0); // dummy for gdb
        }
      if (verbosity == 1)
        {
          pout() << endl;
        }
      CH_assert(0 == 0); // dummy for gdb
    }
  if (verbosity >= 1)
    {
      pout() << "worst relative difference / precision = " << worstRatio << endl;
    }
  
  // Clean up.
  delete coordSysCoarse;
  delete coordSysFine;
  delete cssFactory;

  // Git.
  return status;
}

//-----------------------------------------------------------------------

int
cubedSphereVectorScalarTest()
{
  // Create a cubed sphere coordinate system with six N-by-N blocks.
#if CH_SPACEDIM >= 3
  // Avoid irritating unused variable warning if dimension < 3.
  int nlayers = 8;
#endif
  int N = 16;
  int nRefRatio = 2;

  int interpDims = SpaceDim;
  if (SpaceDim == 3) interpDims = SpaceDim-1;

  //  int layer = 0; // which layer in 3D
  // CubedSphereShellCS* coordSys = new CubedSphereShellCS(N, nlayers);

  IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
  ProblemDomain coarseDomain(IntVect::Zero, domainLen-IntVect::Unit);

  IntVect refRatioVect = IntVect::Unit;
  for (int idir = 0; idir < interpDims; idir++)
    {
      refRatioVect[idir] = nRefRatio;
    }
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatioVect);

  CubedSphereShellCSFactory* cssFactory = new CubedSphereShellCSFactory();
  cssFactory->setHeight(1.);
  cssFactory->setRadius(1.);
  RealVect probLen = RealVect(D_DECL(11*(M_PI/2.), M_PI/2., 1.));
  RealVect dxCoarse = probLen / RealVect(domainLen);
  RealVect dxFine = dxCoarse / RealVect(refRatioVect);
  MultiBlockCoordSys* coordSysCoarse =
    cssFactory->getCoordSys(coarseDomain, dxCoarse);
  MultiBlockCoordSys* coordSysFine =
    cssFactory->getCoordSys(fineDomain, dxFine);

  const Vector<Box>& blockBoxesCoarse =
    coordSysCoarse->mappingBlocks();
  const Vector<Box>& blockBoxesFine =
    coordSysFine->mappingBlocks();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSysCoarse->boundaries();

  int nblocks = blockBoxesCoarse.size();

  Box refinementBox(IntVect::Zero, refRatioVect-IntVect::Unit);
  IntVect baseLen = IntVect(D_DECL(N, N, nlayers));
  Box baseBox(IntVect::Zero, baseLen-IntVect::Unit);

  Vector< List<int> > coarseBlockIndices = getCubedSphereCoarseBlockIndices();

  Interval vectorIntv(0, SpaceDim-1);
  int iscalar = SpaceDim;
  Interval scalarIntv(iscalar, iscalar);
  int nscalar = 1;
  int ncomp = SpaceDim + nscalar;
  Interval allIntv(0, ncomp-1);

  Real worstRatio = 0.;

  int status = 0;
  for (int icbi = 0; icbi < coarseBlockIndices.size(); icbi++)
    {
      Vector<Box> coarseBoxes;
      Vector<Box> fineBoxes;
      const List<int>& coarseList = coarseBlockIndices[icbi];

      if (verbosity >= 1)
        {
          pout() << "Coarse blocks";
          for (ListIterator<int> li(coarseList); li.ok(); ++li)
            {
              int blockNum = coarseList[li];
              pout() << " " << blockNum;
            }
          if (verbosity == 1)
            pout() << ": maxnorms";
          else
            pout() << endl;
          CH_assert(0 == 0); // dummy for gdb
        }

      for (int iblock = 0; iblock < nblocks; iblock++)
        {
          const Box& blockBox = blockBoxesCoarse[iblock];
          const IntVect& blockBoxLo = blockBox.smallEnd();
          if (coarseList.includes(iblock))
            {
              coarseBoxes.push_back(blockBox);
            }
          else
            {
              IntVect fineBlockBoxLo = refRatioVect * blockBoxLo;
              for (BoxIterator bitRef(refinementBox); bitRef.ok(); ++bitRef)
                {
                  IntVect br = bitRef();
                  IntVect thisFineShift = baseLen * br + fineBlockBoxLo;
                  Box thisFineBox(baseBox);
                  thisFineBox.shift(thisFineShift);
                  fineBoxes.push_back(thisFineBox);
                }
            }
        }

      mortonOrdering(coarseBoxes);
      Vector<int> coarseProcs(coarseBoxes.size());
      LoadBalance(coarseProcs, coarseBoxes);
      DisjointBoxLayout coarseGrids(coarseBoxes, coarseProcs, coarseDomain);
      DataIterator ditC = coarseGrids.dataIterator();

      mortonOrdering(fineBoxes);
      Vector<int> fineProcs(fineBoxes.size());
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
      DataIterator ditF = fineGrids.dataIterator();

      // Create a solution Ucoarse with startVal on
      // the coarse grids with a single ghost layer.
      Real startVal = 0.1;
      LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          Ucoarse[ditC].setVal(startVal);
        }

      // Create coarse and fine fluxes that straddle the
      // coarse-fine boundary.
      LevelData<FluxBox> coarseFlux(coarseGrids, ncomp);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& blockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysCoarsePtr =
            coordSysCoarse->getCoordSys(thisBlock);

          FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& coarseFluxFab = coarseFlub[idir];
              const Box& faceBox = coarseFluxFab.box();
              FArrayBox coarseVectorFluxFab(vectorIntv, coarseFluxFab);
              evalVectorFunction(coarseVectorFluxFab, faceBox, idir,
                                 thisCoordSysCoarsePtr, blockBox);
              FArrayBox coarseScalarFluxFab(scalarIntv, coarseFluxFab);
              evalScalarFunction(coarseScalarFluxFab, faceBox, idir,
                                 thisCoordSysCoarsePtr, blockBox);
            }
        }

      LevelData<FluxBox> fineFlux(fineGrids, ncomp);
      // Don't just set to norm(fineFlux), because most of it is E+298.
      Real fluxMaxNorm = 0.; // norm(fineFlux, allIntv, MAXNORM);
      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const Box& baseBox = fineGrids[ditF];
          Box coarseBaseBox = coarsen(baseBox, refRatioVect);
          int thisBlock = coordSysCoarse->whichBlock(coarseBaseBox);
          const Box& blockBox = blockBoxesFine[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& fineFluxFab = fineFlub[idir];
              const Box& faceBox = fineFluxFab.box();
              FArrayBox fineVectorFluxFab(vectorIntv, fineFluxFab);
              evalVectorFunction(fineVectorFluxFab, faceBox, idir,
                                 thisCoordSysFinePtr, blockBox);
              FArrayBox fineScalarFluxFab(scalarIntv, fineFluxFab);
              evalScalarFunction(fineScalarFluxFab, faceBox, idir,
                                 thisCoordSysFinePtr, blockBox);
              fluxMaxNorm =
                Max(fluxMaxNorm, maxnormFaceFab(fineFluxFab, blockBox));
            }
        }
      reduceMax(fluxMaxNorm); // Get global max across processors.

      // Create the flux register.
      MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
      // If SpaceDim = 3, don't refine in the r (last) direction (anisotropic refinement)
      int unrefDir = (SpaceDim == 3) ? SpaceDim-1 : -1;
      bool scaleFlag = true;
      fluxReg.define(fineGrids, coarseGrids, coarseDomain,
                     nRefRatio, ncomp, scaleFlag, unrefDir);
      fluxReg.setToZero();
                  
      // Stick the coarse and fine fluxes into the register
      // to compute the defect.
      Real scale1 = 1.0;
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementCoarse(coarseFlub[idir], scale1, ditC(),
                                      allIntv, allIntv, idir);
            }
        }

      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementFineCentered(fineFlub[idir], scale1, ditF(),
                                            vectorIntv, allIntv, allIntv, idir);
                                           
            }
        }
      
      // Apply the defect to the solution Ucoarse.
      fluxReg.reflux(Ucoarse, vectorIntv, scale1);

      // Find what the defect and Ucoarse should be.
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];

          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& blockBox = blockBoxesCoarse[thisBlock];
          const Box& blockFineBox = blockBoxesFine[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          const Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries =
            boundaries[thisBlock];

          FArrayBox defectFab(baseBox, ncomp);
          defectFab.setVal(0.);

          int faceID = 0;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (idir < interpDims)
                    {
                      // Check that Ucoarse is correct on this face.
                      const BlockBoundary& bb = blockBoundaries[faceID];
                      int otherBlock = bb.neighbor();
                      if ( !coarseList.includes(otherBlock) )
                        {
                          // There is fine patch on other side of this face.

                          IndicesTransformation it = bb.getTransformation();

                          int idirOther = bb.dirOther(idir);
                          int reorient = bb.reorientFace(idir);
                          Side::LoHiSide sideOther =
                            (reorient == -1) ? side : flip(side);

                          Box adjBox = adjCellBox(blockBox, idir, side, 1);
                          Box adjOtherBox = it.transform(adjBox);

                          // coarse flux on fullFaceBox
                          Box fullFaceBox(adjBox);
                          fullFaceBox.shiftHalf(idir, -sign(side));

                          // fine flux on fullFaceOtherBox
                          Box fullFaceOtherBox(adjOtherBox);
                          fullFaceOtherBox.shiftHalf(idirOther, sideOther);

                          // Ucoarse adjusted on fullCellBox
                          Box fullCellBox(adjBox);
                          fullCellBox.shift(idir, -sign(side));

                          const FArrayBox& coarseFluxFab = coarseFlub[idir];
                          defectFab.shiftHalf(idir, sign(side));
                          defectFab.plus(coarseFluxFab, fullFaceBox, fullFaceBox, 
                                         -sign(side)*1., 0, 0, ncomp);

                          // Add in sign(side) * reorient * fine flux.

                          // face-centered box
                          Box fullFaceFineBox =
                            refine(fullFaceBox, refRatioVect);
                          FArrayBox fineFluxFab(fullFaceFineBox, ncomp);
                          FArrayBox fineVectorFluxFab(vectorIntv, fineFluxFab);
                          evalVectorFunction(fineVectorFluxFab, fullFaceFineBox,
                                             idir,
                                             thisCoordSysFinePtr, blockFineBox);
                          FArrayBox fineScalarFluxFab(scalarIntv, fineFluxFab);
                          evalScalarFunction(fineScalarFluxFab, fullFaceFineBox,
                                             idir,
                                             thisCoordSysFinePtr, blockFineBox);
                          Box refinementFaceBox(refinementBox);
                          refinementFaceBox.setRange(idir, 0);
                          Real scaling =
                             1.0 / Real(refinementFaceBox.numPts());
                          FArrayBox otherFluxFab(fullFaceBox, ncomp);
                          otherFluxFab.setVal(0.);
                          for (BoxIterator bit(fullFaceBox); bit.ok(); ++bit)
                            {
                              IntVect iv = bit();
                              for (BoxIterator bitRef(refinementFaceBox); bitRef.ok(); ++bitRef)
                                {
                                  IntVect br = bitRef();
                                  for (int comp = 0; comp < ncomp; comp++)
                                    {
                                      otherFluxFab(iv, comp) +=
                                        scaling *
                                        fineFluxFab(refRatioVect*iv + br, comp);
                                    }
                                }
                            }
                          // otherFluxFab *= reorient;
                          defectFab.plus(otherFluxFab, fullFaceBox, fullFaceBox, 
                                         sign(side)*1., 0, 0, ncomp);
                          defectFab.shiftHalf(idir, -sign(side));
                        }
                    }
                  faceID++;
                }
            }
          const FArrayBox& UFab = Ucoarse[ditC];

          FArrayBox UshouldFab(baseBox, ncomp);
          UshouldFab.setVal(startVal);
          UshouldFab += defectFab;

          FArrayBox diffFab(baseBox, ncomp);
          diffFab.copy(UshouldFab);
          diffFab -= UFab;
          Vector<Real> diffFabMaxNorm(ncomp);
          Real diffFabMaxNormAll = 0.;
          for (int comp = 0; comp < SpaceDim; comp++)
            {
              diffFabMaxNorm[comp] = diffFab.norm(baseBox, MAXNORM, comp, 1);
              if (diffFabMaxNorm[comp] > diffFabMaxNormAll)
                diffFabMaxNormAll = diffFabMaxNorm[comp];
            }
          if (verbosity >= 1)
            {
              if (verbosity >= 3) pout() << "block=" << thisBlock
                                         << " maxnorm";
              pout() << " (" << diffFabMaxNorm << ")";
            }
          if (verbosity >= 3)
            {
              // Look at difference on each face.
              int faceIDdebug = 0;
              for (SideIterator sitdebug; sitdebug.ok(); ++sitdebug)
                {
                  Side::LoHiSide sidedebug = sitdebug();
                  for (int idirdebug = 0; idirdebug < SpaceDim; idirdebug++)
                    {
                      Box adjWithinBox =
                        adjCellBox(baseBox, idirdebug, sidedebug, -1);
                      // Omit corners.
                      adjWithinBox.grow(BASISV(idirdebug) - IntVect::Unit);
                      Vector<Real> diffFabFaceMaxNorm(ncomp);
                      for (int comp = 0; comp < ncomp; comp++)
                        {
                          diffFabFaceMaxNorm[comp] =
                            diffFab.norm(adjWithinBox, 0, comp, 1);
                        }
                      char letts[] = {'x', 'y', 'z'};
                      pout() << " " << letts[idirdebug]
                             << ((sidedebug == Side::Lo) ? "Lo" : "Hi")
                             << "=(" << diffFabFaceMaxNorm << ")";
                      faceIDdebug++;
                    }
                }
            }
          if (verbosity >= 3) pout() << endl;

          Real relativeMaxNorm = diffFabMaxNormAll / fluxMaxNorm;
          Real thisRatio = relativeMaxNorm / precision;
          worstRatio = Max(worstRatio, thisRatio);
          if (thisRatio > 16)
            {
              status--;
            }
          CH_assert(0 == 0); // dummy for gdb
        }
      if (verbosity == 1)
        {
          pout() << endl;
        }
      CH_assert(0 == 0); // dummy for gdb
    }
  if (verbosity >= 1)
    {
      pout() << "worst relative difference / precision = " << worstRatio << endl;
    }
  
  // Clean up.
  delete coordSysCoarse;
  delete coordSysFine;
  delete cssFactory;

  // Git.
  return status;
}

//-----------------------------------------------------------------------
int
cubedSphereScalarAvgTest()
{
  // Create a cubed sphere coordinate system with six N-by-N blocks.
#if CH_SPACEDIM >= 3
  // Avoid irritating unused variable warning if dimension < 3.
  int nlayers = 8;
#endif
  int N = 32;
  int nRefRatio = 2;

  int interpDims = SpaceDim;
  if (SpaceDim == 3) interpDims = SpaceDim-1;

  //  int layer = 0; // which layer in 3D
  // CubedSphereShellCS* coordSys = new CubedSphereShellCS(N, nlayers);

  IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
  ProblemDomain coarseDomain(IntVect::Zero, domainLen-IntVect::Unit);

  IntVect refRatioVect = IntVect::Unit;
  for (int idir = 0; idir < interpDims; idir++)
    {
      refRatioVect[idir] = nRefRatio;
    }
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refRatioVect);

  CubedSphereShellCSFactory* cssFactory = new CubedSphereShellCSFactory();
  cssFactory->setHeight(1.);
  cssFactory->setRadius(1.);
  RealVect probLen = RealVect(D_DECL(11*(M_PI/2.), M_PI/2., 1.));
  RealVect dxCoarse = probLen / RealVect(domainLen);
  RealVect dxFine = dxCoarse / RealVect(refRatioVect);
  MultiBlockCoordSys* coordSysCoarse =
    cssFactory->getCoordSys(coarseDomain, dxCoarse);
  MultiBlockCoordSys* coordSysFine =
    cssFactory->getCoordSys(fineDomain, dxFine);

  const Vector<Box>& blockBoxesCoarse =
    coordSysCoarse->mappingBlocks();
  const Vector<Box>& blockBoxesFine =
    coordSysFine->mappingBlocks();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSysCoarse->boundaries();

  int nblocks = blockBoxesCoarse.size();

  Box refinementBox(IntVect::Zero, refRatioVect-IntVect::Unit);
  IntVect baseLen = IntVect(D_DECL(N, N, nlayers));
  Box baseBox(IntVect::Zero, baseLen-IntVect::Unit);

  Vector< List<int> > coarseBlockIndices = getCubedSphereCoarseBlockIndices();

  int ncomp = 1;
  Interval allIntv(0, ncomp-1);
  Interval vectorIntv; // empty

  Real worstRatio = 0.;

  int status = 0;
  for (int icbi = 0; icbi < coarseBlockIndices.size(); icbi++)
    {
      Vector<Box> coarseBoxes;
      Vector<Box> fineBoxes;
      const List<int>& coarseList = coarseBlockIndices[icbi];

      if (verbosity >= 2)
        {
          pout() << "Coarse blocks";
          for (ListIterator<int> li(coarseList); li.ok(); ++li)
            {
              int blockNum = coarseList[li];
              pout() << " " << blockNum;
            }
          pout() << ": maxnorms";
        }

      for (int iblock = 0; iblock < nblocks; iblock++)
        {
          const Box& blockBox = blockBoxesCoarse[iblock];
          const IntVect& blockBoxLo = blockBox.smallEnd();
          if (coarseList.includes(iblock))
            {
              coarseBoxes.push_back(blockBox);
            }
          else
            {
              IntVect fineBlockBoxLo = refRatioVect * blockBoxLo;
              for (BoxIterator bitRef(refinementBox); bitRef.ok(); ++bitRef)
                {
                  IntVect br = bitRef();
                  IntVect thisFineShift = baseLen * br + fineBlockBoxLo;
                  Box thisFineBox(baseBox);
                  thisFineBox.shift(thisFineShift);
                  fineBoxes.push_back(thisFineBox);
                }
            }
        }

      mortonOrdering(coarseBoxes);
      Vector<int> coarseProcs(coarseBoxes.size());
      LoadBalance(coarseProcs, coarseBoxes);
      DisjointBoxLayout coarseGrids(coarseBoxes, coarseProcs, coarseDomain);
      DataIterator ditC = coarseGrids.dataIterator();

      mortonOrdering(fineBoxes);
      Vector<int> fineProcs(fineBoxes.size());
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
      DataIterator ditF = fineGrids.dataIterator();

      // Create a solution Ucoarse with startVal on
      // the coarse grids with a single ghost layer.
      Real startVal = 0.1;
      LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          Ucoarse[ditC].setVal(startVal);
        }

      // Create coarse and fine fluxes that straddle the
      // coarse-fine boundary.
      LevelData<FluxBox> coarseFlux(coarseGrids, ncomp);
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& thisBlockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysCoarsePtr =
            coordSysCoarse->getCoordSys(thisBlock);

          FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& coarseFluxFab = coarseFlub[idir];
              const Box& faceBox = coarseFluxFab.box();
              evalScalarFunction(coarseFluxFab, faceBox, idir,
                                 thisCoordSysCoarsePtr, thisBlockBox);
            }
        }

      LevelData<FluxBox> fineFlux(fineGrids, ncomp);
      // Don't just set to norm(fineFlux), because most of it is E+298.
      Real fluxMaxNorm = 0.; // norm(fineFlux, allIntv, MAXNORM);
      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const Box& baseBox = fineGrids[ditF];
          int thisBlock = coordSysFine->whichBlock(baseBox);
          const Box& thisBlockBox = blockBoxesFine[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox& fineFluxFab = fineFlub[idir];
              const Box& faceBox = fineFluxFab.box();
              evalScalarFunction(fineFluxFab, faceBox, idir,
                                 thisCoordSysFinePtr, thisBlockBox);
              fluxMaxNorm =
                Max(fluxMaxNorm, maxnormFaceFab(fineFluxFab, thisBlockBox));
            }
        }
      reduceMax(fluxMaxNorm); // Get global max across processors.

      // Create the flux register.
      MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
      // If SpaceDim = 3, don't refine in the r (last) direction (anisotropic refinement)
      int unrefDir = (SpaceDim == 3) ? SpaceDim-1 : -1;
      bool scaleFlag = true;
      fluxReg.define(fineGrids, coarseGrids, coarseDomain,
                     nRefRatio, ncomp, scaleFlag, unrefDir);
      fluxReg.setToZero();
                  
      // Stick the coarse and fine fluxes into the register
      // to compute the defect.
      Real scale1 = 1.0;
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementCoarse(coarseFlub[idir], scale1, ditC(),
                                      allIntv, allIntv, idir);
            }
        }

      for (ditF.begin(); ditF.ok(); ++ditF)
        {
          const FluxBox& fineFlub = fineFlux[ditF];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxReg.incrementFineAveraged(fineFlub[idir], fineFlub[idir],
                                            scale1, ditF(),
                                            vectorIntv, allIntv, allIntv, idir);
            }
        }
      
      // Apply the defect to the solution Ucoarse.
      fluxReg.reflux(Ucoarse, vectorIntv, scale1);

      // Find what the defect and Ucoarse should be.
      for (ditC.begin(); ditC.ok(); ++ditC)
        {
          const FluxBox& coarseFlub = coarseFlux[ditC];

          const Box& baseBox = coarseGrids[ditC];
          int thisBlock = coordSysCoarse->whichBlock(baseBox);
          const Box& blockBox = blockBoxesCoarse[thisBlock];
          const NewCoordSys* thisCoordSysFinePtr =
            coordSysFine->getCoordSys(thisBlock);
          const Box& blockFineBox = blockBoxesFine[thisBlock];
          const Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries =
            boundaries[thisBlock];

          FArrayBox defectFab(baseBox, ncomp);
          defectFab.setVal(0.);

          int faceID = 0;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (idir < interpDims)
                    {
                      // Check that Ucoarse is correct on this face.
                      const BlockBoundary& bb = blockBoundaries[faceID];
                      int otherBlock = bb.neighbor();
                      if ( !coarseList.includes(otherBlock) )
                        {
                          // There is fine patch on other side of this face.

                          IndicesTransformation it = bb.getTransformation();

                          int idirOther = bb.dirOther(idir);
                          int reorient = bb.reorientFace(idir);
                          Side::LoHiSide sideOther =
                            (reorient == -1) ? side : flip(side);

                          Box adjBox = adjCellBox(blockBox, idir, side, 1);
                          Box adjOtherBox = it.transform(adjBox);

                          // coarse flux on fullFaceBox
                          Box fullFaceBox(adjBox);
                          fullFaceBox.shiftHalf(idir, -sign(side));

                          // fine flux on fullFaceOtherBox
                          Box fullFaceOtherBox(adjOtherBox);
                          fullFaceOtherBox.shiftHalf(idirOther, sideOther);

                          // Ucoarse adjusted on fullCellBox
                          Box fullCellBox(adjBox);
                          fullCellBox.shift(idir, -sign(side));

                          const FArrayBox& coarseFluxFab = coarseFlub[idir];
                          defectFab.shiftHalf(idir, sign(side));
                          defectFab.plus(coarseFluxFab, fullFaceBox, fullFaceBox, 
                                         -sign(side)*1., 0, 0, ncomp);

                          // Add in sign(side) * reorient * fine flux.

                          // face-centered box
                          Box fullFaceFineBox =
                            refine(fullFaceBox, refRatioVect);
                          FArrayBox fineFluxFab(fullFaceFineBox, ncomp);
                          evalScalarFunction(fineFluxFab, fullFaceFineBox, idir,
                                             thisCoordSysFinePtr, blockFineBox);
                          Box refinementFaceBox(refinementBox);
                          refinementFaceBox.setRange(idir, 0);
                          Real scaling =
                             1.0 / Real(refinementFaceBox.numPts());
                          FArrayBox otherFluxFab(fullFaceBox, ncomp);
                          otherFluxFab.setVal(0.);
                          for (BoxIterator bit(fullFaceBox); bit.ok(); ++bit)
                            {
                              IntVect iv = bit();
                              for (BoxIterator bitRef(refinementFaceBox); bitRef.ok(); ++bitRef)
                                {
                                  IntVect br = bitRef();
                                  for (int comp = 0; comp < ncomp; comp++)
                                    {
                                      otherFluxFab(iv, comp) +=
                                        scaling *
                                        fineFluxFab(refRatioVect*iv + br, comp);
                                    }
                                }
                            }
                          // otherFluxFab *= reorient;
                          defectFab.plus(otherFluxFab, fullFaceBox, fullFaceBox, 
                                         sign(side)*1., 0, 0, ncomp);
                          defectFab.shiftHalf(idir, -sign(side));
                        }
                    }
                  faceID++;
                }
            }
          const FArrayBox& UFab = Ucoarse[ditC];

          FArrayBox UshouldFab(baseBox, ncomp);
          UshouldFab.setVal(startVal);
          UshouldFab += defectFab;

          FArrayBox diffFab(baseBox, ncomp);
          diffFab.copy(UshouldFab);
          diffFab -= UFab;
          Real diffFabMaxNorm = diffFab.norm(baseBox, MAXNORM, 0, ncomp);
          if (verbosity >= 2)
            {
              pout() << " " << diffFabMaxNorm;
            }
          Real relativeMaxNorm = diffFabMaxNorm / fluxMaxNorm;
          Real thisRatio = relativeMaxNorm / precision;
          worstRatio = Max(worstRatio, thisRatio);
          if (thisRatio > 16)
            {
              status--;
            }
          CH_assert(0 == 0); // dummy for gdb
        }
      if (verbosity >= 2)
        {
          pout() << endl;
        }
      CH_assert(0 == 0); // dummy for gdb
    }
  if (verbosity >= 2)
    {
      pout() << "worst relative difference / precision = " << worstRatio << endl;
    }
  
  // Clean up.
  delete coordSysCoarse;
  delete coordSysFine;
  delete cssFactory;

  // Git.
  return status;
}

//-----------------------------------------------------------------------

int
cubedSphereVectorAvgTest()
{
  int interpDims = SpaceDim;
  if (SpaceDim == 3) interpDims = SpaceDim-1;

  int nRefRatio = 2;
  IntVect refRatioVect = IntVect::Unit;
  for (int idir = 0; idir < interpDims; idir++)
    {
      refRatioVect[idir] = nRefRatio;
    }
  Box refinementBox(IntVect::Zero, refRatioVect-IntVect::Unit);

  CubedSphereShellCSFactory* cssFactory = new CubedSphereShellCSFactory();
  cssFactory->setHeight(1.);
  cssFactory->setRadius(1.);
  RealVect probLen = RealVect(D_DECL(11*(M_PI/2.), M_PI/2., 1.));
  int nblocks = 6;

  Vector< List<int> > coarseBlockIndices = getCubedSphereCoarseBlockIndices();

  int ncomp = SpaceDim;
  Interval allIntv(0, ncomp-1);
  Interval vectorIntv = allIntv;

  int Nmin = 16;
  int Nmax = 256;
#if CH_SPACEDIM >= 3
  Nmax = 128;
#endif
  int status = 0;
  for (int icbi = 0; icbi < coarseBlockIndices.size(); icbi++)
    {
      const List<int>& coarseList = coarseBlockIndices[icbi];

      if (verbosity >= 1)
        {
          pout() << "========== Coarse blocks";
          for (ListIterator<int> li(coarseList); li.ok(); ++li)
            {
              int blockNum = coarseList[li];
              pout() << " " << blockNum;
            }
          pout() << endl;
        }

      Vector<Real> allErr;
      for (int N = Nmin; N < Nmax; N *= 2)
        {
          Real maxErr = 0.;

          // Create a cubed sphere coordinate system with six N-by-N blocks.
#if CH_SPACEDIM >= 3
          // Avoid irritating unused variable warning if dimension < 3.
          int nlayers = N/2;
#endif

          IntVect domainLen = IntVect(D_DECL(11*N, N, nlayers));
          ProblemDomain coarseDomain(IntVect::Zero, domainLen-IntVect::Unit);

          ProblemDomain fineDomain = refine(coarseDomain, refRatioVect);

          RealVect dxCoarse = probLen / RealVect(domainLen);
          RealVect dxFine = dxCoarse / RealVect(refRatioVect);
          MultiBlockCoordSys* coordSysCoarse =
            cssFactory->getCoordSys(coarseDomain, dxCoarse);
          MultiBlockCoordSys* coordSysFine =
            cssFactory->getCoordSys(fineDomain, dxFine);

          MultiBlockCoordSys* coordSysCoarseSub =
            cssFactory->getCoordSys(refine(coarseDomain, subRatioVect),
                                    dxCoarse / RealVect(subRatioVect));
          MultiBlockCoordSys* coordSysFineSub =
            cssFactory->getCoordSys(refine(fineDomain, subRatioVect),
                                    dxFine / RealVect(subRatioVect));

          const Vector<Box>& blockBoxesCoarse =
            coordSysCoarse->mappingBlocks();
          const Vector<Box>& blockBoxesFine =
            coordSysFine->mappingBlocks();

          const Vector<Box>& blockBoxesCoarseSub =
            coordSysCoarseSub->mappingBlocks();
          const Vector<Box>& blockBoxesFineSub =
            coordSysFineSub->mappingBlocks();

          const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
            coordSysCoarse->boundaries();

          IntVect baseLen = IntVect(D_DECL(N, N, nlayers));
          Box baseBox(IntVect::Zero, baseLen-IntVect::Unit);

          Vector<Box> coarseBoxes;
          Vector<Box> fineBoxes;

          for (int iblock = 0; iblock < nblocks; iblock++)
            {
              const Box& blockBox = blockBoxesCoarse[iblock];
              const IntVect& blockBoxLo = blockBox.smallEnd();
              if (coarseList.includes(iblock))
                {
                  coarseBoxes.push_back(blockBox);
                }
              else
                {
                  IntVect fineBlockBoxLo = refRatioVect * blockBoxLo;
                  for (BoxIterator bitRef(refinementBox); bitRef.ok(); ++bitRef)
                    {
                      IntVect br = bitRef();
                      IntVect thisFineShift = baseLen * br + fineBlockBoxLo;
                      Box thisFineBox(baseBox);
                      thisFineBox.shift(thisFineShift);
                      fineBoxes.push_back(thisFineBox);
                    }
                }
            }

          mortonOrdering(coarseBoxes);
          Vector<int> coarseProcs(coarseBoxes.size());
          LoadBalance(coarseProcs, coarseBoxes);
          DisjointBoxLayout coarseGrids(coarseBoxes, coarseProcs, coarseDomain);
          DataIterator ditC = coarseGrids.dataIterator();

          mortonOrdering(fineBoxes);
          Vector<int> fineProcs(fineBoxes.size());
          LoadBalance(fineProcs, fineBoxes);
          DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);
          DataIterator ditF = fineGrids.dataIterator();

          // Create a solution Ucoarse with startVal on
          // the coarse grids with a single ghost layer.
          Real startVal = 0.1;
          LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
          for (ditC.begin(); ditC.ok(); ++ditC)
            {
              Ucoarse[ditC].setVal(startVal);
            }

          // Create coarse and fine fluxes that straddle the
          // coarse-fine boundary.
          LevelData<FluxBox> coarseFlux(coarseGrids, ncomp);
          for (ditC.begin(); ditC.ok(); ++ditC)
            {
              const Box& baseBox = coarseGrids[ditC];
              int thisBlock = coordSysCoarse->whichBlock(baseBox);
              const Box& blockSubBox = blockBoxesCoarseSub[thisBlock];
              const NewCoordSys* thisCoordSysCoarseSubPtr =
                coordSysCoarseSub->getCoordSys(thisBlock);

              FluxBox& coarseFlub = coarseFlux[ditC];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  FArrayBox& coarseFluxFab = coarseFlub[idir];
                  const Box& faceBox = coarseFluxFab.box();
                  evalVectorAvgFunction(coarseFluxFab, faceBox, idir,
                                        thisCoordSysCoarseSubPtr, blockSubBox);
                }
            }

          LevelData<FluxBox> fineFlux(fineGrids, ncomp);
          LevelData<FluxBox> fineFlux2(fineGrids, ncomp, IntVect::Unit);
          // Don't just set to norm(fineFlux), because most of it is E+298.
          Real fluxMaxNorm = 0.; // norm(fineFlux, allIntv, MAXNORM);
          for (ditF.begin(); ditF.ok(); ++ditF)
            {
              const Box& baseBox = fineGrids[ditF];
              Box coarseBaseBox = coarsen(baseBox, refRatioVect);
              int thisBlock = coordSysCoarse->whichBlock(coarseBaseBox);
              const NewCoordSys* thisCoordSysFineSubPtr =
                coordSysFineSub->getCoordSys(thisBlock);
              const Box& fineBlockSubBox = blockBoxesFineSub[thisBlock];
              const Box& fineBlockBox = blockBoxesFine[thisBlock];
              FluxBox& fineFlub = fineFlux[ditF];
              FluxBox& fineFlub2 = fineFlux2[ditF];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  FArrayBox& fineFluxFab = fineFlub[idir];
                  FArrayBox& fineFlux2Fab = fineFlub2[idir];
                  // const Box& faceBox = fineFluxFab.box();
                  const Box& face2Box = fineFlux2Fab.box();
                  // Don't fill in fineFlux2Fab on normal-direction faces
                  // outside block.
                  Box fineBlockFaces = surroundingNodes(fineBlockBox, idir);
                  fineBlockFaces.grow(IntVect::Unit - BASISV(idir));
                  Box face2BoxIn = face2Box & fineBlockFaces;
                  // fineFluxFab.setVal(0.);
                  // evalVectorAvgFunction(fineFluxFab, faceBox, idir,
                  // thisCoordSysFineSubPtr, fineBlockSubBox);
                  evalVectorAvgFunction(fineFlux2Fab, face2BoxIn, idir,
                                        thisCoordSysFineSubPtr, fineBlockSubBox);
                  fineFluxFab.copy(fineFlux2Fab);
                  fluxMaxNorm =
                    Max(fluxMaxNorm, maxnormFaceFab(fineFluxFab, fineBlockBox));
                }
            }
          reduceMax(fluxMaxNorm); // Get global max across processors.

          // Create the flux register.
          MultiBlockFluxRegister fluxReg(coordSysCoarse, coordSysFine);
          // If SpaceDim = 3, don't refine in the r (last) direction (anisotropic refinement)
          int unrefDir = (SpaceDim == 3) ? SpaceDim-1 : -1;
          bool scaleFlag = true;
          fluxReg.define(fineGrids, coarseGrids, coarseDomain,
                         nRefRatio, ncomp, scaleFlag, unrefDir);
          fluxReg.setToZero();
                  
          // Stick the coarse and fine fluxes into the register
          // to compute the defect.
          Real scale1 = 1.0;
          for (ditC.begin(); ditC.ok(); ++ditC)
            {
              const FluxBox& coarseFlub = coarseFlux[ditC];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  fluxReg.incrementCoarse(coarseFlub[idir], scale1, ditC(),
                                          allIntv, allIntv, idir);
                }
            }

          for (ditF.begin(); ditF.ok(); ++ditF)
            {
              const FluxBox& fineFlub = fineFlux[ditF];
              const FluxBox& fineFlub2 = fineFlux2[ditF];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  fluxReg.incrementFineAveraged(fineFlub[idir], fineFlub2[idir],
                                                scale1, ditF(),
                                                vectorIntv, allIntv, allIntv, idir);
                }
            }
      
          // Apply the defect to the solution Ucoarse.
          fluxReg.reflux(Ucoarse, vectorIntv, scale1);

          // Find what the defect and Ucoarse should be.
          for (ditC.begin(); ditC.ok(); ++ditC)
            {
              const FluxBox& coarseFlub = coarseFlux[ditC];
              
              const Box& baseBox = coarseGrids[ditC];
              int thisBlock = coordSysCoarse->whichBlock(baseBox);
              const Box& blockBox = blockBoxesCoarse[thisBlock];
              const Box& fineBlockSubBox = blockBoxesFineSub[thisBlock];
              const NewCoordSys* thisCoordSysFineSubPtr =
                coordSysFineSub->getCoordSys(thisBlock);
              const Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries =
                boundaries[thisBlock];

              FArrayBox defectFab(baseBox, ncomp);
              defectFab.setVal(0.);

              int faceID = 0;
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Side::LoHiSide side = sit();
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      if (idir < interpDims)
                        {
                          // Check that Ucoarse is correct on this face.
                          const BlockBoundary& bb = blockBoundaries[faceID];
                          int otherBlock = bb.neighbor();
                          if ( !coarseList.includes(otherBlock) )
                            {
                              // There is fine patch on other side of this face.
                              
                              IndicesTransformation it = bb.getTransformation();
                              
                              int idirOther = bb.dirOther(idir);
                              int reorient = bb.reorientFace(idir);
                              Side::LoHiSide sideOther =
                                (reorient == -1) ? side : flip(side);
                              
                              Box adjBox = adjCellBox(blockBox, idir, side, 1);
                              Box adjOtherBox = it.transform(adjBox);
                              
                              // coarse flux on fullFaceBox
                              Box fullFaceBox(adjBox);
                              fullFaceBox.shiftHalf(idir, -sign(side));
                              
                              // fine flux on fullFaceOtherBox
                              Box fullFaceOtherBox(adjOtherBox);
                              fullFaceOtherBox.shiftHalf(idirOther, sideOther);
                              
                              // Ucoarse adjusted on fullCellBox
                              Box fullCellBox(adjBox);
                              fullCellBox.shift(idir, -sign(side));

                              // FIXME: This is average flux on coarse face.
                              // How accurate is it?
                              const FArrayBox& coarseFluxFab = coarseFlub[idir];
                              defectFab.shiftHalf(idir, sign(side));
                              defectFab.plus(coarseFluxFab, fullFaceBox, fullFaceBox, 
                                             -sign(side)*1., 0, 0, ncomp);
                              
                              // Add in sign(side) * reorient * fine flux.
                              
                              // face-centered box
                              Box fullFaceFineBox =
                                refine(fullFaceBox, refRatioVect);
                              FArrayBox fineFluxFab(fullFaceFineBox, ncomp);
                              evalVectorAvgFunction(fineFluxFab,
                                                    fullFaceFineBox,
                                                    idir,
                                                    thisCoordSysFineSubPtr,
                                                    fineBlockSubBox);
                              Box refinementFaceBox(refinementBox);
                              refinementFaceBox.setRange(idir, 0);
                              Real scaling =
                                1.0 / Real(refinementFaceBox.numPts());
                              FArrayBox otherFluxFab(fullFaceBox, ncomp);
                              otherFluxFab.setVal(0.);
                              for (BoxIterator bit(fullFaceBox); bit.ok(); ++bit)
                                {
                                  IntVect iv = bit();
                                  for (BoxIterator bitRef(refinementFaceBox); bitRef.ok(); ++bitRef)
                                    {
                                      IntVect br = bitRef();
                                      for (int comp = 0; comp < ncomp; comp++)
                                        {
                                          otherFluxFab(iv, comp) +=
                                            scaling *
                                            fineFluxFab(refRatioVect*iv + br, comp);
                                        }
                                    }
                                }
                              // otherFluxFab *= reorient;
                              defectFab.plus(otherFluxFab, fullFaceBox, fullFaceBox, 
                                             sign(side)*1., 0, 0, ncomp);
                              defectFab.shiftHalf(idir, -sign(side));
                            }
                        }
                      faceID++;
                    }
                }
              const FArrayBox& UFab = Ucoarse[ditC];
              
              FArrayBox UshouldFab(baseBox, ncomp);
              UshouldFab.setVal(startVal);
              UshouldFab += defectFab;
              
              FArrayBox diffFab(baseBox, ncomp);
              diffFab.copy(UshouldFab);
              diffFab -= UFab;
              RealVect diffFabMaxNorm;
              Real diffFabMaxNormAll = 0.;
              for (int comp = 0; comp < SpaceDim; comp++)
                {
                  diffFabMaxNorm[comp] =
                    diffFab.norm(baseBox, MAXNORM, comp, 1);
                  diffFabMaxNormAll =
                    Max(diffFabMaxNormAll, diffFabMaxNorm[comp]);
                }
              if (verbosity >= 1)
                {
                  if (verbosity >= 3) pout() << "block=" << thisBlock
                                             << " maxnorm";
                  pout() << " " << diffFabMaxNorm;
                }
              if (verbosity >= 3)
                {
                  // Look at difference on each face.
                  int faceIDdebug = 0;
                  for (SideIterator sitdebug; sitdebug.ok(); ++sitdebug)
                    {
                      Side::LoHiSide sidedebug = sitdebug();
                      for (int idirdebug = 0; idirdebug < SpaceDim; idirdebug++)
                        {
                          Box adjWithinBox =
                            adjCellBox(baseBox, idirdebug, sidedebug, -1);
                          // Omit corners.
                          adjWithinBox.grow(BASISV(idirdebug) - IntVect::Unit);
                          RealVect diffFabFaceMaxNorm;
                          for (int comp = 0; comp < SpaceDim; comp++)
                            {
                              diffFabFaceMaxNorm[comp] =
                                diffFab.norm(adjWithinBox, 0, comp, 1);
                            }
                          char letts[] = {'x', 'y', 'z'};
                          pout() << " " << letts[idirdebug]
                                 << ((sidedebug == Side::Lo) ? "Lo" : "Hi")
                                 << "=" << diffFabFaceMaxNorm;
                          faceIDdebug++;
                        }
                    }
                }
              if (verbosity >= 3) pout() << endl;
              
              maxErr = Max(maxErr, diffFabMaxNormAll/fluxMaxNorm);
            }
          reduceMax(maxErr); // Get global max across processors.

          if (verbosity == 1)
            {
              pout() << endl;
            }

          if (verbose) pout() << " " << N << " error = " << maxErr;
          if (verbosity >= 3)
            {
              pout() << " (relative; abs = " << (maxErr * fluxMaxNorm) << ")";
            }
          int allErrSize = allErr.size();
          if (allErrSize > 0)
            {
              Real prevErr = allErr[allErrSize-1];
              Real rate = log2(prevErr/maxErr);
              if (verbose) pout() << " rate " << rate;
              // Check for 4th-order convergence.
              if (rate < 0.9 * 4) status--;
            }
          if (verbose) pout() << endl;
          allErr.push_back(maxErr);

          // Clean up.
          delete coordSysCoarse;
          delete coordSysFine;
          delete coordSysCoarseSub;
          delete coordSysFineSub;
        }
    }
  
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

  int statusAll = 0;
  int status = 0;

  status = singleBlockTest();
  if ( status == 0 )
    cout << indent << pgmname << ": singleBlockTest passed." << endl ;
  else
    cout << indent << pgmname << ": singleBlockTest failed with return code " << status << endl ;
  statusAll += status;

  status = tripleCartesianTest();
  if ( status == 0 )
    cout << indent << pgmname << ": tripleCartesianTest passed." << endl ;
  else
    cout << indent << pgmname << ": tripleCartesianTest failed with return code " << status << endl ;
  statusAll += status;

  status = multiCartesianTest();
  if ( status == 0 )
    cout << indent << pgmname << ": multiCartesianTest passed." << endl ;
  else
    cout << indent << pgmname << ": multiCartesianTest failed with return code " << status << endl ;
  statusAll += status;

#if CH_SPACEDIM >= 2
  status = cubedSphereCellTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereCellTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereCellTest failed with return code " 
      << status << endl ;
  statusAll += status;

  // Test MultiBlockFluxRegister on CubedSphereShellCS
  // calling incrementFine() (NOT incrementFine{Centered|Averaged})
  // with face-CENTERED fluxes, scalar fluxes only.
  status = cubedSphereWholeTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereWholeTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereWholeTest failed with return code " 
      << status << endl ;
  statusAll += status;

  // Test MultiBlockFluxRegister on CubedSphereShellCS
  // calling incrementFineCentered(),
  // with face-CENTERED fluxes, scalar fluxes only.
  status = cubedSphereScalarTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereScalarTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereScalarTest failed with return code " 
      << status << endl ;
  statusAll += status;

  // Test MultiBlockFluxRegister on CubedSphereShellCS
  // calling incrementFineCentered(),
  // with face-CENTERED fluxes, vector fluxes only.
  status = cubedSphereVectorTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereVectorTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereVectorTest failed with return code " 
      << status << endl ;
  statusAll += status;

  // Test MultiBlockFluxRegister on CubedSphereShellCS
  // calling incrementFineCentered(),
  // with face-CENTERED fluxes, both vector and scalar fluxes.
  status = cubedSphereVectorScalarTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereVectorScalarTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereVectorScalarTest failed with return code " 
      << status << endl ;
  statusAll += status;

  // Test MultiBlockFluxRegister on CubedSphereShellCS
  // calling incrementFineAveraged(),
  // with face-AVERAGED fluxes, scalar fluxes only.
  status = cubedSphereScalarAvgTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereScalarAvgTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereScalarAvgTest failed with return code " 
      << status << endl ;
  statusAll += status;

#ifdef CH_USE_DOUBLE
  // This is a convergence test that is useless in single precision
  // because error on coarsest example is already close to machine precision,
  // so it won't get much better by refining.
  // Test MultiBlockFluxRegister on CubedSphereShellCS
  // calling incrementFineAveraged(),
  // with face-AVERAGED fluxes, vector fluxes only.
  // This is by necessity a convergence test, because it won't be exact.
  status = cubedSphereVectorAvgTest();
  if ( status == 0 )
    cout << indent << pgmname << ": cubedSphereVectorAvgTest passed." << endl ;
  else
    cout << indent << pgmname << ": cubedSphereVectorAvgTest failed with return code " 
      << status << endl ;
  statusAll += status;
#endif

#endif

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize ();
#endif
  return statusAll ;
}
//-----------------------------------------------------------------------
