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

#include "parstream.H"
#include "MultiBlockFluxRegister.H"
#include "SingleBlockCSAdaptor.H"
#include "CartesianBlockCS.H"
#include "DoubleCartesianCS.H"
#include "TripleCartesianCS.H"
#include "CubedSphere2DCS.H"
#include "CubedSphereShellCS.H"
#include "RefCountedPtr.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "DebugOut.H"
#include "DebugDump.H"
#include "FABView.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testMultiBlockFluxRegister";
static const char* indent = "   ";
static const char* indent2 = "      ";
//static bool verbose = false ;
static bool verbose = true ;

#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif

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
int
singleBlockTest(bool anisotropic, const IntVect& a_refRat, const int& a_unrefDir)
{
  // NOTE: This test demonstrates that the multi-block flux register can
  // NOTE: perform as a single-block flux register. Note that
  // NOTE: the coarse-fine grids here are properly nested in this test.
  // NOTE: This is not so with the multi-block tests, since the boxes
  // NOTE: in multi-block calculations are necessarily disjoint. So the
  // NOTE: logic is different for single- and multi-block cases.

  // Create a single-block Cartesian coordinate system.
  int ncomp = 1;
  IntVect zero;
  RealVect dx = 1.0*RealVect::Unit;
  Box domainBox(IntVect::Zero, 16*IntVect::Unit);
  ProblemDomain coarseDomain(domainBox);
  NewCoordSys* cartesianCS = new CartesianBlockCS(0, zero, dx, domainBox);
  // Note that coordSys takes control of cartesianCS!
  RefCountedPtr<MultiBlockCoordSys > coordSys(new SingleBlockCSAdaptor(cartesianCS, coarseDomain));

  // Create the coarse-level box (two grid cells at (7,7), (8,7)).
  Vector<Box> coarseBoxes;
  Box coarseBox = domainBox;
  coarseBoxes.push_back(coarseBox);
  Vector<int> procs(1,0);
  DisjointBoxLayout coarseGrids(coarseBoxes, procs, coarseDomain);

  // Create the fine-level box: a 2 x 2 square spanning (14, 14) x (15, 15).
  Vector<Box> fineBoxes;
  IntVect If = 14*IntVect::Unit;
  Box fineBox(If, If + IntVect::Unit);
  fineBoxes.push_back(fineBox);
  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(a_refRat);
  DisjointBoxLayout fineGrids(fineBoxes, procs, fineDomain);

  RealVect dxFine = dx / RealVect(a_refRat);
  Box domainBoxFine = refine(domainBox, a_refRat);
  NewCoordSys* cartesianCSFine = new CartesianBlockCS(0, zero, dxFine, domainBoxFine);
  RefCountedPtr<MultiBlockCoordSys > coordSysFine(new SingleBlockCSAdaptor(cartesianCSFine, fineDomain));

  // Create the flux register.
  MultiBlockFluxRegister fluxReg(&(*coordSys), &(*coordSysFine));
  fluxReg.define(fineGrids, coarseGrids, fineDomain, 2, ncomp, true, a_unrefDir);
  fluxReg.setToZero();

  // Create a zero solution U the coarse grids with a single ghost layer.
  LevelData<FArrayBox> Ucoarse(coarseGrids, ncomp, IntVect::Unit);
  Real uval = 4;
  Real fval = 7;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    Ucoarse[dit()].setVal(uval);

  // Create coarse and fine fluxes that straddle the coarse-fine boundary.
  LevelData<FluxBox> coarseFlux(coarseGrids, ncomp),
    fineFlux(fineGrids, ncomp);
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    coarseFlux[dit()].setVal(fval);
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    fineFlux[dit()].setVal(fval);

  fluxReg.setToZero();
  // Stick the coarse and fine fluxes into the register to compute the defect.
  Interval i(0, ncomp-1);
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          fluxReg.incrementCoarse(coarseFlux[dit()][0], 1.0, dit(), i, i, idir,Side::Lo);
          fluxReg.incrementCoarse(coarseFlux[dit()][0], 1.0, dit(), i, i, idir,Side::Hi);
        }
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          fluxReg.incrementFine(fineFlux[dit()][idir], 1.0, dit(), i, i, idir, Side::Lo);
          fluxReg.incrementFine(fineFlux[dit()][idir], 1.0, dit(), i, i, idir, Side::Hi);
        }
    }

//  // Apply the defect to the solution.
//  fluxReg.poutCoarseRegisters();
//  fluxReg.poutFineRegisters();

  fluxReg.reflux(Ucoarse, Interval(), 1.0);

  // Check the answer.
  int status = 0;
  for (DataIterator dit = coarseGrids.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox& U = Ucoarse[dit()];
    for (BoxIterator bit(coarseGrids[dit()]); bit.ok(); ++bit)
      {
        IntVect Ig =  bit();
        Ig[0] += 1;
        for (int icomp = 0; icomp < ncomp; icomp++)
          {
            if (Abs(U(Ig, icomp) - uval) > precision)
              {
                status = -1;
                break;
              }
          }
      }

  }

  // Git.
  return status;
}
//-----------------------------------------------------------------------
int getFineBoxes(Vector<Box>& a_fineBoxes,
                 const Vector<Box>& a_coarBoxes,
                 const IntVect& a_refRat,
                 const int & a_boxsize,
                 const int & a_unrefDir,
                 const int & a_case)
{
  if (a_case < a_coarBoxes.size())
    {
      pout() << "refining and shrinking the box low end of box " << a_case << endl;
      Box oneBox  = a_coarBoxes[a_case];
      oneBox.refine(a_refRat);
      for (int idir = 0; idir < SpaceDim; idir++) 
        {
          if (idir != a_unrefDir)
            {
              oneBox.growHi(idir, -a_boxsize);
            }
        }
      a_fineBoxes.resize(1, oneBox);
      return a_case+1;
    }
  else if (a_case <= (2*a_coarBoxes.size()-1))
    {
      int ibox = a_case - a_coarBoxes.size();
      pout() << "refining and shrinking  the box high end of box " << ibox << endl;
      Box oneBox  = a_coarBoxes[ibox];
      oneBox.refine(a_refRat);
      for (int idir = 0; idir < SpaceDim; idir++) 
        {
          if (idir != a_unrefDir)
            {
              oneBox.growLo(idir, -a_boxsize);
            }
        }
      a_fineBoxes.resize(1, oneBox);
      return a_case+1;
    }
  else if (a_case == (2*a_coarBoxes.size()))
    {
      pout() << "refining and shrinking  the box low end all boxes" << endl;
      a_fineBoxes.resize(a_coarBoxes.size());
      
      for (int ibox = 0; ibox < a_coarBoxes.size(); ibox++)
        {
          Box oneBox  = a_coarBoxes[ibox];
          oneBox.refine(a_refRat);
          for (int idir = 0; idir < SpaceDim; idir++) 
            {
              if (idir != a_unrefDir)
                {
                  oneBox.growLo(idir, -a_boxsize);
                }
            }
          a_fineBoxes[ibox]= oneBox;
        }
      return a_case+1;
    }
  else if (a_case == (2*a_coarBoxes.size()+1))
    {
      pout() << "refining and shrinking  the box high end all boxes" << endl;
      a_fineBoxes.resize(a_coarBoxes.size());
      
      for (int ibox = 0; ibox < a_coarBoxes.size(); ibox++)
        {
          Box oneBox  = a_coarBoxes[ibox];
          oneBox.refine(a_refRat);
          for (int idir = 0; idir < SpaceDim; idir++) 
            {
              if (idir != a_unrefDir)
                {
                  oneBox.growHi(idir, -a_boxsize);
                }
            }
          a_fineBoxes[ibox]= oneBox;
        }
      return -1; //no more cases
    }
  else
    {
      MayDay::Error("bad case number for cubed sphere test");
    }
  return -1;
}
//-----------------------------------------------------------------------
int
doubleCartesianTest(bool anisotropic, const IntVect& a_refRat, const int& a_unrefDir)
{
  // Create a double-Cartesian coordinate system.
  // This creates a 2x2 set of N/3 x N/3 patches in indexspace,
  // separated by buffers of width N/3.

  int boxsize = 16;
  int ncomp =4;
  RealVect dx = 1.0*RealVect::Unit;


  Box domainBox(IntVect::Zero, (3*boxsize-1)*IntVect::Unit);
  ProblemDomain coarDomain(domainBox);
  ProblemDomain fineDomain = coarDomain;
  fineDomain.refine(a_refRat);
  DoubleCartesianCS* coordSysDC = new DoubleCartesianCS();
  coordSysDC->define(coarDomain, dx);
  RefCountedPtr<MultiBlockCoordSys > coordSys(dynamic_cast<MultiBlockCoordSys* > (coordSysDC));

  // Hm, in anisotropic case, it's complicated, because you can't
  // refine the coordinate system by a_refRat:
  // DoubleCartesianCS requires a square domain.
  IntVect refRatFull = a_refRat;
  ProblemDomain fineDomainFull = fineDomain;
  if (anisotropic)
    {
      int refDir = (a_unrefDir == 0) ? 1 : 0;
      refRatFull[a_unrefDir] = refRatFull[refDir];
      fineDomainFull = refine(coarDomain, refRatFull);
    }
  RealVect dxFine = dx / RealVect(refRatFull);
  DoubleCartesianCS* coordSysDCFine = new DoubleCartesianCS();
  coordSysDCFine->define(fineDomainFull, dxFine);
  RefCountedPtr<MultiBlockCoordSys > coordSysFine(dynamic_cast<MultiBlockCoordSys* > (coordSysDCFine));

  //create a coarse dbl over the entire grid
  Vector<Box> coarBoxes = coordSys->mappingBlocks();
  Vector<int> coarProcs;
  LoadBalance(coarProcs, coarBoxes);
  DisjointBoxLayout coarGrids(coarBoxes, coarProcs, coarDomain);

  //create the fine level boxes.
  int icase = 0;
  while (icase >= 0)
    {
      Vector<Box> fineBoxes;
      Vector<int> fineProcs;
      //increments icase until no more cases, then sets it 
      int inextcase = getFineBoxes(fineBoxes, coarBoxes, a_refRat, boxsize, a_unrefDir, icase);
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);

      // Create the flux register.
      MultiBlockFluxRegister fluxReg(&(*coordSys), &(*coordSysFine));
      fluxReg.define(fineGrids, coarGrids, fineDomain, 2, ncomp, true);
      fluxReg.setToZero();
      Vector<Real> uval(ncomp), fval(ncomp);
      for (int icomp = 0; icomp < ncomp; icomp++)
        {
          fval[icomp] = 7*(1 + icomp);
          uval[icomp] = 4*(1 + icomp);
        }
      // Create a zero solution U the coarse grids with a single ghost layer.
      LevelData<FArrayBox> Ucoar(coarGrids, ncomp, IntVect::Unit);
      LevelData<FluxBox> coarFlux(coarGrids, ncomp), fineFlux(fineGrids, ncomp);
      for (int icomp = 0; icomp < ncomp; icomp++)
        {
          for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
            {
              Ucoar[dit()].setVal(uval[icomp],icomp);
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  coarFlux[dit()][idir].setVal(fval[icomp],icomp);
                }
            }
          for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
            {
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  fineFlux[dit()][idir].setVal(fval[icomp], icomp);
                }
            }
        }

      // Stick the coarse and fine fluxes into the register to compute the defect.
      Interval interv(0, ncomp-1);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
                {
                  fluxReg.incrementCoarse(coarFlux[dit()][idir], 1.0, dit(), interv, interv, idir, sit());
                }
              for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
                {
                  fluxReg.incrementFine(fineFlux[dit()][idir], 1.0, dit(), interv, interv, idir, sit());
                }
            }
        }

      // Apply the defect to the solution.
//      fluxReg.poutCoarseRegisters();
//      fluxReg.poutFineRegisters();
//      //interval means no vectors
      fluxReg.reflux(Ucoar, Interval(), 1.0);

      // Check the answer.
      int status = 0;
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          const FArrayBox& U = Ucoar[dit()];
          for (BoxIterator bit(coarGrids[dit()]); bit.ok(); ++bit)
            {
              IntVect Ig =  bit();
              for (int ivar = 0; ivar < ncomp; ivar++)
                {
                  Real calc =  U(Ig, ivar);
                  if (Abs(calc - uval[ivar]) > precision)
                    status = -1;
                  if (status != 0)
                    {
                      pout() << "double cartesian case failed for case = " << icase << ", cell = " <<  Ig  << endl;
                      return -1;
                    }
                }
            }
        }
      icase = inextcase;
    }

  // Git.
  return 0;
}
//-----------------------------------------------------------------------
int
tripleCartesianTest(bool anisotropic, const IntVect& a_refRat, const int& a_unrefDir)
{
  // Create a triple-Cartesian coordinate system.
  // This creates a 3 x 3 set of N/5 x N/5 patches in indexspace,
  // separated by buffers of width N/5.
  int boxsize = 16;
  RealVect dx = 1.0*RealVect::Unit;
  Box domainBox(IntVect::Zero, (5*boxsize-1)*IntVect::Unit);
  ProblemDomain coarDomain(domainBox);
  ProblemDomain fineDomain = coarDomain;
  fineDomain.refine(2);
  TripleCartesianCS* coordSysDC = new TripleCartesianCS();
  coordSysDC->define(coarDomain, dx);
  RefCountedPtr<MultiBlockCoordSys > coordSys(dynamic_cast<MultiBlockCoordSys* > (coordSysDC));

  RealVect dxFine = dx / RealVect(a_refRat);
  TripleCartesianCS* coordSysDCFine = new TripleCartesianCS();
  coordSysDCFine->define(fineDomain, dxFine);
  RefCountedPtr<MultiBlockCoordSys > coordSysFine(dynamic_cast<MultiBlockCoordSys* > (coordSysDCFine));

  //create a coarse dbl over the entire grid
  Vector<int> coarProcs;
  Vector<Box> coarBoxes = coordSys->mappingBlocks();
  LoadBalance(coarProcs, coarBoxes);
  DisjointBoxLayout coarGrids(coarBoxes, coarProcs, coarDomain);

  //create the fine level boxes.
  int icase = 0;
  while (icase >= 0)
    {
      Vector<Box> fineBoxes;
      Vector<int> fineProcs;
      //increments icase until no more cases, then sets it 
      int inextcase = getFineBoxes(fineBoxes, coarBoxes, a_refRat, boxsize, a_unrefDir, icase);
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);

      // Create the flux register.
      MultiBlockFluxRegister fluxReg(&(*coordSys), &(*coordSysFine));
      fluxReg.define(fineGrids, coarGrids, fineDomain, 2, 1, true);
      fluxReg.setToZero();
      Real uval = 4;
      Real fval = 7;
      // Create a zero solution U the coarse grids with a single ghost layer.
      LevelData<FArrayBox> Ucoar(coarGrids, 1, IntVect::Unit);
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        Ucoar[dit()].setVal(uval);

      // Create coarse and fine fluxes that straddle the coarse-fine boundary.
      LevelData<FluxBox> coarFlux(coarGrids, 1), fineFlux(fineGrids, 1);

      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        coarFlux[dit()].setVal(fval);
      for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
        fineFlux[dit()].setVal(fval);

      // Stick the coarse and fine fluxes into the register to compute the defect.
      Interval interv(0, 0);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
                {
                  fluxReg.incrementCoarse(coarFlux[dit()][idir], 1.0, dit(), interv, interv, idir, sit());
                }
              for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
                {
                  fluxReg.incrementFine(fineFlux[dit()][idir], 1.0, dit(), interv, interv, idir, sit());
                }
            }
        }

      // Apply the defect to the solution.
//      fluxReg.poutCoarseRegisters();
//      fluxReg.poutFineRegisters();
      //interval means no vector stuff
      fluxReg.reflux(Ucoar, Interval(), 1.0);

      // Check the answer.
      int status = 0;
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          const FArrayBox& U = Ucoar[dit()];
          for (BoxIterator bit(coarGrids[dit()]); bit.ok(); ++bit)
            {
              IntVect Ig =  bit();
              Real calc =  U(Ig, 0);
              if (Abs(calc - uval) > precision)
                status = -1;
              if (status != 0)
                {
                  pout() << "triple cartesian case failed for case = " << icase << ", cell = " <<  Ig  << endl;
                  return -1;
                }
            }
        }
      icase = inextcase;
    }

  // Git.
  return 0;
}


//-----------------------------------------------------------------------
int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  bool justGetOut = false;
#ifdef CH_USE_FLOAT
  justGetOut = true;
#endif
  if (SpaceDim==1) 
  {
    justGetOut = true;
  }

  int overallStatus = 0;

  if (!justGetOut)
  {
    parseTestOptions( argc ,argv ) ;
    if ( verbose )
      pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

    for (int isotropic = 0; isotropic < 2; isotropic++)
      {
        bool anisotropic = (isotropic==0);
        int unrefDir= -1;
        IntVect refRat = 2*IntVect::Unit;
        if (anisotropic)
          {
            refRat[SpaceDim-1] = 1;
            unrefDir =  SpaceDim-1;
          }

        int status = 0;
        status = singleBlockTest(anisotropic, refRat, unrefDir);
        if ( status == 0 )
        {
          cout << indent << pgmname << ": singleBlockTest passed." << endl ;
        }
        else
        {
          cout << indent << pgmname << ": singleBlockTest failed with return code " << status << " for isotropic =  " << isotropic << endl ;
          overallStatus += 1;
        }

        status = doubleCartesianTest(anisotropic, refRat, unrefDir);
        if ( status == 0 )
        {
          cout << indent << pgmname << ": doubleCartesianTest passed." << endl ;
        }
        else
        {
          cout << indent << pgmname << ": doubleCartesianTest failed with return code " << status << " for isotropic =  " << isotropic << endl ;
          overallStatus += 1;
        }


        status = tripleCartesianTest(anisotropic, refRat, unrefDir);
        if ( status == 0 )
        {
          cout << indent << pgmname << ": tripleCartesianTest passed." << endl ;
        }
        else
        {
          cout << indent << pgmname << ": tripleCartesianTest failed with return code " << status << " for isotropic =  " << isotropic << endl ;
          overallStatus += 1;
        }
      }
    if (overallStatus == 0)
    {
      pout() << "all testOldMBFR tests pass " << endl;
    }
  }
#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return overallStatus;
}
//-----------------------------------------------------------------------
