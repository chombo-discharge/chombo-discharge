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
#include "NewMultiBlockFluxRegister.H"
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
static const char* pgmname = "testNewMultiBlockFluxRegister";
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
singleBlockTest()
{
  // NOTE: This test demonstrates that the multi-block flux register can
  // NOTE: perform as a single-block flux register. Note that
  // NOTE: the coarse-fine grids here are properly nested in this test.
  // NOTE: This is not so with the multi-block tests, since the boxes
  // NOTE: in multi-block calculations are necessarily disjoint. So the
  // NOTE: logic is different for single- and multi-block cases.

  // Create a single-block Cartesian coordinate system.
  int ncomp = 4;
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
  fineDomain.refine(2);
  DisjointBoxLayout fineGrids(fineBoxes, procs, fineDomain);

  // Create the flux register.
  NewMultiBlockFluxRegister fluxReg(coordSys, fineGrids, coarseGrids, 
                                 2, ncomp, true);
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

  // Apply the defect to the solution.
  fluxReg.poutCoarseRegisters();
  fluxReg.poutFineRegisters();

  fluxReg.reflux(Ucoarse, i, 1.0);

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
int getDCSFineBoxes(Vector<Box>& a_fineBoxes,
                    const int  & a_boxsize,
                    const int  & a_case)
{
  a_fineBoxes.clear();
  int retval;
 
  if (a_case == 0)
    {
      pout() << "case 0 box in lower left" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 1;
    }
  else if (a_case == 1)
    {
      pout() << "case 1 box in lower right" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      loiv += (3*a_boxsize*BASISV(0));
      hiiv += (3*a_boxsize*BASISV(0));

      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 2;
    }
  else if (a_case == 2)
    {
      pout() << "case 2 box in upper left" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      loiv += (3*a_boxsize*BASISV(1));
      hiiv += (3*a_boxsize*BASISV(1));

      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 3;
    }
  else if (a_case == 3)
    {
      pout() << "case 3 box in upper right" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      loiv += (3*a_boxsize*BASISV(0));
      hiiv += (3*a_boxsize*BASISV(0));
      loiv += (3*a_boxsize*BASISV(1));
      hiiv += (3*a_boxsize*BASISV(1));
      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 4; //no more cases
    }
  else if (a_case == 4)
    {
      pout() << "case 3-- boxes in three corners" << endl;
      a_fineBoxes.resize(3);
      
      IntVect baseloiv = a_boxsize*IntVect::Unit;
      IntVect basehiiv = (2*a_boxsize-1)*IntVect::Unit;
      {
        IntVect loiv = baseloiv;
        IntVect hiiv = basehiiv;
        a_fineBoxes[0] = Box(loiv, hiiv);
      }

      {
        IntVect loiv = baseloiv;
        IntVect hiiv = basehiiv;
        loiv += (3*a_boxsize*BASISV(0));
        hiiv += (3*a_boxsize*BASISV(0));
        a_fineBoxes[1] = Box(loiv, hiiv);
      }

      {
        IntVect loiv = baseloiv;
        IntVect hiiv = basehiiv;
        loiv += (3*a_boxsize*BASISV(1));
        hiiv += (3*a_boxsize*BASISV(1));
        a_fineBoxes[2] = Box(loiv, hiiv);
      }

      retval = -1; //no more cases
    }
  else
    {
      MayDay::Error("bogus case in getDCSFineBoxes");
    }
  return retval;
}
//-----------------------------------------------------------------------
int
doubleCartesianTest()
{
  // Create a double-Cartesian coordinate system.
  // This creates a 2x2 set of N/3 x N/3 patches in indexspace,
  // separated by buffers of width N/3.

  int boxsize = 16;
  int ncomp =4;
  RealVect dx = 1.0*RealVect::Unit;
  Box domainBox(IntVect::Zero, 3*boxsize*IntVect::Unit-1);
  ProblemDomain coarDomain(domainBox);
  ProblemDomain fineDomain = coarDomain;
  fineDomain.refine(2);
  DoubleCartesianCS* coordSysDC = new DoubleCartesianCS();
  coordSysDC->define(coarDomain, dx);
  RefCountedPtr<MultiBlockCoordSys > coordSys(dynamic_cast<MultiBlockCoordSys* > (coordSysDC));

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
      int inextcase = getDCSFineBoxes(fineBoxes, boxsize, icase);
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);

      // Create the flux register.
      NewMultiBlockFluxRegister fluxReg(coordSys, fineGrids, coarGrids, 2, ncomp, true);
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
      fluxReg.poutCoarseRegisters();
      fluxReg.poutFineRegisters();
      fluxReg.reflux(Ucoar, interv, 1.0);

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
int getTCSFineBoxes(Vector<Box>& a_fineBoxes,
                    const int  & a_boxsize,
                    const int  & a_case)
{
  a_fineBoxes.clear();
  int retval;
 
  if (a_case == 0)
    {
      pout() << "case 0 box in lower left" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 1;
    }
  else if (a_case == 1)
    {
      pout() << "case 1 box in lower right" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      loiv += (3*a_boxsize*BASISV(0));
      hiiv += (3*a_boxsize*BASISV(0));

      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 2;
    }
  else if (a_case == 2)
    {
      pout() << "case 2 box in upper left" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      loiv += (3*a_boxsize*BASISV(1));
      hiiv += (3*a_boxsize*BASISV(1));

      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 3;
    }
  else if (a_case == 3)
    {
      pout() << "case 3 box in upper right" << endl;
      a_fineBoxes.resize(1);
      IntVect loiv = a_boxsize*IntVect::Unit;
      IntVect hiiv = (2*a_boxsize-1)*IntVect::Unit;
      loiv += (3*a_boxsize*BASISV(0));
      hiiv += (3*a_boxsize*BASISV(0));
      loiv += (3*a_boxsize*BASISV(1));
      hiiv += (3*a_boxsize*BASISV(1));
      a_fineBoxes[0] = Box(loiv, hiiv);
      retval = 4; //no more cases
    }
  else if (a_case == 4)
    {
      pout() << "case 3-- boxes in three corners" << endl;
      a_fineBoxes.resize(3);
      
      IntVect baseloiv = a_boxsize*IntVect::Unit;
      IntVect basehiiv = (2*a_boxsize-1)*IntVect::Unit;
      {
        IntVect loiv = baseloiv;
        IntVect hiiv = basehiiv;
        a_fineBoxes[0] = Box(loiv, hiiv);
      }

      {
        IntVect loiv = baseloiv;
        IntVect hiiv = basehiiv;
        loiv += (3*a_boxsize*BASISV(0));
        hiiv += (3*a_boxsize*BASISV(0));
        a_fineBoxes[1] = Box(loiv, hiiv);
      }

      {
        IntVect loiv = baseloiv;
        IntVect hiiv = basehiiv;
        loiv += (3*a_boxsize*BASISV(1));
        hiiv += (3*a_boxsize*BASISV(1));
        a_fineBoxes[2] = Box(loiv, hiiv);
      }

      retval = -1; //no more cases
    }
  else
    {
      MayDay::Error("bogus case in getDCSFineBoxes");
    }
  return retval;
}
//-----------------------------------------------------------------------
int
tripleCartesianTest()
{
  // Create a triple-Cartesian coordinate system.
  // This creates a 3 x 3 set of N/5 x N/5 patches in indexspace,
  // separated by buffers of width N/5.
  int boxsize = 16;
  RealVect dx = 1.0*RealVect::Unit;
  Box domainBox(IntVect::Zero, 5*boxsize*IntVect::Unit-1);
  ProblemDomain coarDomain(domainBox);
  ProblemDomain fineDomain = coarDomain;
  fineDomain.refine(2);
  TripleCartesianCS* coordSysDC = new TripleCartesianCS();
  coordSysDC->define(coarDomain, dx);
  RefCountedPtr<MultiBlockCoordSys > coordSys(dynamic_cast<MultiBlockCoordSys* > (coordSysDC));

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
      int inextcase = getTCSFineBoxes(fineBoxes, boxsize, icase);
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);

      // Create the flux register.
      NewMultiBlockFluxRegister fluxReg(coordSys, fineGrids, coarGrids, 2, 1, true);
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
      fluxReg.poutCoarseRegisters();
      fluxReg.poutFineRegisters();
      fluxReg.reflux(Ucoar, interv, 1.0);

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
int getCSTFineBoxes(Vector<Box>& a_fineBoxes,
                    const Vector<Box>& a_coarBoxes,
                    const int  & a_boxsize,
                    const int  & a_case,
                    const int & a_nlayers = 0)
{
  if (a_case < a_coarBoxes.size())
    {
      pout() << "refining and shrinking the box low end of box " << a_case << endl;
      Box oneBox  = a_coarBoxes[a_case];
      oneBox.refine(IntVect(D_DECL6(2,2,1,1,1,1)));
      for (int idir = 0; idir < 2; idir++) // Don't shrink dim 2 for aniso
        {
          oneBox.growHi(idir, -a_boxsize);
        }
      a_fineBoxes.resize(1, oneBox);
      return a_case+1;
    }
  else if (a_case <= (2*a_coarBoxes.size()-1))
    {
      int ibox = a_case - a_coarBoxes.size();
      pout() << "refining and shrinking  the box high end of box " << ibox << endl;
      Box oneBox  = a_coarBoxes[ibox];
      oneBox.refine(IntVect(D_DECL6(2,2,1,1,1,1)));
      for (int idir = 0; idir < 2; idir++) // Don't shrink dim 2 for aniso
        {
          oneBox.growLo(idir, -a_boxsize);
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
          oneBox.refine(IntVect(D_DECL6(2,2,1,1,1,1)));
          for (int idir = 0; idir < 2; idir++) // Don't shrink dim 2 for aniso
            {
              oneBox.growLo(idir, -a_boxsize);
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
          oneBox.refine(IntVect(D_DECL6(2,2,1,1,1,1)));
          for (int idir = 0; idir < 2; idir++) // Don't shrink dim 2 for aniso
            {
              oneBox.growHi(idir, -a_boxsize);
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
#if CH_SPACEDIM == 2
int
cubedSphereTest()
{
  int boxsize = 16;
  RealVect dx = 1.0*RealVect::Unit;
  // Create a cubed sphere coordinate system with six  blocks.
  CubedSphere2DCS* coordSysCS2DCS = new CubedSphere2DCS(boxsize);
  IntVect ivhi = (boxsize-1)*IntVect::Unit;
  ivhi[0] = (11*boxsize-1);
  Box domainBox(IntVect::Zero, ivhi);
  ProblemDomain coarDomain(domainBox);
  ProblemDomain fineDomain = coarDomain;
  fineDomain.refine(2);
  coordSysCS2DCS->define(coarDomain, dx);
  RefCountedPtr<MultiBlockCoordSys > coordSys(dynamic_cast<MultiBlockCoordSys* >(coordSysCS2DCS));

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
      int inextcase = getCSTFineBoxes(fineBoxes, coarBoxes, boxsize, icase);
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);

      // Create the flux register.
      NewMultiBlockFluxRegister fluxReg(coordSys, fineGrids, coarGrids, 2, 1, true);
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
      fluxReg.poutCoarseRegisters();
      fluxReg.poutFineRegisters();
      fluxReg.reflux(Ucoar, interv, 1.0);

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
                  pout() << "cubed sphere test failed for case = " << icase << ", cell = " <<  Ig  << endl;
                  return -1;
                }
            }
        }
      icase = inextcase;
    }

  // Git.
  return 0;
}
#endif

//-----------------------------------------------------------------------
#if CH_SPACEDIM >= 2
int
anisoCubedSphereShellTest()
{
  int boxsize = 16;
  int nlayers = (SpaceDim == 3) ? 8 : 1;
  RealVect dx = 1.0*RealVect::Unit;
  // Create a cubed sphere coordinate system with six  blocks.
  CubedSphereShellCS* css = new CubedSphereShellCS(boxsize,nlayers);
  // Cubed sphere shell has a specified size of index space 
  IntVect size(D_DECL6(11*boxsize, boxsize, nlayers, 1, 1, 1));
  Box domainBox(IntVect::Zero, size-IntVect::Unit);
  ProblemDomain coarDomain(domainBox);
  ProblemDomain fineDomain = coarDomain;
  fineDomain.refine(IntVect(D_DECL6(2,2,1,1,1,1)));
  css->define(coarDomain, dx);
  RefCountedPtr<MultiBlockCoordSys > coordSys(
      dynamic_cast<MultiBlockCoordSys* >(css));

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
      int inextcase = getCSTFineBoxes(fineBoxes, coarBoxes, boxsize, icase, nlayers);
      LoadBalance(fineProcs, fineBoxes);
      DisjointBoxLayout fineGrids(fineBoxes, fineProcs, fineDomain);

      // Create the flux register.
      NewMultiBlockFluxRegister fluxReg(coordSys, fineGrids, coarGrids, 2, 1, true);
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
      fluxReg.poutCoarseRegisters();
      fluxReg.poutFineRegisters();
      fluxReg.reflux(Ucoar, interv, 1.0);

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
                  pout() << "aniso sphere test failed for case = " << icase << ", cell = " <<  Ig  << endl;
                  return -1;
                }
            }
        }
      icase = inextcase;
    }

  // Git.
  return 0;
}
#endif

//-----------------------------------------------------------------------

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
    int status = 0;
    parseTestOptions( argc ,argv ) ;
    if ( verbose )
      pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

    status = singleBlockTest();
    if ( status == 0 )
    {
      cout << indent << pgmname << ": singleBlockTest passed." << endl ;
    }
    else
    {
      cout << indent << pgmname << ": singleBlockTest failed with return code " << status << endl ;
      overallStatus += 1;
    }

#if CH_SPACEDIM >= 2
    status = doubleCartesianTest();
    if ( status == 0 )
    {
      cout << indent << pgmname << ": doubleCartesianTest passed." << endl ;
    }
    else
    {
      cout << indent << pgmname << ": doubleCartesianTest failed with return code " << status << endl ;
      overallStatus += 1;
    }


    status = tripleCartesianTest();
    if ( status == 0 )
    {
      cout << indent << pgmname << ": tripleCartesianTest passed." << endl ;
    }
    else
    {
      cout << indent << pgmname << ": tripleCartesianTest failed with return code " << status << endl ;
      overallStatus += 1;
    }

#if CH_SPACEDIM == 2
    status = cubedSphereTest();
    if ( status == 0 )
    {
      cout << indent << pgmname << ": cubedSphereTest passed." << endl ;
    }
    else
    {
      cout << indent << pgmname << ": cubedSphereTest failed with return code " << status << endl ;
      overallStatus += 1;
    }
#endif

    status = anisoCubedSphereShellTest();
    if ( status == 0 )
    {
      cout << indent << pgmname << ": anisoCubedSphereShellTest passed." << endl ;
    }
    else
    {
      cout << indent << pgmname << 
        ": anisoCubedSphereShellTest failed with return code " << status << endl ;
      overallStatus += 1;
    }
#endif
  }

#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return overallStatus;
}
//-----------------------------------------------------------------------
