#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
using std::string;
#include "MBAggStencil.H"
#include "FArrayBox.H"
#include "CONSTANTS.H"
#include "RealVect.H"
#include "Box.H"
#include "IVSFAB.H"
#include "CH_Attach.H"
#include "CH_Timer.H"
#include "AMRIO.H"
#include "FABView.H"
#include "UsingNamespace.H"
#if (CH_SPACEDIM == 2) || (CH_SPACEDIM==3)
///
void  
getSlowStencil(IVSFAB<MBStencil>  & a_slowStenc, 
               const IntVectSet   & a_ivs, 
               const Real         & a_dx,
               const int          & a_numBlocks)
{
  CH_TIME("getslowstencil");
  int whichblock = 0;
  Real weightcen =  -4.0/a_dx/a_dx;
  Real weightout =   1.0/a_dx/a_dx;
  for(IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      RefCountedPtr<Vector<MBStencilElement> > ptsten(new Vector<MBStencilElement>());
      //deliberate non-use of SpaceDim (want 2D sten in this case)
      for(int idir = 0; idir < 2; idir++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
            {
              int isign = sign(sit());
              IntVect ivshift = iv + isign*BASISV(idir);

              MBStencilElement elemout(ivshift, whichblock, weightout);
              ptsten->push_back(elemout);

              // move blocks around because that is what we are testing
              whichblock++;
              if(whichblock >= a_numBlocks)
                {
                  whichblock = 0;
                }
            }
        }
      MBStencilElement elemcen(iv, whichblock, weightcen);
      ptsten->push_back(elemcen);
      
      a_slowStenc(iv, 0).define(ptsten);
    }
}
///
void
fillSrcData(FArrayBox  & a_srcData,
            const int  & a_block,
            const Real & a_dx)
{
  CH_TIME("fillsrcdata");
  for(BoxIterator bit(a_srcData.box()); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      RealVect loc;
      Real value = 1;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          loc[idir] = a_dx*(iv[idir] + 0.5);
          value *= sin(PI*loc[idir]);
        }
      for(int icomp = 0; icomp < a_srcData.nComp(); icomp++)
        {
          a_srcData(iv, icomp) = value + Real(icomp + a_block);
        }
    }
}
///
void
getSlowAns(FArrayBox                & a_slowAns, 
           const IVSFAB<MBStencil>  & a_slowStenc, 
           const Vector<FArrayBox*> & a_srcData)
{
  CH_TIME("getslowans");
  Box ansbox = a_slowAns.box();
  int zstart = 0;
  int zend   = 0;
  if(SpaceDim ==3)
    {
      zstart  = ansbox.smallEnd()[2];
      zend    = ansbox.bigEnd()  [2];
    }
  for(int zind = zstart; zind <= zend; zind++)
    {
      for(IVSIterator ivsit(a_slowStenc.getIVS()); ivsit.ok(); ++ivsit)
        {
          const IntVect& iv = ivsit();
          const MBStencil& stenc = a_slowStenc(iv, 0);
          for(int icomp = 0; icomp < a_slowAns.nComp(); icomp++)
            {
              Real value  = 0;
              for(int isten = 0; isten < stenc.size(); isten++)
                {
                  IntVect  stencell = stenc[isten].cell  ();
                  if(SpaceDim == 3)
                    {
                      stencell[2] = zind;
                    }
              
                  const Real      & weight = stenc[isten].weight();
                  const Real      & block  = stenc[isten].block ();
                  const FArrayBox & srcdat = *(a_srcData[block]);
                  const Real      & srcval = srcdat(stencell, icomp);
                  value += weight*srcval;
                }
              IntVect dstcell = iv;
              if(SpaceDim == 3)
                {
                  dstcell[2] = zind;
                }
          
              a_slowAns(dstcell, icomp) = value;
            }
        }
    }
}
///
void
getFastAns(FArrayBox                & a_fastAns, 
           const MBAggStencil       & a_fastStenc, 
           const Vector<FArrayBox*> & a_srcData)
{
  CH_TIME("getfastans");

  int isrc = 0; int idst = 0; int inco = a_fastAns.nComp();
  //the fast is to not do incrementonly (set to zero first)
  a_fastStenc.apply(a_fastAns, a_srcData, isrc, idst, inco, false);
}
///
int 
compareAns(const FArrayBox & a_slowAns, 
           const FArrayBox & a_fastAns, 
           const Box       & a_interBox)
{

#ifdef CH_USE_FLOAT
  Real tolerance = 1.0e-2;
#else
  Real tolerance = 1.0e-6;
#endif

  for(BoxIterator bit(a_interBox); bit.ok(); ++bit)
    {
      for(int icomp  = 0; icomp < a_slowAns.nComp(); icomp++)
        {
          Real fastans = a_fastAns(bit(), icomp);
          Real slowans = a_slowAns(bit(), icomp);
          if(Abs(fastans-slowans) > tolerance)
            {
              pout() << "answers different at iv = " << bit() << ", comp = " << icomp<< endl;
              pout() << "fast ans = " << fastans << endl;
              pout() << "slow ans = " << slowans << endl;
              return -1;
            }
        }
    }

  return 0;
}
///
int
testMBAggStencil()
{
  CH_TIME("testmbaggstencil");
  int ncomp = 1;
  int nx    = 32;
  int nz    = 8;
  int numBlocks = 6;
  IntVect ivlo =        IntVect::Zero;
  IntVect ivhi = (nx-1)*IntVect::Unit;
  Real dx = 1./nx;
  if(SpaceDim == 3) 
    {
      ivhi[2] = nz-1;
    }

  Box interBox(ivlo, ivhi);
  Box ghostBox = interBox;
  //not growing in z deliberately
  for(int idir = 0; idir < 2; idir++)
    {
      ghostBox.grow(idir, 1);
    }

  IntVectSet ivssten(interBox);
  RefCountedPtr<IVSFAB<MBStencil> > slowStenc(new IVSFAB<MBStencil>(ivssten, 1));
  getSlowStencil(*slowStenc, ivssten, dx, numBlocks);
  int stencilDim = 2;  //using 2D stencils in both 2 and 3 dimensions
  RefCountedPtr<MBAggStencil> fastStenc(new MBAggStencil(slowStenc, ivssten, ghostBox, interBox, stencilDim));

  Vector<FArrayBox*> srcData(numBlocks, NULL);
  for(int iblock = 0; iblock < numBlocks; iblock++)
    {
      srcData[iblock] = new FArrayBox(ghostBox, ncomp);      
      fillSrcData(*srcData[iblock], iblock, dx);
    }

  FArrayBox slowAns(interBox, ncomp);
  FArrayBox fastAns(interBox, ncomp);
  slowAns.setVal(0.);
  fastAns.setVal(0.);

  getSlowAns(slowAns, *slowStenc, srcData);
  getFastAns(fastAns, *fastStenc, srcData);

  for(int iblock = 0; iblock < numBlocks; iblock++)
    {
      delete srcData[iblock];
    }

  int retval = compareAns(slowAns, fastAns, interBox);
  return retval;
}

int main(int argc, char* argv[])
{
  int status = 0; // number of errors detected.
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  //scoping trick
  {
   status =  testMBAggStencil();
  if(status == 0)
    {
      pout() << "testMBAggStencil passed"  << endl;
    }
  else
    {
      pout() << "testMBAggStencil failed with error code "  << status << endl;
    }
  } // end scoping trick
#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
  return status;
}
#else
//spacedim != 1 or 2
int main(int argc, char* argv[])
{
  int status = 0; // number of errors detected.
  return status;
}
#endif
