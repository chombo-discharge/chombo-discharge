#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereIF.H"
#include "PlaneIF.H"

#include "MFIndexSpace.H"
#include "MFCellFAB.H"
#include "MFRemapper.H"
#include "MFAliasFactory.H"
#include "MFPoissonOp.H"

#include "EBISLayout.H"
#include "EBAMRIO.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "CH_Attach.H"
#include "EBAMRPoissonOp.H"
#include "MFStencil.H"

#include "UsingNamespace.H"

int readGeometryInfo(Box& a_domain,
                     Real& a_dx,
                     RealVect& a_origin,
                     RealVect& a_center,
                     Real& a_radius)
{

  // parse input file
  ParmParse pp;

  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);

  IntVect lo = IntVect::Zero;
  IntVect hi;

  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  a_dx = (prob_hi-prob_lo[0])/n_cell[0];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_origin[idir] = prob_lo[idir];
    }

  pp.get("radius",a_radius);

  // ParmParse doesn't get RealVects, so work-around with Vector<Real>
  Vector<Real> vectorCenter;
  pp.getarr("center",vectorCenter,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_center[idir] = vectorCenter[idir];
    }

  return 0;
} // end read from file

int makeGeometry(MFIndexSpace& mfIndexSpace,
                 const Box& a_domain,
                 const Real& a_dx,
                 const RealVect& a_origin,
                 const RealVect& a_center,
                 const Real& a_radius)
{

  int eekflag = 0;
  Vector<GeometryService*> geometry(2, NULL);

  RealVect normal = BASISV(0);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  SphereIF is(a_radius, a_center, true);
  //PlaneIF   is(normal, a_center, true);
  GeometryShop workshop0(is,0,vectDx);
  workshop0.m_phase = 0;
  geometry[0] = &workshop0;

  SphereIF os(a_radius, a_center, false);
  //PlaneIF   os(normal, a_center, false);
  GeometryShop workshop1(os,0,vectDx);
  workshop1.m_phase=1;
  geometry[1] = &workshop1;

  mfIndexSpace.define(a_domain,a_origin,a_dx,geometry);

  return eekflag;
}
/***/
int getStencils(LayoutData< Vector<MFStencil::agg_t> >&  a_stencil,
                const Vector<EBISLayout>&                a_ebisl,
                const ProblemDomain&                     a_domain,
                const Real&                              a_dx)
{
  int eekflag = 0;
  RealVect rvdx = a_dx*RealVect::Unit;
  const BoxLayout& dbl = a_stencil.boxLayout();
  IntVectSet cfivs; //I magically know there is no cfivs
  for (DataIterator dit = a_stencil.dataIterator(); dit.ok(); ++dit)
    {
      //We want only irregular cells that are in both fluids.
      //because we are going to be adding some off-fluid  terms just for grins.
      IntVectSet ivsIrregIntersect = a_ebisl[0][dit()].getIrregIVS(dbl[dit()]);
      ivsIrregIntersect &= a_ebisl[1][dit()].getIrregIVS(dbl[dit()]);
      a_stencil[dit()].resize(0);
      for (int ifluid = 0; ifluid < 2; ifluid++)
        {
          int iother = 0;
          if (ifluid == 0) iother = 1;
          for (VoFIterator vofit(ivsIrregIntersect, a_ebisl[ifluid][dit()].getEBGraph()); vofit.ok(); ++vofit)
            {
              MFStencil::agg_t thisAgg;
              thisAgg.destVoF   = vofit();
              thisAgg.destFluid = ifluid;
              EBAMRPoissonOp::getDivFStencil(thisAgg.stenFluid[ifluid], vofit(),a_ebisl[ifluid][dit()], cfivs, rvdx);
              //just to test that off fluid terms work...
              int numfaces = a_ebisl[ifluid][dit()].numFacePhase(vofit());
              for (int iface = 0; iface < numfaces; iface++)
                {
                  VolIndex otherVoF = a_ebisl[ifluid][dit()].faceIndex(vofit(), iface);
                  thisAgg.stenFluid[iother].add(otherVoF, 1.0/a_dx);
                }
              a_stencil[dit()].push_back(thisAgg);
            }
        }
    }
  return eekflag;
}
/***/
Real getValue(const IntVect& a_iv, Real a_dx)
{
  RealVect xval(a_iv);
  xval += 0.5*(RealVect::Unit);
  xval *= a_dx;
  Real value = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      value += xval[idir]*xval[idir] + xval[idir];
    }
  return value;
}
/***/
int fillInputData(LevelData<MFCellFAB>&       a_data,
                  const DisjointBoxLayout&    a_dbl,
                  const Vector<EBISLayout>&   a_ebisl,
                  const ProblemDomain&        a_domain,
                  const Real&                 a_dx)
{
  int eekflag = 0;
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivs(a_dbl.get(dit()));
      for (int ifluid = 0; ifluid < 2; ifluid++)
        {
          EBCellFAB& fab = a_data[dit()].getPhase(ifluid);
          for (VoFIterator vofit(ivs, a_ebisl[ifluid][dit()].getEBGraph()); vofit.ok(); ++vofit)
            {
              fab(vofit(), 0) = getValue(vofit().gridIndex(), a_dx);
            }
        }
    }
  a_data.exchange();
  return eekflag;
}
/***/
int compareFastAndSlowEval(const LayoutData< Vector<MFStencil::agg_t> >&  a_stencil,
                           const DisjointBoxLayout &                      a_dbl,
                           const Vector<EBISLayout>&                      a_ebisl,
                           const ProblemDomain&                           a_domain,
                           const Real&                                    a_dx)
{
  int eekflag = 0;
  Vector<int> ncomp(2, 1); //each fluid gets one comp

  MFCellFactory mfcffact((Vector<EBISLayout>&)a_ebisl, ncomp);
  IntVect ghostIVPhi = 4*IntVect::Unit;
  IntVect ghostIVLph =   IntVect::Zero;
  LevelData<MFCellFAB>  inputDat(a_dbl, 1, ghostIVPhi, mfcffact);
  LevelData<MFCellFAB> outputDat(a_dbl, 1, ghostIVLph, mfcffact);

  eekflag = fillInputData(inputDat, a_dbl, a_ebisl, a_domain, a_dx);
  if (eekflag != 0)
    {
      pout() << "problem in fillinputdata " << endl;
    }
  for (DataIterator dit = a_stencil.dataIterator(); dit.ok(); ++dit)
    {
      outputDat[dit()].setVal(0.);
      Vector<EBISBox> vectEBISBox;
      for (int ifluid = 0; ifluid < 2; ifluid++)
        {
          vectEBISBox.push_back(a_ebisl[ifluid][dit()]);
        }
      MFStencil aggSten(a_stencil[dit()], a_dbl[dit()], vectEBISBox,
                        ghostIVLph, ghostIVPhi, 0);
      //this evaluates stuff the fast way
      aggSten.apply(outputDat[dit()], inputDat[dit()], false);

      //now let us compare the answer to the slow way
      for (int isten = 0; isten < a_stencil[dit()].size(); isten++)
        {
          const MFStencil::agg_t& pointSten = a_stencil[dit()][isten];
          const EBCellFAB& inputfab  =  inputDat[dit()].getPhase(pointSten.destFluid);
          const EBCellFAB& outputfab = outputDat[dit()].getPhase(pointSten.destFluid);
          Real slowAns = 0;
          for (int ifluid = 0; ifluid < 2; ifluid++)
            {
              const VoFStencil& vofsten = pointSten.stenFluid[ifluid];
              for (int isten = 0; isten < vofsten.size(); isten++)
                {
                  const Real& weight  = vofsten.weight(isten);
                  const VolIndex& vof = vofsten.vof(isten);
                  slowAns += weight*inputfab(vof, 0);
                }
            }
          Real fastAns = outputfab(pointSten.destVoF, 0);
          Real tol = 1.0e-4;
          if (Abs(fastAns-slowAns) > tol)
            {
              pout() << "answers do not match at vof " << pointSten.destVoF << endl;
              eekflag = -7;
              return eekflag;
            }
        }
    }
  return eekflag;
}


/***/
int testMFStencil(MFIndexSpace& a_mfIndexSpace,
                  const Box&  a_domain,
                  const Real& a_dx)
{
  int eekflag = 0;
  ParmParse pp;
  int maxBoxSize;
  pp.get("maxboxsize", maxBoxSize);
  Vector<Box> vbox;
  Vector<int> proc;

  domainSplit(a_domain, vbox, maxBoxSize, maxBoxSize);
  LoadBalance(proc, vbox);
  DisjointBoxLayout dbl(vbox, proc);

  Vector<EBISLayout> ebisl(2);
  a_mfIndexSpace.fillEBISLayout(ebisl, dbl, a_domain, 4);

  LayoutData< Vector<MFStencil::agg_t> >  stencil(dbl);
  eekflag = getStencils(stencil, ebisl, a_domain, a_dx);
  if (eekflag != 0)
    {
      pout() << "problem in getStencils" << endl;
      return eekflag;
    }

  eekflag = compareFastAndSlowEval(stencil, dbl, ebisl, a_domain, a_dx);
  if (eekflag != 0)
    {
      pout() << "problem in compareFastAndSlowEval" << endl;
      return eekflag;
    }

  return eekflag;
}
/***/
int testStuff()
{
  ParmParse pp;
  RealVect center;
  Real radius;
  RealVect origin;
  Real dx;
  Box domain;


  readGeometryInfo(domain,
                   dx,
                   origin,
                   center,
                   radius);

  int nsteps = 5;
  int eekflag = 0;
  //jitter geometry around and run a test.
  for (int istep = 0; istep < nsteps; istep++)
    {
      pout() << "testing step " << istep << endl;
      pout() << "center  = " << center << endl;
      pout() << "radius  = " << radius << endl;
      MFIndexSpace mfIndexSpace;
      eekflag = makeGeometry(mfIndexSpace,
                             domain,
                             dx,
                             origin,
                             center,
                             radius);

      if (eekflag != 0)
        {
          pout() << "problem in makegeom" << endl;
          return eekflag;
        }

      eekflag = testMFStencil(mfIndexSpace, domain, dx);
      if (eekflag != 0)
        {
          pout() << "problem in testmfstencil" << endl;
          return eekflag;
        }

      center[0]-= dx/3.0;
      center[1]-= dx/2.0;
      radius += dx/6.0;
    }
  return eekflag;
}
char iter_str[80];

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag = 0;
  // begin forever present scoping trick
  {
    const char* in_file = "sphere.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    int eekflag = testStuff();
    if (eekflag != 0)
      {
        pout() << "mfstencil test failed with error code" << eekflag << endl;
      }
    else
      {
        pout() << "mfstencil test passed " << endl;
      }
  } // end scoping trick


#ifdef CH_MPI
    MPI_Finalize();
#endif

    return eekflag;
}



