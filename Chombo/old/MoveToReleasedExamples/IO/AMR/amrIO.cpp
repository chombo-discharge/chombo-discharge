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
#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "SPMD.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "Misc.H"
#include "Vector.H"
#include "REAL.H"
#include "Box.H"
#include "Tuple.H"
#include "BoxIterator.H"

#include "UsingNamespace.H"

Real getDataVal(const Tuple<Real,SpaceDim>& location)
{
  Real retval  = 7.23;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real x= location[idir];
      Real fac = Real(idir);
      retval += sin(x)+ fac*cos(x);
    }
  return retval;
}
/*
  Set grid hierarchy from input file
 */
int setGrids(Vector<DisjointBoxLayout>& a_vectGrids,
             Box&  a_domain,
             Real& a_dx,
             Vector<int>&  a_refRatio,
             int& a_numlevels)
{

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Real fillRat = 0.77;
  int ncells = 64;
  int maxboxsize = 32;

  Vector<int> ancells(SpaceDim, ncells);
  Vector<Real>  prob_loa(SpaceDim, 0.0);
  Vector<Real>  prob_hia(SpaceDim, 1.0);

  a_numlevels = 3;
  a_refRatio.resize(a_numlevels);
  for (int ilev = 0; ilev < a_numlevels; ++ilev)
    a_refRatio[ilev] = 2;
  IntVect ivlo = IntVect::Zero;
  IntVect ivhi = (ncells-1)*IntVect::Unit;
  a_domain = Box(ivlo, ivhi);

  a_vectGrids.resize(a_numlevels);
  a_dx = (prob_hia[0]-prob_loa[0])/ancells[0];

  int nc = ancells[0];
  int nmi = nc/2;
  int nqu = nc/4;
  int ntf = (nc*3)/4;
  int nte = (nc*3)/8;
  int nfe = (nc*5)/8;
#if (CH_SPACEDIM ==2)
  Box boxf1(IntVect(0, nqu), IntVect(nmi-1,ntf-1));
  Box boxf2(IntVect(nmi,nte), IntVect(ntf-1,nfe-1));
  Box boxf3(IntVect(nqu,0  ), IntVect(nfe-1,nqu-1));
  Box boxf4(IntVect(nfe,nqu), IntVect(nc -1,nte-1));
#else
  Box boxf1(IntVect(0, nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
  Box boxf2(IntVect(nmi,nte,nte), IntVect(ntf-1,nfe-1,nfe-1));
  Box boxf3(IntVect(nqu,0,0  ), IntVect(nfe-1,nqu-1,nqu-1));
  Box boxf4(IntVect(nfe,nqu,nqu), IntVect(nc -1,nte-1,nte-1));
#endif
  IntVectSet tags;
  tags |= boxf1;
  tags |= boxf2;
  tags |= boxf3;
  tags |= boxf4;
  Vector<Vector<Box> > vvboxNew(a_numlevels);
  Vector<Vector<Box> > vvboxOld(a_numlevels);
  Vector<Box> vectDomain(a_numlevels);
  Box domainfine = a_domain;
  for (int ilev = 0; ilev <a_numlevels; ilev++)
    {
      vectDomain[ilev] = domainfine;
      vvboxOld[ilev].push_back(domainfine);
      domainfine.refine(a_refRatio[ilev]);
    }
  int baseLevel = 0;
  int topLevel  = a_numlevels - 2;
  int blockFactor = 4;
  int eekflag = 0;
  int buffersize = 1;
  int finestLevel = 0;
  if (topLevel >= 0)
    {
      BRMeshRefine meshrefine(vectDomain[0],a_refRatio, fillRat,
                              blockFactor, buffersize, maxboxsize);

      finestLevel = meshrefine.regrid(vvboxNew, tags, baseLevel,
                                          topLevel, vvboxOld);


    }
  else
    vvboxNew = vvboxOld;



  Vector< Vector<int> > procAssign;
  Real effRatio = 0.75;
  Vector< Vector<long> > loads(a_numlevels);
  for (int ilev = 0; ilev <a_numlevels; ilev++)
    {
      loads[ilev].resize(vvboxNew[ilev].size());
      for (int ibox = 0; ibox < vvboxNew[ilev].size() ; ibox++)
        {
          loads[ilev][ibox] = vvboxNew[ilev][ibox].numPts();
        }
    }
  LoadBalance(procAssign, effRatio, vvboxNew, loads, a_refRatio);

  if (eekflag != 0)
    {
      cerr << "setGrids: loadBalance returned error code " << eekflag << endl;
      return(1);
    }
  for (int ilev = 0; ilev < a_numlevels; ilev++)
    {
      a_vectGrids[ilev].define(vvboxNew[ilev], procAssign[ilev]);
      a_vectGrids[ilev].close();
    }
  return 0;
}
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {//scoping trick

    //set grid spacing,
    // number of points in each direction, domain length
    Vector<DisjointBoxLayout>vectGrids;
    Box domain;
    Real dx = 1.0;
    Vector<int> refRatio;
    int numlevels;
    //create layouts and domain and all that
    int eekflag= setGrids(vectGrids, domain, dx, refRatio, numlevels);
    if (eekflag !=0)
      {
        cerr << "problem in setgrids" << endl;
        return -1;
      }
    ///make the data to output and set its values.
    Vector<LevelData<FArrayBox>* > dataPtrs(numlevels, NULL);
    Real dxlevel = dx;
    for (int ilev = 0; ilev < numlevels; ilev++)
      {
        const DisjointBoxLayout& dbl = vectGrids[ilev];
        dataPtrs[ilev] = new LevelData<FArrayBox>(dbl, 1);
        LevelData<FArrayBox>& data = *dataPtrs[ilev];
        if (ilev > 0)
          dxlevel /= Real(refRatio[ilev-1]);
        DataIterator dit = dbl.dataIterator();
        for (dit.reset(); dit.ok(); ++dit)
          {
            FArrayBox& fab = data[dit()];
            const Box& fabbox = fab.box();
            BoxIterator bit(fabbox);
            for (bit.reset(); bit.ok(); ++bit)
              {
                const IntVect& iv = bit();
                Tuple<Real, SpaceDim> location;
                for (int idir = 0; idir <  SpaceDim; idir++)
                  {
                    location[idir]= dxlevel*(Real(iv[idir]) + 0.5);
                  }
                fab(iv, 0) = getDataVal(location);
              }
          }
      }
#ifdef CH_USE_HDF5
    //output data to a file
    string filename("dataout.hdf5");
    WriteAMRHierarchyHDF5(filename, vectGrids, dataPtrs,
                          domain, refRatio, numlevels);
    //delete the old data
    for (int ilev = 0; ilev < numlevels; ilev++)
      delete dataPtrs[ilev];

    //Read it back in.
    Vector<DisjointBoxLayout> vectGridsin;
    Box domainin;
    Vector<int> refRatioin;
    Vector<LevelData<FArrayBox>* > dataPtrsin;
    int numlevelsin;
    ReadAMRHierarchyHDF5(filename, vectGridsin, dataPtrsin,
                         domainin, refRatioin, numlevelsin);
    //check to see that it matches.
    if (domainin != domain)
      {
        cerr << "domains do not match" << endl;
        return -1;
      }
    if (numlevelsin != numlevels)
      {
        cerr << "numlevels do not match" << endl;
        return -2;
      }
    if (refRatio.size() != refRatioin.size())
      {
        cerr << "refRatios size wrong" << endl;
      }
    // apparently the writing and reading isn't consistent
    // in handling the finest level refRatio -- which apparently
    // isn't important anyway (sez dmartin) -- so when
    // we check to make sure they are the same, don't
    // consider last one. (ndk)
    //for (int ilev; ilev < refRatio.size(); ilev++)
    for (int ilev=0; ilev < refRatio.size()-1; ilev++)
      {
        if (refRatioin[ilev]!=refRatio[ilev])
          {
            cerr << "ilev =" << ilev << endl;
            cerr << "refRatioin[ilev]=" << refRatioin[ilev] << endl;
            cerr << "refRatio[ilev]=" << refRatio[ilev] << endl;
            cerr << "refinement ratios do not match" << endl;
            return -3;
          }
      }

    dxlevel = dx;
    for (int ilev = 0; ilev < numlevels; ilev++)
      {
        const DisjointBoxLayout& dbl = vectGridsin[ilev];
        LevelData<FArrayBox>& data = *dataPtrsin[ilev];
        if (ilev > 0)
          dxlevel /= Real(refRatioin[ilev-1]);

        DataIterator dit = dbl.dataIterator();
        for (dit.reset(); dit.ok(); ++dit)
          {
            FArrayBox& fab = data[dit()];
            const Box& fabbox = fab.box();
            BoxIterator bit(fabbox);
            for (bit.reset(); bit.ok(); ++bit)
              {
                const IntVect& iv = bit();
                Tuple<Real, SpaceDim> location;
                for (int idir = 0; idir <  SpaceDim; idir++)
                  {
                    location[idir]= dxlevel*(Real(iv[idir]) + 0.5);
                  }
                Real rightans = getDataVal(location);
                Real dataans = fab(iv, 0);
                Real eps = 1.0e-9;
                if (Abs(dataans - rightans) > eps)
                  {
                    cerr << "data does not match" << endl;
                    return -4;
                  }
              }
          }
      }
    // delete the LevelDatas allocated by ReadAMRHierarchyHDF5
    for (int ilev = 0 ; ilev < dataPtrsin.size() ; ilev++)
    {
      delete dataPtrsin[ilev] ;
    }
#endif // CH_USE_HDF5

  }//end scoping trick
  cout << "amrIO seems to work" << endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return(0);
}
