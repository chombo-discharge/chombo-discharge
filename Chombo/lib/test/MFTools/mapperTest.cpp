#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;
using std::endl;

#include "SphereIF.H"
#include "MFIndexSpace.H"
#include "MFCellFAB.H"
#include "MFRemapper.H"
#include "MFAliasFactory.H"

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
#include "UsingNamespace.H"

int readGeometryInfo(Box& a_domain,
                     Real& a_dx,
                     RealVect& a_origin,
                     RealVect& a_center,
                     Real& a_radius);

void makeHierarchy(Vector<DisjointBoxLayout>& dbl,
                   const ProblemDomain& baseDomain,
                   const IntVectSet&    baseTags);

void swapvector(Vector<LevelData<MFCellFAB>* >& left,
                Vector<LevelData<MFCellFAB>* >& right);

/***************/
// define a sphere EBIS.
/***************/
int makeGeometry(MFIndexSpace& mfIndexSpace,
                 const Box& domain,
                 const Real& dx,
                 const RealVect& origin,
                 const RealVect& center,
                 const Real& radius);

void makeLevelSet(Vector<LevelData<FArrayBox>* >& levelSet,
                 const Box& domain,
                 const Real& dx,
                 const RealVect& origin,
                 const RealVect& center,
                 const Real& radius);

/***************/
//make the corresponding layout
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& nghost);
/***************/
/***************/
/*
int checkEBISL(const EBISLayout& a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box& a_domain,
               const Real& dx,
               const RealVect& origin,
               const Real& radius,
               const RealVect& center,
               const bool& insideRegular);
*/

/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/

char iter_str[80];

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // registerDebugger();

  // begin forever present scoping trick
  {
    pout() << std::endl;

    Vector<std::string> names0(1, "phi0");
    Vector<std::string> names1(1, "phi1");
    Vector<std::string> names(2);
    names[0] = "phi0";
    names[1] = "phi1";

    Vector<int> refRatio(3,2);
    Vector<Real> coveredVal(1,7.0);

    const char* in_file = "mapperTest.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    RealVect center;
    Real radius;
    RealVect origin;
    RealVect aspect = RealVect::Unit;
    Real dx;
    Box domain;
    int eekflag = 0;
    MFIndexSpace a, b;
    MFIndexSpace *mfIndexSpace, *oldMFIndexSpace;
    mfIndexSpace = &a;
    oldMFIndexSpace = &b;

    LevelData<MFCellFAB> fine, med, coarse;
    LevelData<MFCellFAB> fineOld, medOld, coarseOld;
    LevelData<FArrayBox> lsCoarse, lsMedium, lsFine;

    Vector<LevelData<MFCellFAB>* > mfvector(3,NULL);
    mfvector[0]=&coarse;
    mfvector[1]=&med;
    mfvector[2]=&fine;
    Vector<LevelData<MFCellFAB>* > mfvector_old(3,NULL);
    mfvector_old[0]=&coarseOld;
    mfvector_old[1]=&medOld;
    mfvector_old[2]=&fineOld;
    Vector<LevelData<FArrayBox>* > levelset(3, NULL);
    levelset[0] = &lsCoarse;
    levelset[1] = &lsMedium;
    levelset[2] = &lsFine;

    readGeometryInfo(domain,
                     dx,
                     origin,
                     center,
                     radius);

    Box domainFine(domain), domainMedi, domainCoar;
    ProblemDomain pFine(domain);
    Real dxFine(dx), dxMedi, dxCoar;


    CH_assert(eekflag == 0);

    domainMedi = coarsen(domainFine, 2);
    domainCoar = coarsen(domainMedi, 2);
    dxMedi = 2.0*dxFine;
    dxCoar = 2.0*dxMedi;
    Vector<DisjointBoxLayout> grids(3);
    ProblemDomain baseDomain(domainCoar);
    ProblemDomain pMed(domainMedi);

    RealVect dxvec(D_DECL(dxCoar,dxCoar,dxCoar));

    MFRemapper remapper;

    // make data holders
    Vector<int> comps(2,1);
    int ghost = 1;

    int steps = 2;
    int step  = 0;


    Vector<LevelData<EBCellFAB>* > phaseA(3);
    Vector<LevelData<EBCellFAB>* > phaseB(3);

    for (int i=0; i<3; i++)
    {
      phaseA[i] = new LevelData<EBCellFAB>();
      phaseB[i] = new LevelData<EBCellFAB>();
    }

    while (step < steps)
    {
      eekflag = makeGeometry(*mfIndexSpace,
                             domain,
                             dx,
                             origin,
                             center,
                             radius);


      if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        // MayDay::Error("problem in makeGeometry");
      }

      if (step > 0)
      {
        // same grids, but new topology
        mfvector[2]->define(grids[2], 1, IntVect::Unit,
                            MFCellFactory(*mfIndexSpace, grids[2],
                                          domainFine,
                                          comps, ghost));

        mfvector[1]->define(grids[1], 1, IntVect::Unit,
                            MFCellFactory(*mfIndexSpace, grids[1],
                                          domainMedi,
                                          comps, ghost));

        mfvector[0]->define(grids[0], 1, IntVect::Unit,
                            MFCellFactory(*mfIndexSpace, grids[0],
                                          domainCoar,
                                          comps, ghost));
        DisjointBoxLayout d1 = mfvector_old[0]->getBoxes();
        CH_assert(d1 == grids[0]);
        for (int i = 0; i < 3; ++i)
        {
          remapper.remap(*oldMFIndexSpace,
                         *mfvector_old[i],
                         *mfIndexSpace,
                         *mfvector[i]);
        }

        swapvector(mfvector, mfvector_old);
      }

#if 0
//       if (step > 0)
//       {
//         for (int i = 0; i < 3; i++)
//         {
//           phase[i] = new LevelData<EBCellFAB>();
//         }

//         aliasMF(phase, 0, mfvector_old);
//         sprintf(iter_str, "phase0_old.%03d.%dd.hdf5",step, SpaceDim);
//         writeEBHDF5(iter_str, grids, phase, names0, domainCoar,
//                     dxCoar, 1, step, refRatio, 3, true, coveredVal);

//         aliasMF(phase, 1, mfvector_old);
//         sprintf(iter_str, "phase1_old.%03d.%dd.hdf5",step, SpaceDim);
//         writeEBHDF5(iter_str, grids, phase, names1, domainCoar,
//                     dxCoar, 1, step, refRatio, 3, true, coveredVal);
//       }
#endif

      // make grids
      IntVectSet tags = mfIndexSpace->interfaceRegion(2);
      // tags.grow(1);

      makeHierarchy(grids, baseDomain, tags);

      MFCellFactory fineFactory(*mfIndexSpace, grids[2],domainFine,
                                comps, ghost);
      mfvector[2]->define(grids[2], 1, IntVect::Unit, fineFactory);
      levelset[2]->define(grids[2], 1, IntVect::Unit);

      MFCellFactory  medFactory(*mfIndexSpace, grids[1], domainMedi,
                                comps, ghost);
      mfvector[1]->define(grids[1], 1, IntVect::Unit, medFactory);
      levelset[1]->define(grids[1], 1, IntVect::Unit);

      MFCellFactory coarseFactory(*mfIndexSpace, grids[0],domainCoar,
                                  comps, ghost);
      mfvector[0]->define(grids[0], 1, IntVect::Unit, coarseFactory);
      levelset[0]->define(grids[0], 1, IntVect::Unit);


      makeLevelSet(levelset,
                   domainCoar,
                   dxCoar,
                   origin,
                   center,
                   radius);


        //initialize data holders
      LevelData<EBCellFAB> phase0, phase1;
      for (int lev = 0; lev < 3; lev++)
        {
          aliasMF(phase0, 0, *mfvector[lev]);
          aliasMF(phase1, 1, *mfvector[lev]);
          DataIterator dit = phase0.dataIterator();

          for (dit.begin(); dit.ok(); ++dit)
            {
              phase0[dit].setVal(0.0);
              phase1[dit].setVal(1.0);
            }
        }

      if (step != 0)
      {

        // perform second type of remap.  constant MFIndexSpace, but
        // different grids.
        //
        // coarse grids in this case remain identical.
        mfvector_old[0] ->copyTo(mfvector[0]->interval(),*mfvector[0],
                                 mfvector_old[0]->interval());

        remapper.remap(*mfIndexSpace, baseDomain, pMed,
                       *mfvector_old[1], *mfvector[0],
                       2, 1, *mfvector[1]);
        remapper.remap(*mfIndexSpace, pMed, pFine,
                       *mfvector_old[2], *mfvector[1],
                       2, 1, *mfvector[2]);
      }


      aliasMF(phaseA, 0, mfvector);
#ifdef CH_USE_HDF5


      sprintf(iter_str, "phase0.%03d.%dd.hdf5",step, SpaceDim);
      //       writeEBHDF5(iter_str, grids, phaseA, names0, domainCoar,
      //                    dxCoar, 1, step, refRatio, 3, true, coveredVal);
      cerr << "\x1b[31m\x1b[48mwriting to " << iter_str << endl;


      HDF5Handle handle;
      createEBFile(handle, iter_str, 3, refRatio, domainCoar, origin, dxvec, aspect, ghost*IntVect::Unit);
      writeCellCenteredNames(handle, names0);
      writeCellCentered(handle, 0, phaseA[0]);
      writeCellCentered(handle, 1, phaseA[1]);
      writeCellCentered(handle, 2, phaseA[2]);

      handle.close();

      LevelData<EBCellFAB> readData;

      handle.open(iter_str, HDF5Handle::OPEN_RDONLY, "Chombo_global");

      readCellCentered(handle, 2, mfIndexSpace->EBIS(0), ghost, &readData);

      DataIterator dit_before =  phaseA[2]->dataIterator();
      DataIterator dit_after  =  readData.dataIterator();
      for (; dit_before.ok(); ++dit_before, ++dit_after)
        {
          EBCellFAB& ef = readData[dit_after];
          const EBCellFAB& before = phaseA[2]->operator[](dit_before);
          //pout()<<phase[2]->box(dit_before())<<"  "<<readData.box(dit_after())<<std::endl;
          CH_assert(phaseA[2]->box(dit_before()) == readData.box(dit_after()));
          ef -= before;
          Real norm = ef.norm(0,0,1);
          if (norm != 0)
          {
            MayDay::Error("Error `reading LevelData<EBCellFAB>");
          }
        }

      handle.close();

      aliasMF(phaseB, 1, mfvector);

      sprintf(iter_str, "phase1.%03d.%dd.hdf5",step, SpaceDim);
      //       writeEBHDF5(iter_str, grids, phaseB, names1, domainCoar,
      //                    dxCoar, 1, step, refRatio, 3, true, coveredVal);
      cerr << "\x1b[33m\x1b[48mwriting to " << iter_str << endl;

      //openEBFile(handle, iter_str, 3, refRatio, domainCoar, dxCoar, ghost*IntVect::Unit);
      createEBFile(handle, iter_str, 3, refRatio, domainCoar, origin, dxvec, aspect, IntVect::Zero);
      writeCellCenteredNames(handle, names1);
      writeCellCentered(handle, 0, phaseB[0]);
      writeCellCentered(handle, 1, phaseB[1]);
      writeCellCentered(handle, 2, phaseB[2]);

      handle.close();

      handle.open(iter_str, HDF5Handle::OPEN_RDONLY, "Chombo_global");

      readCellCentered(handle, 1, mfIndexSpace->EBIS(1), ghost, &readData);

      handle.close();

      dit_before =  phaseB[1]->dataIterator();
      dit_after  =  readData.dataIterator();
      for (; dit_before.ok(); ++dit_before, ++dit_after)
        {
          CH_assert(phaseB[1]->box(dit_before()) == readData.box(dit_after()));
          readData[dit_after] -= phaseB[1]->operator[](dit_before);
          Real norm = readData[dit_after].norm(0,0,1);
          if (norm != 0)
          {
            MayDay::Error("Error reading LevelData<EBCellFAB>");
          }
        }

      sprintf(iter_str, "mPhase.%03d.%dd.hdf5",step, SpaceDim);
      // put writeMFAMR function call here.
 //      writeEBHDF5(iter_str, grids, &phaseA, &phaseB, &levelset, names, domainCoar,
//                   dxCoar, 1, step, refRatio, 3, true, coveredVal, IntVect::Unit);
#endif

      step++;

      center[0]+= dx/3.0;
      center[1]+= dx/2.0;
      radius += dx/5.0;

      {
        MFIndexSpace* swap = oldMFIndexSpace;
        oldMFIndexSpace = mfIndexSpace;
        mfIndexSpace = swap;
      }

      swapvector(mfvector, mfvector_old);
      pout()<<step<<std::endl;
    }


    for (int i=0; i<3; i++)
      {
        delete phaseA[i];
        delete phaseB[i];
      }
    pout() << "\n mapper test passed" << endl;
    } // end scoping trick


#ifdef CH_MPI
    MPI_Finalize();
#endif

    return 0;
  }

  void swapvector(Vector<LevelData<MFCellFAB>* >& left,
                Vector<LevelData<MFCellFAB>* >& right)

{

  LevelData<MFCellFAB>* swap;
  for (int i = 0; i < left.size(); ++i)
  {
    swap     = left[i];
    left[i]  = right[i];
    right[i] = swap;
  }
}

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

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  SphereIF   inside(a_radius, a_center, true);
  GeometryShop workshop0(inside,0,vectDx);
  workshop0.m_phase = 0;
  geometry[0] = &workshop0;

  SphereIF  outside(a_radius, a_center, false);
  GeometryShop workshop1(outside,0,vectDx);
  workshop1.m_phase=1;
  geometry[1] = &workshop1;

  mfIndexSpace.define(a_domain,a_origin,a_dx,geometry);

  return eekflag;
}

void makeLevelSet(Vector<LevelData<FArrayBox>* >& levelSet,
                 const Box& domain,
                 const Real& dx,
                 const RealVect& origin,
                 const RealVect& center,
                 const Real& radius)
{
  Real d = dx;
  Box  box = domain;
  RealVect loc;
  RealVect newCenter = center;
  newCenter -=origin;
  for (int i=0; i<levelSet.size(); i++)
    {
      LevelData<FArrayBox>& data = *(levelSet[i]);
      for (DataIterator dit = data.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& fab = data[dit];
          ForAllX(Real, fab)
            {
              D_EXPR(loc[0]=iR+nR-nR, loc[1]=jR, loc[2]=kR);
              loc+=0.5;
              loc*=d;
              loc-=newCenter;
              Real r = D_TERM(loc[0]*loc[0],+loc[1]*loc[1], +loc[2]*loc[2]);
              r = sqrt(r);
              fabR = radius - r;
            } EndFor ;
        }

      d/=2;
      box.refine(2);
    }


}
void makeHierarchy(Vector<DisjointBoxLayout>& dbl,
                   const ProblemDomain& baseDomain,
                   const IntVectSet&    baseTags)
{

  dbl.resize(3);
  Vector<int> refRatios(2,2);
  Real fillRatio     = 0.85;
  int blockingFactor = 4;
  int bufferSize     = 2;
  int maxSize        = 32;
  BRMeshRefine regridder(baseDomain, refRatios, fillRatio, blockingFactor,
                         bufferSize, maxSize);

  Vector<Vector<Box> > oldGrids(3,1), newGrids(3);
  oldGrids[0][0]=baseDomain.domainBox();
  oldGrids[1][0]=coarsen(oldGrids[0][0], 2);
  oldGrids[2][0]=coarsen(oldGrids[1][0], 2);

  regridder.regrid(newGrids, baseTags, 0, 1, oldGrids);

  Vector<int> procs;
  for (int i = 0; i < 3; i++)
  {
    newGrids[i].sort();
    LoadBalance(procs, newGrids[i]);
    // pout()<<newGrids[i]<<" : "<<procs<<std::endl;
    dbl[i] = DisjointBoxLayout(newGrids[i], procs);
  }
}
