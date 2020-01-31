#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TestUtilities.H"

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

#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "DirichletPoissonDomainBC.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testMFpoissonOp" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;


class InsideDirichletBC: public BaseBCValue
{
public:
  Real value (const RealVect &a_point,
              const RealVect &a_n,
              const Real&     a_time,
              const int&      a_comp
              ) const;
};

class OutsideDirichletBC: public BaseBCValue
{
public:
  Real value (const RealVect &a_point,
              const RealVect &a_n,
              const Real&     a_time,
              const int&      a_comp
              ) const;
};
class InsideRHS: public BaseBCValue
{
public:
  Real value (const RealVect &a_point,
              const RealVect &a_normal,
              const Real&     a_time,
              const int&      a_comp
              ) const;
};

class OutsideRHS: public BaseBCValue
{
public:
  Real value (const RealVect &a_point,
              const RealVect &a_normal,
              const Real&     a_time,
              const int&      a_comp
              ) const;
};

Real InsideDirichletBC::value (const RealVect &a_point,
                               const RealVect &a_normal,
                               const Real&     a_time,
                               const int&      a_comp) const
{
  return 3*a_point[0]*a_point[0]+2;
  //return a_point[0]*2.0/3.0 + a_point[1]*1.5 + 2;
  //return 2*a_point[0];
  //return 1;
}

Real OutsideDirichletBC::value (const RealVect &a_point,
                                const RealVect &a_normal,
                                const Real&     a_time,
                                const int&      a_comp
                                ) const
{
  return 3*a_point[0]*a_point[0];
  //return a_point[0]*2.0/3.0 + a_point[1]*1.5;
  // return 2*a_point[0];
  //return 0;
}

Real InsideRHS::value (const RealVect &a_point,
                       const RealVect &a_normal,
                       const Real&     a_time,
                       const int&      a_comp
                       ) const
{
  return 6;
  //return 0;

}

Real OutsideRHS::value (const RealVect &a_point,
                        const RealVect &a_normal,
                        const Real&     a_time,
                        const int&      a_comp
                        ) const
{
  return 6;
  //return 0;
}


int readGeometryInfo(Box& a_domain,
                     Real& a_dx,
                     RealVect& a_origin,
                     RealVect& a_center,
                     Real& a_radius);


void makeHierarchy(Vector<DisjointBoxLayout>& dbl,
                   const ProblemDomain& baseDomain,
                   const IntVectSet&    baseTags,
                   int a_maxSize);

void makeLevelSet(Vector<LevelData<FArrayBox>* >& levelSet,
                 const Box& domain,
                 const Real& dx,
                 const RealVect& origin,
                 const RealVect& center,
                 const Real& radius);

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
/***************/

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

  //registerDebugger();

  // begin forever present scoping trick
  {
    pout()<<std::endl;

    Vector<std::string> names0(1, "phi0");
    Vector<std::string> names1(1, "phi1");
    Vector<int> refRatio(3,2);
    Vector<Real> coveredVal(1,3.0);

    const char* in_file = "sphere.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    RealVect center;
    Real radius;
    RealVect origin;
    Real dx;
    Box domain;
    ProblemDomain pDomain(domain);
    int maxSize = 32;
    pp.get("maxboxsize",maxSize);

    int eekflag = 0;
    //MFIndexSpace a, b;
    //MFIndexSpace *mfIndexSpace, *oldMFIndexSpace;
    //mfIndexSpace = &a;
    //oldMFIndexSpace = &b;
    RefCountedPtr<MFIndexSpace> mfIndexSpace(new MFIndexSpace());


    LevelData<MFCellFAB> fine, med, coarse;

    LevelData<MFCellFAB> fineRHS, medRHS, coarseRHS;
    LevelData<MFCellFAB> fineResidual, mediumResidual, coarseResidual;

    LevelData<FArrayBox> lsCoarse, lsMedium, lsFine;
    Vector<LevelData<MFCellFAB>* > mfvector(3,NULL);
    Vector<LevelData<MFCellFAB>* > vresidual(3,NULL);
    Vector<LevelData<MFCellFAB>* > rhsvector(3,NULL);
    mfvector[0]=&coarse;
    mfvector[1]=&med;
    mfvector[2]=&fine;
    vresidual[0]=&coarseResidual;
    vresidual[1]=&mediumResidual;
    vresidual[2]=&fineResidual;
    rhsvector[0] = &coarseRHS;
    rhsvector[1] = &medRHS;
    rhsvector[2] = &fineRHS;

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
    Vector<RealVect> xVec(3, IntVect::Unit);
    xVec[0]*= dxCoar;
    xVec[1]*= dxMedi;
    xVec[2]*= dxFine;

    Vector<DisjointBoxLayout> grids(3);
    ProblemDomain baseDomain(domainCoar);
    ProblemDomain pMed(domainMedi);
    Vector<ProblemDomain> pd(3);
    pd[0] = baseDomain;
    pd[1] = pMed;
    pd[2] = ProblemDomain(domainFine);


    Vector<RefCountedPtr<BaseDomainBC> > bc;
    RefCountedPtr<DirichletPoissonDomainBC> insideBC( new DirichletPoissonDomainBC() );
    RefCountedPtr<DirichletPoissonDomainBC> outsideBC( new DirichletPoissonDomainBC() );
    bc.push_back(insideBC);
    bc.push_back(outsideBC);
//     insideBC->setValue(1.0);
//     outsideBC->setValue(0.0);
    RefCountedPtr<BaseBCValue> bcv( new InsideDirichletBC() );
    insideBC->setFunction(bcv);
    bcv = RefCountedPtr<BaseBCValue>(new OutsideDirichletBC());
    outsideBC->setFunction(bcv);



    //make data holders
    Vector<int> comps(2,1);
    int ghost = 3;

    int steps= 5;
    int step = 0;

    Vector<LevelData<EBCellFAB>* > phaseA(3, NULL), phaseB(3,NULL);
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


      //make grids
      //IntVectSet tags = mfIndexSpace->interfaceRegion(2);
      IntVectSet   tags(domainCoar);
      tags.grow(1);
      makeHierarchy(grids, baseDomain, tags, maxSize);






      MFCellFactory fineFactory(*mfIndexSpace, grids[2],domainFine,
                                comps, ghost);
      mfvector[2]->define(grids[2], 1, 3*IntVect::Unit, fineFactory);
      rhsvector[2]->define(grids[2], 1, IntVect::Zero, fineFactory);
      levelset[2]->define(grids[2], 1, IntVect::Unit);

      MFCellFactory  medFactory(*mfIndexSpace, grids[1], domainMedi,
                                comps, ghost);
      mfvector[1]->define(grids[1], 1, 3*IntVect::Unit, medFactory);
      rhsvector[1]->define(grids[1], 1, IntVect::Zero, medFactory);
      levelset[1]->define(grids[1], 1, IntVect::Unit);

      MFCellFactory coarseFactory(*mfIndexSpace, grids[0],domainCoar,
                                  comps, ghost);
      mfvector[0]->define(grids[0], 1, 3*IntVect::Unit, coarseFactory);
      rhsvector[0]->define(grids[0], 1, IntVect::Zero, coarseFactory);
      levelset[0]->define(grids[0], 1, IntVect::Unit);

      makeLevelSet(levelset,
                   domainCoar,
                   dxCoar,
                   origin,
                   center,
                   radius);

      //initialize data holders
      LevelData<EBCellFAB> phase0, phase1;
      for (int lev=0; lev<3; lev++)
        {
          aliasMF(phase0, 0, *mfvector[lev]);
          aliasMF(phase1, 1, *mfvector[lev]);
          setValue(phase0, InsideDirichletBC(), pd[lev].domainBox()  , xVec[lev], origin);
          setValue(phase1, OutsideDirichletBC(),pd[lev].domainBox()  , xVec[lev], origin);
          aliasMF(phase0, 0, *rhsvector[lev]);
          aliasMF(phase1, 1, *rhsvector[lev]);
          setValue(phase0, InsideRHS(), pd[lev].domainBox()  , xVec[lev], origin, true);
          setValue(phase1, OutsideRHS(),pd[lev].domainBox()  , xVec[lev], origin, true);
        }



      Vector<int> refRatio(3,2);


      for (int i=0; i<3; i++)
      {
        MFPoissonOp op;

        LevelData<MFCellFAB> &phi=*mfvector[i], &rhs=*rhsvector[i], &residual=*vresidual[i];
        Vector<Real> alpha(2, 0.);
        Vector<Real> beta(2, 1.0);
        DisjointBoxLayout dumdbl;
        op.define(*mfIndexSpace, 1, grids[i], dumdbl, false, true, dumdbl, dumdbl,
                  xVec[i], 2, 2,pd[i], bc, 3*IntVect::Unit, IntVect::Zero, false, false,
                  alpha, beta);


        op.setJump(-2,0);
        //op.create(phi_exact, phi);
        //Vector<Real> exact(2,1.0);
        //op.setVal(phi_exact, exact);
        op.create(residual, rhs);
        op.setToZero(residual);



        op.residual(residual, phi, rhs);
        zeroAdjacentToBoundary(residual, pd[i]);
        Real r2norm = op.norm(residual, 2);
        Real r0norm = op.norm(residual, 0);


        pout()<<indent<<"Residual L2 norm "<<r2norm<<"  Residual max norm = "
              <<r0norm<<std::endl;
      }



      aliasMF(phaseA, 0, vresidual);
      aliasMF(phaseB, 1, vresidual);

#ifdef CH_USE_HDF5
      sprintf(iter_str, "residual.%03d.%dd.hdf5",step, SpaceDim);
      Vector<std::string> names(2); names[0]="residual0"; names[1]="residual1";
      writeEBHDF5(iter_str, grids, &phaseA, &phaseB, &levelset, names, domainCoar,
                  dxCoar, 1, step, refRatio, 3, true, coveredVal);
#endif
      aliasMF(phaseA, 0, mfvector);
      aliasMF(phaseB, 1, mfvector);
#ifdef CH_USE_HDF5
      sprintf(iter_str, "phi.%03d.%dd.hdf5",step, SpaceDim);
      names[0]="phi0"; names[1]="phi1";
      writeEBHDF5(iter_str, grids, &phaseA, &phaseB, &levelset, names, domainCoar,
                  dxCoar, 1, step, refRatio, 3, true, coveredVal);
#endif

      step++;

      center[0]-= dx/3.0;
      center[1]-= dx/2.0;
      radius += dx/6.0;

      pout()<<step<<std::endl;
    }

    pout() <<"\n "<<indent2<<pgmname<<"  test passed" << endl;

    for (int i=0; i<3; i++)
      {
        delete phaseA[i];
        delete phaseB[i];
      }

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
  for (int i=0; i<left.size(); ++i)
    {
      swap = left[i];
      left[i]=right[i];
      right[i]=swap;
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
                   const IntVectSet&    baseTags,
                   int a_maxSize)
{

  dbl.resize(3);
  Vector<int> refRatios(2,2);
  Real fillRatio = 0.85;
  int blockingFactor= 4;
  int bufferSize =   2;
  BRMeshRefine regridder(baseDomain, refRatios, fillRatio, blockingFactor,
                         bufferSize, a_maxSize);

  Vector<Vector<Box> > oldGrids(3,1), newGrids(3);
  oldGrids[0][0]=baseDomain.domainBox();
  oldGrids[1][0]=coarsen(oldGrids[0][0], 2);
  oldGrids[2][0]=coarsen(oldGrids[1][0], 2);

  regridder.regrid(newGrids, baseTags, 0, 1, oldGrids);

  Vector<int> procs;
  for (int i=0; i<3; i++)
  {
    LoadBalance(procs, newGrids[i]);
    dbl[i] = DisjointBoxLayout(newGrids[i], procs);
  }

}
