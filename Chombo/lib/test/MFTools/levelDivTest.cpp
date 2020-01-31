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

#include "FABView.H"
#include "SphereIF.H"
#include "MFIndexSpace.H"
#include "MFCellFAB.H"
#include "MFRemapper.H"
#include "MFAliasFactory.H"

#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "EBAMRIO.H"
#include "BoxIterator.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "EBArith.H"
#include "EBPatchGodunovF_F.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "MFFABView.H"
#include "UsingNamespace.H"

Real pi = 2*atan2((Real)1,(Real)0);

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



void
compareError(const LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const LevelData<EBCellFAB>& a_errorMed,
             const EBISLayout& a_ebislMed,
             const DisjointBoxLayout& a_gridMed,
             const Box& a_domainMed,
             const LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar);

void
getError(LevelData<EBCellFAB>&       a_error,
         const EBISLayout&           a_ebisl,
         const DisjointBoxLayout&    a_dbl,
         const Real&                 a_dx);

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

  // begin forever present scoping trick
  {
    const char* in_file = "sphere.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    RealVect center;
    Real radius;
    RealVect origin;
    Real dx;
    Box domain;
    int eekflag = 0;
    MFIndexSpace a, b;
    MFIndexSpace *mfIndexSpace, *oldMFIndexSpace;
    mfIndexSpace = &a;
    oldMFIndexSpace = &b;

    Vector<std::string> names0(1, "phi0");
    Vector<std::string> names1(1, "phi1");
    Vector<std::string> names(2);
    names[0] = "phi0";
    names[1] = "phi1";

    Vector<int> refRatio(3,2);
    Vector<Real> coveredVal(1,7.0);

    LevelData<MFCellFAB> fine, med, coarse;
    Vector<LevelData<MFCellFAB>* > mfvector(3,NULL);
    mfvector[0]=&coarse;
    mfvector[1]=&med;
    mfvector[2]=&fine;

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

    int maxsize = 16;
    Vector<Box> boxCoarse;
    domainSplit(domainCoar, boxCoarse, maxsize);

    Vector<int> procs(boxCoarse.size(), 0);
    LoadBalance(procs, boxCoarse);
    grids[0] = DisjointBoxLayout(boxCoarse, procs);
    refine(grids[1], grids[0], 2);
    refine(grids[2], grids[1], 2);

    //make data holders
    Vector<int> comps(2,1);
    int ghost = 1;

    int steps= 10;
    int step = 0;

    Vector<LevelData<EBCellFAB>* > phase(3);
    for (int i=0; i<3; i++) phase[i] = new LevelData<EBCellFAB>();

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
          MayDay::Error("problem in makeGeometry");
        }

      MFCellFactory fineFactory(*mfIndexSpace, grids[2],
                                domainFine,
                                comps, ghost);
      MFCellFactory medFactory(*mfIndexSpace, grids[1],
                               domainMedi,
                               comps, ghost);
      MFCellFactory coarseFactory(*mfIndexSpace, grids[0],
                                  domainCoar,
                                  comps, ghost);

      mfvector[2]->define(grids[2], 1, IntVect::Unit,
                          fineFactory);

      mfvector[1]->define(grids[1], 1, IntVect::Unit,
                          medFactory);

      mfvector[0]->define(grids[0], 1, IntVect::Unit,
                          coarseFactory);

      // initialize data holders

      aliasMF(phase, 0, mfvector);

      getError(*phase[2], fineFactory.getEBISLayout(0),    grids[2], dxFine);
      getError(*phase[1] , medFactory.getEBISLayout(0),     grids[1], dxMedi);
      getError(*phase[0] , coarseFactory.getEBISLayout(0),  grids[0], dxCoar);
//       sprintf(iter_str, "phase0.%03d.%dd.hdf5",step, SpaceDim);
//       writeEBHDF5(iter_str, grids, phase, names0, domainCoar,
//                   dxCoar, 1, step, refRatio, 3, true, coveredVal);


      compareError(*phase[2], fineFactory.getEBISLayout(0), grids[2], domainFine,
                   *phase[1] , medFactory.getEBISLayout(0),  grids[1], domainMedi,
                *phase[0] , coarseFactory.getEBISLayout(0),  grids[0], domainCoar);

      aliasMF(phase, 1, mfvector);

      getError(*phase[2], fineFactory.getEBISLayout(1),    grids[2], dxFine);
      getError(*phase[1] , medFactory.getEBISLayout(1),     grids[1], dxMedi);
      getError(*phase[0] , coarseFactory.getEBISLayout(1),  grids[0], dxCoar);
//       sprintf(iter_str, "phase1.%03d.%dd.hdf5",step, SpaceDim);
//       writeEBHDF5(iter_str, grids, phase, names1, domainCoar,
//                   dxCoar, 1, step, refRatio, 3, true, coveredVal);


      compareError(*phase[2], fineFactory.getEBISLayout(1), grids[2], domainFine,
                   *phase[1] , medFactory.getEBISLayout(1),  grids[1], domainMedi,
                *phase[0] , coarseFactory.getEBISLayout(1),  grids[0], domainCoar);


      step++;

      center[0]+= dx/2.0;
      center[1]+= dx/2.0;
      radius -= dx/4.0;
      {
        MFIndexSpace* swap = oldMFIndexSpace;
        oldMFIndexSpace = mfIndexSpace;
        mfIndexSpace = swap;
      }
      pout()<<step<<std::endl;
    }
    for (int i=0; i<3; i++) delete phase[i];

    pout() << "\n levelDivTest test passed" << endl;
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

void makeHierarchy(Vector<DisjointBoxLayout>& dbl,
                   const ProblemDomain& baseDomain,
                   const IntVectSet&    baseTags)
{

  dbl.resize(3);
  Vector<int> refRatios(2,2);
  Real fillRatio = 0.85;
  int blockingFactor= 4;
  int bufferSize =   2;
  int maxSize    =   32;
  BRMeshRefine regridder(baseDomain, refRatios, fillRatio, blockingFactor,
                         bufferSize, maxSize);

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
/************/
Real
exactDivergence(const RealVect& a_xval)
{
  Real retval;
  Real x = a_xval[0];
  Real y = a_xval[1];
  //retval = 0.0;
  retval = 3*x*x + pi*cos(y*pi);

  return retval;
}
/************/
Real
exactFlux(const RealVect& a_xval, int facedir)
{
  RealVect retval = RealVect::Zero;
  Real x = a_xval[0];
  Real y = a_xval[1];
  //retval[0] = 2;
  //retval[1] = 5;
   retval[0] = x*x*x;
   retval[1] = sin(y*pi);

  return retval[facedir];
}

/************/
void
setToExactDivF(EBCellFAB&     a_exactDivF,
               const EBISBox& a_ebisBox,
               const Box&     a_region,
               const Real&    a_dx)
{
  a_exactDivF.setVal(0.);
  IntVectSet ivsregion(a_region);
  for (VoFIterator vofit(ivsregion, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect xval;
      IntVect iv = vof.gridIndex();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xval[idir] = (Real(iv[idir]) + 0.5)*a_dx;
        }
      Real solnrv = exactDivergence(xval);
      Real kappa = a_ebisBox.volFrac(vof);
      a_exactDivF(vof,0) = kappa*solnrv;
    }
}
/************/
void
setToExactFlux(EBFluxFAB&             a_flux,
               const EBISBox&         a_ebisBox,
               const Box&             a_region,
               const Real&            a_dx)
{
  IntVectSet ivsregion(a_region);
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (FaceIterator faceit(ivsregion, a_ebisBox.getEBGraph(), faceDir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
        {
          RealVect xval;
          IntVect iv = faceit().gridIndex(Side::Hi);
          RealVect centroid = a_ebisBox.centroid(faceit());
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir == faceDir)
                {
                  xval[idir] = (Real(iv[idir]))*a_dx;
                }
              else
                {
                  xval[idir] = (Real(iv[idir]) + 0.5 + centroid[idir])*a_dx;
                }
            }
          Real fluxDir = exactFlux(xval, faceDir);
          a_flux[faceDir](faceit(), 0) = fluxDir;
        }
    } //end loop over face directions

}
/************/
void
kappaDivergence(EBCellFAB&             a_divF,
                const EBFluxFAB&       a_flux,
                const EBISBox&         a_ebisBox,
                const Box&             a_box,
                const Real&            a_dx)
{
  //set the divergence initially to zero
  //then loop through directions and increment the divergence
  //with each directions flux difference.
  a_divF.setVal(0.0);
  BaseFab<Real>&       regDivF = a_divF.getSingleValuedFAB();
  regDivF.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //update for the regular vofs in the nonconservative
      //case  works for all single valued vofs.
      /* do the regular vofs */
      /**/
      const EBFaceFAB& fluxDir = a_flux[idir];
      const BaseFab<Real>& regFluxDir = fluxDir.getSingleValuedFAB();
      int ncons = 1;
      FORT_DIVERGEF( CHF_BOX(a_box),
                     CHF_FRA(regDivF),
                     CHF_CONST_FRA(regFluxDir),
                     CHF_CONST_INT(idir),
                     CHF_CONST_INT(ncons),
                     CHF_CONST_REAL(a_dx));
      /**/
    }
  //update the irregular vofs using conservative diff
  IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(a_box);
  for (VoFIterator vofit(ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //divergence was set in regular update.  we reset it
      // to zero and recalc.
      Real update = 0.;
      for ( int idir = 0; idir < SpaceDim; idir++)
        {
          const EBFaceFAB& fluxDir = a_flux[idir];
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int isign = sign(sit());
              Vector<FaceIndex> faces =
                a_ebisBox.getFaces(vof, idir, sit());
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& face = faces[iface];
                  Real areaFrac = a_ebisBox.areaFrac(face);
                  Real faceFlux =fluxDir(face, 0);
                  update += isign*areaFrac*faceFlux;

                }
            }
        }
      //add EB boundary condtions in divergence
      const IntVect& iv = vof.gridIndex();
      Real bndryArea = a_ebisBox.bndryArea(vof);
      RealVect bndryCent = a_ebisBox.bndryCentroid(vof);
      RealVect normal = a_ebisBox.normal(vof);
      RealVect bndryLoc;
      RealVect exactF;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          bndryLoc[idir] = a_dx*(iv[idir] + 0.5 + bndryCent[idir]);
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          exactF[idir] = exactFlux(bndryLoc, idir);
        }
      Real bndryFlux = PolyGeom::dot(exactF, normal);

      update -= bndryFlux*bndryArea;
      update /= a_dx; //note NOT divided by volfrac

      a_divF(vof, 0) = update;
    }
}
/***************/
void
kappaDivergenceLD(LevelData<EBCellFAB>&       a_divF,
                  const LevelData<EBFluxFAB>& a_flux,
                  const EBISLayout&           a_ebisl,
                  const DisjointBoxLayout&    a_dbl,
                  const Real&                 a_dx)
{
  for (DataIterator dit= a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      kappaDivergence(a_divF[dit()],
                      a_flux[dit()],
                      a_ebisl[dit()],
                      a_dbl.get(dit()),
                      a_dx);
    }

}

/**************************/
void
setToExactFluxLD(LevelData<EBFluxFAB>&       a_flux,
                 const EBISLayout&           a_ebisl,
                 const DisjointBoxLayout&    a_dbl,
                 const Real&                 a_dx)
{
  for (DataIterator dit= a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      setToExactFlux(a_flux[dit()],
                     a_ebisl[dit()],
                     a_dbl.get(dit()),
                     a_dx);
    }
}
/**************************/
void
setToExactDivFLD(LevelData<EBCellFAB>&       a_soln,
                 const EBISLayout&           a_ebisl,
                 const DisjointBoxLayout&    a_dbl,
                 const Real&                 a_dx)
{
  for (DataIterator dit= a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      setToExactDivF(a_soln[dit()],
                     a_ebisl[dit()],
                     a_dbl.get(dit()),
                     a_dx);
    }
}
/**************************/
void
compareError(const LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const LevelData<EBCellFAB>& a_errorMed,
             const EBISLayout& a_ebislMed,
             const DisjointBoxLayout& a_gridsMed,
             const Box& a_domainMed,
             const LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar)
{
  CH_assert(a_errorFine.nComp() == 1);
  CH_assert(a_errorCoar.nComp() == 1);

  pout() << "==============================================" << endl;

  EBNormType::NormMode normtype = EBNormType::OverBoth;
  pout() << endl << "Using all uncovered cells." << endl  ;

  for (int inorm = 0; inorm <= 2; inorm++)
    {

      if (inorm == 0)
        {
          pout() << endl << "Using max norm." << endl;
        }
      else
        {
          pout() << endl << "Using L-" << inorm << "norm." << endl;
        }
      int comp = 0;
      Real coarnorm = EBArith::norm(a_errorCoar,
                                    a_gridsCoar, a_ebislCoar,
                                    comp, inorm, normtype);
      Real mednorm  = EBArith::norm(a_errorMed,
                                    a_gridsMed, a_ebislMed,
                                    comp, inorm, normtype);

      Real finenorm = EBArith::norm(a_errorFine,
                                    a_gridsFine, a_ebislFine,
                                    comp, inorm, normtype);

      pout() << "Coarse Error Norm = " << coarnorm << endl;
      pout() << "Medium Error Norm = " << mednorm << endl;
      pout() << "Fine   Error Norm = " << finenorm << endl;

      if ((Abs(finenorm) > 1.0e-12) && (Abs(coarnorm) > 1.0e-12))
        {
          Real order = log(Abs(coarnorm/mednorm))/log(2.0);
          pout() << "Order of scheme (coar to med)= " << order << endl;

          order = log(Abs(mednorm/finenorm))/log(2.0);
          pout() << "Order of scheme (med to fine)= " << order << endl;
        }
    }
  pout() << "==============================================" << endl ;;
}

/**************************/
void
getError(LevelData<EBCellFAB>&       a_error,
         const EBISLayout&           a_ebisl,
         const DisjointBoxLayout&    a_dbl,
         const Real&                 a_dx)
{
  EBCellFactory ebcellfact(a_ebisl);
  EBFluxFactory ebfluxfact(a_ebisl);
  a_error.define(a_dbl, 1, IntVect::Zero,   ebcellfact);
  LevelData<EBCellFAB> divFCalc(a_dbl, 1, IntVect::Zero,   ebcellfact);
  LevelData<EBCellFAB> divFExac(a_dbl, 1, IntVect::Zero,   ebcellfact);
  LevelData<EBFluxFAB> faceFlux(a_dbl, 1, IntVect::Zero,   ebfluxfact);

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      a_error[dit()].setVal(0.);
      divFCalc[dit()].setVal(0.);
      divFExac[dit()].setVal(0.);
    }

  setToExactDivFLD(divFExac,  a_ebisl, a_dbl, a_dx);
  setToExactFluxLD(faceFlux,  a_ebisl, a_dbl, a_dx);

  Interval interv(0, 0);
  faceFlux.exchange(interv);

  kappaDivergenceLD(divFCalc, faceFlux, a_ebisl, a_dbl, a_dx);

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& errorFAB = a_error[dit()];
      EBCellFAB& exactFAB = divFExac[dit()];
      EBCellFAB& calcuFAB = divFCalc[dit()];

      errorFAB += calcuFAB;
      errorFAB -= exactFAB;
    }
}
