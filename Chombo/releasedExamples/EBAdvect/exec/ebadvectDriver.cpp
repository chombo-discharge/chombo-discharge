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

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"

#include "EBFABView.H"
#include "EBDebugDump.H"

#include "EBLevelDataOps.H"
#include "BaseIF.H"
#include "SphereIF.H"
#include "PlaneIF.H"
#include "GeometryShop.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLevelDataOps.H"
#include "LoadBalance.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"
#include "EBPhysIBC.H"
#include "EBArith.H"
#include "ExtrapAdvectBC.H"
#include "EBLevelRedist.H"
#include "EBAdvectPatchIntegrator.H"
#include "EBNormalizeByVolumeFraction.H"
#include "IrregFABFactory.H"
#include "SetupFunctions.H"
#include "memusage.H"
#include "memtrack.H"

#define PI 3.141592653589793

ExtrapAdvectBCFactory s_physibc;
/********/
void
getExtrapolatedScalar(LevelData<EBFluxFAB> & a_extrap,
                      LevelData<EBCellFAB> & a_scala,
                      LevelData<EBCellFAB> & a_ccVel,
                      LevelData<EBFluxFAB> & a_fcVel,
                      const EBLevelGrid    & a_eblg,  
                      Real a_dx, Real a_dt)
{
  Real time = 0; //
  DisjointBoxLayout grids = a_eblg.getDBL();
  DataIterator dit = grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBCellFAB& scala = a_scala[dit[ibox]];
    EBCellFAB& ccVel = a_ccVel[dit[ibox]];
    EBFluxFAB& fcVel = a_fcVel[dit[ibox]];

    Box     grid = a_eblg.getDBL()  [dit[ibox]];
    EBISBox ebis = a_eblg.getEBISL()[dit[ibox]];
    ProblemDomain domain = a_eblg.getDomain();
    RealVect dxvect = a_dx*RealVect::Unit;
    IntVectSet cfivs;
    EBAdvectPatchIntegrator ebpatchad(grid, ebis, cfivs, domain, dxvect, true);


    ebpatchad.setVelocities(ccVel, fcVel);
    ebpatchad.setEBPhysIBC(s_physibc);
    ebpatchad.setDoingVel(0);
    ebpatchad.setCurComp(0);

    EBCellFAB source(ebis, ccVel.box(), 1);
    source.setVal(0.);
    EBFluxFAB& extrapState = a_extrap[dit[ibox]];
    ebpatchad.extrapolateBCG(extrapState,
                             scala,
                             source,
                             dit[ibox], time, a_dt);
    ideb++;

  }
}
/********/
void
getFaceCenteredFlux(LevelData<EBFluxFAB> & a_fcFlx,
                    LevelData<EBFluxFAB> & a_fcSca,
                    LevelData<EBFluxFAB> & a_fcVel,
                    const EBLevelGrid    & a_eblg)
{
  DisjointBoxLayout grids = a_eblg.getDBL();

  DataIterator dit = grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBFluxFAB& fcVel  = a_fcVel[dit[ibox]];
    EBFluxFAB& fcSca  = a_fcSca[dit[ibox]];
    EBFluxFAB& fcFlx  = a_fcFlx[dit[ibox]];

    fcFlx.setVal(0.);
    fcFlx += fcSca;
    fcFlx *= fcVel;
    ideb++;
  }
}
/********/
void
getCentroidFlux(LevelData<EBFluxFAB> & a_centroidFlux,
                LevelData<EBFluxFAB> & a_faceCentFlux,
                const EBLevelGrid    & a_eblg)
{
  DisjointBoxLayout grids = a_eblg.getDBL();
  EBISLayout        ebisl = a_eblg.getEBISL();
  ProblemDomain    domain = a_eblg.getDomain();

  DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBISBox  ebis = ebisl[dit[ibox]];
    Box      grid = grids[dit[ibox]];

    EBFluxFAB& centroidFlux  = a_centroidFlux[dit[ibox]];
    EBFluxFAB& faceCentFlux  = a_faceCentFlux[dit[ibox]];
    IntVectSet ivs = ebis.getIrregIVS(grid);
    IntVectSet cfivs;
    FaceStop::WhichFaces facestop = FaceStop::SurroundingWithBoundary;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      EBFaceFAB& roid = centroidFlux[idir];
      EBFaceFAB& cent = faceCentFlux[idir];
      roid.copy(cent);
      for (FaceIterator faceit(ivs, ebis.getEBGraph(), idir, facestop); faceit.ok(); ++faceit)
      {
        FaceIndex face = faceit();
        FaceStencil sten = EBArith::getInterpStencil(face, cfivs, ebis, domain);
        Real roidval = 0;
        IntVect ivhi = face.getVoF(Side::Hi).gridIndex();
        int ideb = 0;
        if((ivhi[0] == 16) && (ivhi[1] == 18))
        {
          ideb = 1;
        }
        for(int isten = 0; isten < sten.size(); isten++)
        {
          const FaceIndex& stenFace= sten.face(isten);
          Real weight  = sten.weight(isten);
          Real centval = cent(stenFace, 0);
          roidval += weight*centval;
        }
        roid(face, 0) = roidval;
      }// Attention all planets of the Solar Federation.
    }//   Attention all planets of the Solar Federation.
  }//     We have assumed control.   
} //      We have assumed control.   


/********/
void
getKappaDivFlux(LevelData<EBCellFAB> & a_kappaDivFlux,
                LevelData<EBFluxFAB> & a_centroidFlux,
                LevelData<EBCellFAB> & a_ccVel,
                LevelData<EBFluxFAB> & a_fcVel,
                const EBLevelGrid    & a_eblg,
                Real a_dx, Real a_dt)
{
  DisjointBoxLayout grids = a_eblg.getDBL();
  EBISLayout        ebisl = a_eblg.getEBISL();
  ProblemDomain    domain = a_eblg.getDomain();

  DataIterator dit = grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box     grid = a_eblg.getDBL()  [dit[ibox]];
    EBISBox ebis = a_eblg.getEBISL()[dit[ibox]];

    IntVectSet cfivs;
    RealVect dxvect = a_dx*RealVect::Unit;
    EBAdvectPatchIntegrator ebpatchad(grid, ebis, cfivs, domain, dxvect, true);


    EBCellFAB& ccVel = a_ccVel[dit[ibox]];
    EBFluxFAB& fcVel = a_fcVel[dit[ibox]];
    ebpatchad.setVelocities(ccVel, fcVel);
    ebpatchad.setEBPhysIBC(s_physibc);
    ebpatchad.setDoingVel(0);
    ebpatchad.setCurComp(0);


    EBCellFAB& kappaDivF = a_kappaDivFlux[dit[ibox]];
    EBFluxFAB& centroidF = a_centroidFlux[dit[ibox]];
    ebpatchad.kappaDivergenceFlux(kappaDivF,
                                  centroidF,
                                  grid);
    ideb++;
  }
} 
/********/
void
hybridizeDivergence(LevelData<EBCellFAB> & a_hybridDivF,
                    LevelData<EBCellFAB> & a_kappaDivF,
                    LevelData<EBCellFAB> & a_nonConsDivF,
                    LevelData<IrregFAB>  & a_deltaM,
                    const EBLevelGrid    & a_eblg)
{
  DisjointBoxLayout grids = a_eblg.getDBL();
  EBISLayout        ebisl = a_eblg.getEBISL();
  ProblemDomain    domain = a_eblg.getDomain();

  int ideb = 0;
  DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBISBox  ebis = ebisl[dit[ibox]];
    Box      grid = grids[dit[ibox]];

    EBCellFAB& hybridDivF   =  a_hybridDivF[dit[ibox]];
    EBCellFAB& kappaDivF    =   a_kappaDivF[dit[ibox]];
    EBCellFAB& nonConsDivF  = a_nonConsDivF[dit[ibox]];
    IrregFAB & deltaM       =      a_deltaM[dit[ibox]];

    hybridDivF.copy(kappaDivF);
    IntVectSet ivs = ebis.getIrregIVS(grid);
    for (VoFIterator vofit(ivs, ebis.getEBGraph()); vofit.ok(); ++vofit)
    {
      VolIndex vof = vofit();
      Real kappa  = ebis.volFrac(vof);
      Real kdiv   = kappaDivF(vof, 0);
      Real ncdiv  = nonConsDivF(vof, 0);
      Real divval = kdiv + (1. - kappa)*ncdiv;
      Real massval= (1.-kappa)*(kdiv - kappa*ncdiv);

      hybridDivF(vof, 0)  = divval;
      deltaM(    vof, 0)  = massval;
    }
    ideb++;
  }
}

/********/
void
advance(LevelData<EBCellFAB> & a_scala,
        LevelData<EBCellFAB> & a_ccVel,
        LevelData<EBFluxFAB> & a_fcVel,
        const EBLevelGrid    & a_eblg,  
        Real a_dx, Real a_dt)
{
  DisjointBoxLayout grids = a_eblg.getDBL();
  EBISLayout        ebisl = a_eblg.getEBISL();
  EBCellFactory    cellfact(ebisl);
  EBFluxFactory    fluxfact(ebisl);
  IrregFABFactory irregfact(ebisl);
  LevelData<EBCellFAB>    kappaDivF(grids, 1, IntVect::Unit, cellfact);
  LevelData<EBCellFAB>  nonConsDivF(grids, 1, IntVect::Unit, cellfact);
  LevelData<EBCellFAB>   hybridDivF(grids, 1, IntVect::Unit, cellfact);
  LevelData<EBFluxFAB>    macScalar(grids, 1, IntVect::Unit, fluxfact);
  LevelData<EBFluxFAB> faceCentFlux(grids, 1, IntVect::Unit, fluxfact);
  LevelData<EBFluxFAB> centroidFlux(grids, 1, IntVect::Unit, fluxfact);
  LevelData<IrregFAB>        deltaM(grids, 1, IntVect::Unit,irregfact);
  getExtrapolatedScalar(macScalar, a_scala, a_ccVel, a_fcVel, a_eblg, a_dx, a_dt);
  macScalar.exchange();
  getFaceCenteredFlux(faceCentFlux, macScalar, a_fcVel, a_eblg);
  faceCentFlux.exchange();
  getCentroidFlux(centroidFlux, faceCentFlux, a_eblg);
  getKappaDivFlux(kappaDivF,    centroidFlux, a_ccVel, a_fcVel, a_eblg, a_dx, a_dt);

  Interval interv(0,0);
  kappaDivF.copyTo(interv, nonConsDivF, interv);
  EBNormalizeByVolumeFraction normalizor(a_eblg);
  normalizor(nonConsDivF);

  hybridizeDivergence(hybridDivF, kappaDivF, nonConsDivF, deltaM, a_eblg);

  EBLevelRedist smush(a_eblg, 1);
  smush.setToZero();
  smush.increment(deltaM, interv);
  smush.redistribute(hybridDivF, interv);
  int ideb = 0;
  DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBISBox ebis = ebisl[dit[ibox]];
    Box     grid = grids[dit[ibox]];
    EBCellFAB incr(ebis, grid, 1);
    EBCellFAB& scal = a_scala[dit[ibox]];
    EBCellFAB& divf = hybridDivF[dit[ibox]];
    incr.copy(divf);
    incr *= -a_dt;
    scal += incr;

    ideb++;
  }
}
/********/
void  computeDt(Real                 & a_dt, 
                LevelData<EBCellFAB> & a_ccVel, 
                Real a_dx, Real a_cfl)
{
  //debug
  dumpEB1(NULL);
  dumpBlob(NULL);
  dumpFace1(NULL);
  //end debug
  Real maxval = 0;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    Real maxloc, minloc;
    EBLevelDataOps::getMaxMin(maxloc, minloc, a_ccVel, idir, true);
    maxval = std::max(maxval, maxloc);
  }
  if(maxval > 0)
    a_dt = a_cfl*a_dx/maxval;

  pout() << "computeDt: dt = " << a_dt << endl;
}
                
/********/
void runAdvection()
{
  CH_TIME("run advection");

  Real coveredval = -1;
  Real cfl    = 0.5;
  int nx      = 32;
  int  max_step   = 10;
  Real max_time   = 1.0;
  Real max_vel_mag = 1.0;
  Real max_vel_rad = 0.25;
  int nStream    = 8;
  int outputInterval = -1;
  ParmParse pp;

  pp.get("nstream", nStream);

    
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);
  pp.get("cfl"  , cfl);
  pp.get("max_vel_mag"  , max_vel_mag);
  pp.get("max_vel_rad"  , max_vel_rad);

  pout() << "num_streams     = " << nStream         << endl;
  pout() << "max_step        = " << max_step        << endl;
  pout() << "max_time        = " << max_time        << endl;
  pout() << "output interval = " << outputInterval  << endl;
  pout() << "cfl             = " << cfl             << endl;
  pout() << "max_vel_mag     = " << max_vel_mag     << endl;
  pout() << "max_vel_rad     = " << max_vel_rad     << endl;

  Real dx;
  DisjointBoxLayout grids;

  pout() << "defining geometry" << endl;


  Real geomCen;
  Real geomRad;
  Real blobCen;
  Real blobRad;
  int whichGeom;
  defineGeometry(grids, dx, geomCen, geomRad, blobCen, blobRad, whichGeom, nx);
  int nghost = 4;
  IntVect dataGhostIV =   nghost*IntVect::Unit;


  Box domain = grids.physDomain().domainBox();
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  EBLevelGrid eblg(grids, domain, nghost, ebisPtr);


  pout() << "inititializing data"   << endl;
  EBCellFactory cellfact(eblg.getEBISL());
  EBFluxFactory fluxfact(eblg.getEBISL());
  LevelData<EBCellFAB>  scala(grids,        1, dataGhostIV,  cellfact);
  LevelData<EBCellFAB>  ccVel(grids, SpaceDim, dataGhostIV,  cellfact);
  LevelData<EBFluxFAB>  fcVel(grids,        1, dataGhostIV,  fluxfact);
  
  InitializeStuff(scala, ccVel, fcVel, grids, dx, geomCen, geomRad, blobCen, blobRad, max_vel_mag, max_vel_rad);


  int step = 0; Real time = 0;
  Real dt = 0;
  pout() << "computing the time step"  << endl;
  computeDt(dt, ccVel, dx, cfl);

  if(outputInterval > 0)
  {
    string filev("velo.hdf5");
    string filep("scal.0.hdf5");
    writeEBLevelName(ccVel, filev);
    writeEBLevelName(scala, filep);
  }

  pout() << "running advection operator " << endl;

  while((step < max_step) && (time < max_time))
  {
    advance(scala, ccVel, fcVel, eblg, dx, dt);
    scala.exchange();

    pout() <<" step = " << step << " time = " << time << " time step = " << dt << endl;
    step++;
    time += dt;

    if((outputInterval > 0) && ( (step%outputInterval == 0) || step == (max_step-1)))
    {
      string filep = string("scal.") + std::to_string(step) + string(".hdf5");
      writeEBLevelName(scala,  filep);
    }
  }
  pout() << "exiting runAdvection" << endl;


}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    //solve the stinking problem and output everything
    runAdvection();

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  }
  // End scoping trick

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
