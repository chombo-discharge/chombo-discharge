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
#include <cstdio>
#include <iostream>


#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "CONSTANTS.H"

#include "PoissonUtilities.H"
#include "EBAMRPoissonOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "EBGeomDebugDump.H"
#include "DebugDump.H"
#include "EBAMRDataOps.H"
#include "MonomialIF.H"
#include "WrappedGShop.H"
#include "MomentIterator.H"
#include "Factorial.H"
#include "CH_Attach.H"
/******/
#include "CutCellMoments.H"

/******/
enum TestType
  {
    VOL_MOM = 0,
    EB_MOM,
    EB_NORMAL_MOM,
    NUM_TESTS
  };


/************/
void
putIDMIntoFAB(EBCellFAB                            &  a_datum,
             IndMomSpaceDim &  a_testOrder, 
              const VolIndex                       &  a_vof)
{
  MomentIterator<SpaceDim, CH_EBIS_ORDER> volmomit;
  for(volmomit.reset(); volmomit.ok(); ++volmomit)
    {
      int ivar             = a_testOrder.indexOf(volmomit());
      a_datum(a_vof, ivar) = a_testOrder[volmomit()];
    }
}
/************/
void
generateData(LevelData<EBCellFAB>       & a_datum,
             const DisjointBoxLayout    & a_grids,
             const EBISLayout           & a_ebisl,
             const GridParameters    & a_params,
             const TestType             & a_type,
             const int                  & a_idir,
             const bool                 & a_shiftToCoar)
{
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivs;
      if(a_type == VOL_MOM)
        {
          ivs = IntVectSet(a_grids[dit()]);
        }
      else
        {
          ivs = a_ebisl[dit()].getIrregIVS(a_grids[dit()]);
        }
      for(VoFIterator vofit(ivs, a_ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          IndMomSpaceDim ebisOrder;
          if(a_type == VOL_MOM)
            {
              ebisOrder = a_ebisl[dit()].getEBData().getVolumeMoments(vofit());
            }
          else if(a_type == EB_MOM)
            {
              ebisOrder = a_ebisl[dit()].getEBData().getEBMoments(vofit());
            }
          else if(a_type == EB_NORMAL_MOM)
            {
              ebisOrder = a_ebisl[dit()].getEBData().getEBNormalMoments(vofit(), a_idir);
            }
          if(a_shiftToCoar)
            {
              Real fineDx = a_params.coarsestDx;
              Real coarDx = 2.*fineDx;
              VolIndex fineVoF= vofit();
              IntVect coarIV = coarsen(fineVoF.gridIndex(), 2);
              VolIndex coarVoF(coarIV, 0);
              RealVect fineLoc = EBArith::getVoFLocation(fineVoF, fineDx*RealVect::Unit, RealVect::Zero);
              RealVect coarLoc = EBArith::getVoFLocation(coarVoF, coarDx*RealVect::Unit, RealVect::Zero);
              RealVect shiftRV = fineLoc - coarLoc;
              RvSpaceDim shiftitm;
              EBArith::convertToITM(shiftitm, shiftRV);
              ebisOrder.shift(shiftitm);
            }
          putIDMIntoFAB(a_datum[dit()], ebisOrder, vofit());
        }
    }
}
/*******/
void
sumFineValues(Vector<Real>           & a_fineSum, 
              const EBCellFAB        & a_solutFine,
              const Vector<VolIndex> & a_fineVoFs)
{
  CH_assert(a_fineSum.size() == a_solutFine.nComp());
  for(int ivar = 0; ivar < a_solutFine.nComp(); ivar++)
    {
      Real value = 0;
      for(int ivof = 0; ivof < a_fineVoFs.size(); ivof++)
        {
          value += a_solutFine(a_fineVoFs[ivof], ivar);
        }
      a_fineSum[ivar] = value;
    }
}
/*******/
void
sumFineMinusCoarse(LevelData<EBCellFAB>       & a_errorMedi,
                   const LevelData<EBCellFAB> & a_solutMedi, 
                   const DisjointBoxLayout    & a_gridsMedi,
                   const EBISLayout           & a_ebislMedi, 
                   const GridParameters    & a_paramMedi,
                   const LevelData<EBCellFAB> & a_solutFine, 
                   const DisjointBoxLayout    & a_gridsFine,
                   const EBISLayout           & a_ebislFine, 
                   const GridParameters    & a_paramFine,
                   const TestType& a_type)
{
  //fine has to be shifted to coarse location
  int ncomp = a_solutMedi.nComp();
  for(DataIterator dit = a_gridsMedi.dataIterator(); dit.ok(); ++dit)
    {
      a_errorMedi[dit()].setVal(0.);
      IntVectSet ivs;
      if(a_type == VOL_MOM)
        {
          ivs = IntVectSet(a_gridsMedi[dit()]);
        }
      else
        {
          ivs = a_ebislMedi[dit()].getIrregIVS(a_gridsMedi[dit()]);
        }
      for(VoFIterator vofit(ivs, a_ebislMedi[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vofMedi = vofit();
          //because this is a screwy geometry test, the graph coarsening can be weird here
          // small cells on the coarse level can be covered where they were not on the finer level
          EBISBox ebisBoxMedi = a_errorMedi[dit()].getEBISBox();
          if(!ebisBoxMedi.isCovered(vofMedi.gridIndex()))
            {
              Vector<VolIndex> fineVoFs = a_ebislMedi[dit()].refine(vofMedi);
              Vector<Real> fineSum(ncomp, 0.);
              //these are moments so we just sum the values (sum of integrals = integral of sum)
              sumFineValues(fineSum, a_solutFine[dit()], fineVoFs);
              for(int ivar = 0; ivar < ncomp; ivar++)
                {
                  a_errorMedi[dit()](vofMedi, ivar) = a_solutMedi[dit()](vofMedi, ivar) - fineSum[ivar];
                }
            }
        }
    }
}

/*****/
void
solutionErrorTest(const GridParameters &  a_paramFine, 
                  const GridParameters &  a_paramMedi, 
                  const GridParameters &  a_paramCoar,
                  const DisjointBoxLayout &  a_gridsFine, 
                  const DisjointBoxLayout &  a_gridsMedi, 
                  const DisjointBoxLayout &  a_gridsCoar, 
                  const EBISLayout        &  a_ebislFine, 
                  const EBISLayout        &  a_ebislMedi, 
                  const EBISLayout        &  a_ebislCoar, 
                  const EBISLayout        &  a_ebislCoFi, 
                  const EBISLayout        &  a_ebislCoMe, 
                  const TestType          &  a_type, 
                  const int               &  a_idir)
{
  IndMomSpaceDim idmproxy;
  int nvar = idmproxy.size();

  string prefix;
  if(a_type == VOL_MOM)
    {
      prefix = string("volume_moment");
    }
  else if(a_type == EB_MOM)
    {
      prefix = string("eb_moment");
    }
  else if(a_type == EB_NORMAL_MOM)
    {
      prefix =  string("ebNormalMoment_") + convertInt(a_idir);
    }
  else
    {
      MayDay::Error("bogus type");
    }

  EBCellFactory       factFine(a_ebislFine);
  EBCellFactory       factMedi(a_ebislMedi);
  EBCellFactory       factCoar(a_ebislCoar);

  LevelData<EBCellFAB>  solutFine(a_gridsFine, nvar, IntVect::Zero, factFine);
  LevelData<EBCellFAB>  solutMedi(a_gridsMedi, nvar, IntVect::Zero, factMedi);
  LevelData<EBCellFAB>  solutCoar(a_gridsCoar, nvar, IntVect::Zero, factCoar);
  LevelData<EBCellFAB>  errorMedi(a_gridsMedi, nvar, IntVect::Zero, factMedi);
  LevelData<EBCellFAB>  errorCoar(a_gridsCoar, nvar, IntVect::Zero, factCoar);

  EBLevelDataOps::setToZero(solutFine);
  EBLevelDataOps::setToZero(solutMedi);
  EBLevelDataOps::setToZero(solutCoar);

  //fine has to be shifted to coarse location
  bool shiftToCoar;
  //need to shift to coarse locations so this will be at the same location as the coarse
  shiftToCoar = true;
  pout() << "generating fine solution" << endl;
  generateData(solutFine, a_gridsFine, a_ebislFine, a_paramFine, a_type, a_idir, shiftToCoar);

  //for this bit, medi is the coarse solution so no shifting
  shiftToCoar = false;
  pout() << "generating medi solution" << endl;
  generateData(solutMedi, a_gridsMedi, a_ebislMedi, a_paramMedi, a_type, a_idir, shiftToCoar);

  pout() << "generating medi error from medi and fine solutions" << endl;
  sumFineMinusCoarse(errorMedi,
                     solutMedi, a_gridsMedi, a_ebislCoFi, a_paramMedi,
                     solutFine, a_gridsFine, a_ebislFine, a_paramFine, a_type);

  //for this bit, medi is the finer solution so it has to get shifted
  shiftToCoar = true;
  pout() << "generating medi solution" << endl;
  generateData(solutMedi, a_gridsMedi, a_ebislMedi, a_paramMedi, a_type, a_idir, shiftToCoar);

  //this *is* the coarse soltuion so no shift
  shiftToCoar = false;
  pout() << "generating coar solution" << endl;
  generateData(solutCoar, a_gridsCoar, a_ebislCoar, a_paramCoar, a_type, a_idir, shiftToCoar);

  pout() << "generating coar error from medi and coar solutions" << endl;
  sumFineMinusCoarse(errorCoar,
                     solutCoar, a_gridsCoar, a_ebislCoMe, a_paramCoar,
                     solutMedi, a_gridsMedi, a_ebislMedi, a_paramMedi, a_type);


  Vector<Real> orders;
  Box domCoar = a_paramCoar.coarsestDomain.domainBox();
  if(a_type == VOL_MOM)
    {
  
      Vector<string> names(nvar);
      getMomentNames<CH_EBIS_ORDER>(names, string("m"));
      //the 1 is a verbosity flag.  leave it at one.  trust me.  
      EBArith::compareError(orders, errorMedi, errorCoar,
                            a_gridsMedi, a_gridsCoar, 
                            a_ebislMedi, a_ebislCoar,
                            domCoar, 1, NULL, names, prefix);
    }
  else if((a_type == EB_MOM) || (a_type == EB_NORMAL_MOM))
    {
      
      Vector<string> names(nvar);
      getMomentNames<CH_EBIS_ORDER>(names, string("b"));
      BaseIVFactory<Real> bivrFactMedi(a_ebislMedi);
      BaseIVFactory<Real> bivrFactCoar(a_ebislCoar);
      LevelData<BaseIVFAB<Real> > sparseErrorMedi(a_gridsMedi, nvar, IntVect::Zero, bivrFactMedi);
      LevelData<BaseIVFAB<Real> > sparseErrorCoar(a_gridsCoar, nvar, IntVect::Zero, bivrFactCoar);

      copyDenseToSparse(sparseErrorMedi, errorMedi);
      copyDenseToSparse(sparseErrorCoar, errorCoar);
      EBArith::compareIrregError(orders, sparseErrorMedi, sparseErrorCoar,
                                 a_gridsMedi, a_gridsCoar, 
                                 a_ebislMedi, a_ebislCoar,
                                 domCoar, prefix, names);
    }
  else
    {
      MayDay::Error("bogus type");
    }
  /**/
  pout() << "Outputting moments to file" << endl;
  string solutFileFine =  prefix + string("_Fine.hdf5");
  string solutFileMedi =  prefix + string("_Medi.hdf5");
  string solutFileCoar =  prefix + string("_Coar.hdf5");
  writeEBLevelName(solutFine, solutFileFine);
  writeEBLevelName(solutMedi, solutFileMedi);
  writeEBLevelName(solutCoar, solutFileCoar);

  pout() << "Outputting error to file" << endl;
  string errorFileMedi =  prefix + string("_Error_Medi.hdf5");
  string errorFileCoar =  prefix + string("_Error_Coar.hdf5");
  writeEBLevelName(errorCoar, errorFileCoar);
  writeEBLevelName(errorMedi, errorFileMedi);
  /**/

}
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif
  {

    // Check for an input file
    char* inFile = NULL;

    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    
  GridParameters    paramFine, paramMedi, paramCoar;
  DisjointBoxLayout gridsFine, gridsMedi, gridsCoar;
  //read params from file
  getGridParameters(paramFine, true);
  paramMedi = paramFine;
  paramMedi.coarsen(2);
  paramCoar = paramMedi;
  paramCoar.coarsen(2);


  Vector<int> procs;
  Vector<Box> boxes; 
  domainSplit(paramFine.coarsestDomain, boxes,
              paramFine.maxGridSize   , paramFine.blockFactor);
  LoadBalance(procs, boxes);
  gridsFine = DisjointBoxLayout(boxes, procs, paramFine.coarsestDomain);
  coarsen(gridsMedi, gridsFine, 2);
  coarsen(gridsCoar, gridsMedi, 2);

  pout() << "rct: defining FINE geometry" << endl;
  definePoissonGeometry(paramFine);
  pout() << "saving fine geometry into eblg" << endl;
  EBLevelGrid eblgFine(gridsFine, paramFine.coarsestDomain, 2, Chombo_EBIS::instance());
  pout() << "making CoFi info into eblg" << endl;
  EBLevelGrid eblgCoFi(gridsMedi, paramMedi.coarsestDomain, 2, Chombo_EBIS::instance());

  barrier();
  pout() << "clearing singleton" << endl;
  Chombo_EBIS::instance()->clear();
  pout() << "rct: defining MEDI geometry" << endl;
  definePoissonGeometry(paramMedi);
  pout() << "saving medi geometry into eblg" << endl;
  EBLevelGrid eblgMedi(gridsMedi, paramMedi.coarsestDomain, 2, Chombo_EBIS::instance());
  pout() << "making CoMe info into eblg" << endl;
  EBLevelGrid eblgCoMe(gridsCoar, paramCoar.coarsestDomain, 2, Chombo_EBIS::instance());

  barrier();
  pout() << "clearing singleton" << endl;
  Chombo_EBIS::instance()->clear();
  pout() << "rct: defining Coar geometry" << endl;
  definePoissonGeometry(paramCoar);
  pout() << "saving medi geometry into eblg" << endl;
  EBLevelGrid eblgCoar(gridsCoar, paramCoar.coarsestDomain, 2, Chombo_EBIS::instance());

  EBISLayout ebislFine = eblgFine.getEBISL();
  EBISLayout ebislMedi = eblgMedi.getEBISL();
  EBISLayout ebislCoar = eblgCoar.getEBISL();
  EBISLayout ebislCoFi = eblgCoFi.getEBISL();
  EBISLayout ebislCoMe = eblgCoMe.getEBISL();
  //all thse calls to geometry because we are testing the 
  //accuracy of geometry generation
  //the CoFi and CoMe stuff is because you cannot call refine on
  //coar and medi stuff because as far as it is concerned, it is the 
  //finest level.  They also might have slightly different graphs so this finesses that
  //problem as well
  barrier();

  pout() << "test of volume moments " << endl;
  solutionErrorTest(paramFine,paramMedi,paramCoar,
                    gridsFine,gridsMedi,gridsCoar,
                    ebislFine,ebislMedi,ebislCoar,
                    ebislCoFi, ebislCoMe,
                    VOL_MOM, 0);

  pout() << "test eb area moments" << endl;
  solutionErrorTest(paramFine,paramMedi,paramCoar,
                    gridsFine,gridsMedi,gridsCoar,
                    ebislFine,ebislMedi,ebislCoar,
                    ebislCoFi, ebislCoMe,
                    EB_MOM , 0);

  pout() << "test eb normal moments" << endl;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      solutionErrorTest(paramFine,paramMedi,paramCoar,
                        gridsFine,gridsMedi,gridsCoar,
                        ebislFine,ebislMedi,ebislCoar,
                        ebislCoFi, ebislCoMe,
                        EB_NORMAL_MOM, idir);
    }

  pout() << "clearing singleton " << endl;
  Chombo_EBIS::instance()->clear();
  }
#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

}

