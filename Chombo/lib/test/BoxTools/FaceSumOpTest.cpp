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

#include "BoxIterator.H"
#include "ParmParse.H"
#include "DebugDump.H"
#include "FABView.H"
#include "LoadBalance.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRIO.H"
#include "BoxLayoutData.H"
#include "ReductionOps.H"
#include "ReductionCopier.H"
#include "CellToEdge.H"
#include "BRMeshRefine.H"

#include <fstream>
#include <float.h>

#include "UsingNamespace.H"

#define LARGEINTVAL -999999999
#define LARGEREALVAL 123456789.0

#if 0
void writeFab(const FArrayBox & a_faceFab,
              const int & a_iDir         ,
              const int & a_sumDir       ,
              const int & a_boxId        )
{
  //string prefix("diagnostic/");
  //std::string suffix(".2d.hdf5");

  char name[1000];

  // dummy default argument
  Vector<string> componentNames;

  Real dummyDx = 1.0;
 
  sprintf(name,"output/sumDir%d.idir%d.boxId%d.%dd.hdf5",a_sumDir,a_iDir,a_boxId,SpaceDim);
  writeFABname(&a_faceFab    ,
               name          ,
               componentNames,
               dummyDx       );
}

void writeLevel(const LevelData <FArrayBox> & a_faceFab,
                const int & a_iDir                     ,
                const int & a_sumDir                   )
{
  //  string prefix("diagnostic/");
  //std::string suffix(".2d.hdf5");

  char name[1000];

  // dummy default argument
  //Vector<string> componentNames;

  Real dummyDx = 1.0;
 
  sprintf(name,"output/idir%dsumDir%d.%dd.hdf5",a_iDir,a_sumDir,SpaceDim);
  writeLevelname(&a_faceFab,name);
}
#endif

void initialize(LevelData<FluxBox >     & a_fluxBox,
                const DisjointBoxLayout & a_dbl    ,
                const Real              & a_initVal)
{
  DataIterator dit = a_dbl.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_fluxBox[dit()].setVal(a_initVal);
    }
}

void useWeights(LevelData<FluxBox >      & a_levelDataFluxBox,
                const LevelData<FluxBox >& a_weights         ,
                DisjointBoxLayout        &  a_dbl            ) 
{
  for (int iDir = 0; iDir< SpaceDim; ++iDir)
    {
      DataIterator dit = a_dbl.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox       & fluxBox   = a_levelDataFluxBox[dit()];
          const FluxBox & weightBox = a_weights         [dit()];
          
          FArrayBox&  faceFab = fluxBox    [iDir];
          const FArrayBox&  weightFab = weightBox[iDir];
          
          // weight the doubled faces by 0.5
          faceFab.mult(weightFab);
        }
    }
}

int serialCountFaces(const Vector<Box> a_boxVec, const Box & a_domain)
{
  int faceTotal = 0;
  for (int iDir = 0; iDir< SpaceDim; ++iDir)
    {
      Box domainDir = a_domain;
      domainDir.surroundingNodes(iDir);
      for (BoxIterator bitDomain(domainDir); bitDomain.ok(); ++bitDomain)
        {
          bool foundIt = false;
          for (int iBox = 0; iBox<a_boxVec.size(); ++iBox)
            {
              Box boxDir = a_boxVec[iBox];
              boxDir.surroundingNodes(iDir);
              
              for (BoxIterator bit(boxDir); bit.ok(); ++bit)
                {
                  if (bit() == bitDomain())
                    {
                      foundIt = true;
                      faceTotal += 1;
                      break;
                    }
                }
              if (foundIt)
                {
                  break;
                }
            }
        }
    }
  return faceTotal;
}

void checkAnswer(LevelData<FluxBox >    & a_levelDataFluxBox,
                 const DisjointBoxLayout & a_dbl             ,
                 const int               & a_sumDir          ,
                 const int               &a_totalNumberFaces) 
{
  //pout()<<"Beginnning checkAnswer"<<endl;
  Real faceTotal = 0.0;
  for (int iDir = 0; iDir< SpaceDim; ++iDir)
    {
      //pout()<<"direction = "<<iDir<<endl;
      // FaceFabDataFactory factory(iDir);
      //LevelData<FArrayBox> ldTemp;
      //      ldTemp.define(a_dbl,nComp,ghost,factory);

      DataIterator dit = a_dbl.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox   & fluxBox = a_levelDataFluxBox[dit()];
          const FArrayBox&  faceFab = fluxBox[iDir];
          //pout()<<"box = "<<faceFab.box()<<endl;
          for (BoxIterator bit(faceFab.box()); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              //pout()<<"faceFab("<<iv<<","<<0<<") = "<<faceFab(iv,0)<<endl;
              
              if (iDir == a_sumDir)
                {
                  faceTotal += 0.5*faceFab(iv,0);
                }
              else
                {
                  faceTotal += faceFab(iv,0); 
                }
            }
        }
    }

  DataIterator dit = a_dbl.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox   & fluxBox = a_levelDataFluxBox[dit()];
      const FArrayBox&  faceFab = fluxBox[a_sumDir];
      
      Box iterBox(faceFab.box().smallEnd(),faceFab.box().bigEnd() - BASISV(a_sumDir));
      for (BoxIterator bit(iterBox); bit.ok(); ++bit)
        {
          IntVect loIv = bit();
          IntVect hiIv = loIv + BASISV(a_sumDir);
          if (faceFab(hiIv,0) != faceFab(loIv,0))
            {
              //pout()<<"hiIv = "<<hiIv<<"(faceFab(hiIv,0) = "<<faceFab(hiIv,0)<<endl;
              //pout()<<"loIv = "<<loIv<<"(faceFab(loIv,0) = "<<faceFab(loIv,0)<<endl;
              MayDay::Abort("check answer found inconsistent values");
            }
        }
    }


  // pout()<<"faceTotal = "<<faceTotal<<endl;
  //pout()<<"a_totalNumberfaces = "<<a_totalNumberFaces<<endl;
  if (faceTotal != a_totalNumberFaces)
    {
      if (procID()==0)
        {
          MayDay::Abort("check answer found a different number of faces");
        }
    }
}

Vector<Box> makeBoxes(const ProblemDomain &a_problemDomain)
{
  Vector<Box> boxVec;
 
  int maxSize = 3;
  
  IntVectSet pnd (a_problemDomain.domainBox());
  IntVectSet tags(a_problemDomain.domainBox());
  int totalBufferSize = LARGEINTVAL;
  
  BRMeshRefine meshRefine;
  meshRefine.makeBoxes(boxVec, 
                       tags, 
                       pnd,
                       a_problemDomain,
                       maxSize,
                       totalBufferSize);
 
#if 0 
  for(int iBox = 0; iBox <boxVec.size(); ++iBox)
    {
      pout()<<"iBox["<<iBox<<"] = " <<boxVec[iBox]<<endl;
    }

  Vector<Box> boxVecLess;
  boxVecLess.push_back(boxVec[2]);
  boxVecLess.push_back(boxVec[3]);
  boxVecLess.push_back(boxVec[10]);
  boxVecLess.push_back(boxVec[11]);
  boxVecLess.push_back(boxVec[18]);
  boxVecLess.push_back(boxVec[19]);
  boxVecLess.push_back(boxVec[26]);
  boxVecLess.push_back(boxVec[27]);
#endif

  return boxVec;
}

#if 0
//makeBoxes 3D
Vector<Box> makeBoxes(const Box &a_domain)
{
  Vector<Box> boxVec;
 
  int offsetDir = 0;
   
  Box box0(a_domain);
  box0.setSmall(IntVect(D_DECL(0,0,0)));
  box0.setBig  (IntVect(D_DECL(1,1,1)) - BASISV(offsetDir));
  
  boxVec.push_back(box0);
  
  Box box1(a_domain);
  box1.setSmall(IntVect(D_DECL(1,1,1)));
  box1.setBig  (IntVect(D_DECL(2,2,2)) - BASISV(offsetDir));
  boxVec.push_back(box1);
  
  for(int iBox = 0; iBox <boxVec.size(); ++iBox)
    {
      //pout()<<"iBox["<<iBox<<"] = " <<boxVec[iBox]<<endl;
    }
  
  return boxVec;
}
#endif

#if 0
//2d boxes
 Vector<Box> makeBoxes(const Box &a_domain)
 {
  
  Vector<Box> boxVec;
  int numBox = 5;
  boxVec.resize(numBox);
  
  Box box0(a_domain);
  box0.setSmall(a_domain.smallEnd());
  box0.setBig(IntVect(7,7)) ;
  boxVec[0] = box0;

  Box box1(a_domain);
  box1.setSmall(IntVect(8,8));
  box1.setBig  (IntVect(12,12));
  boxVec[1] = box1;
 
  Box box2(a_domain);
  box2.setSmall(IntVect(8,2));
  box2.setBig  (IntVect(25,5));
  boxVec[2] = box2;
   
  Box box3(a_domain);
  box3.setSmall(IntVect(14,14));
  box3.setBig  (IntVect(17,17));
  boxVec[3] = box3;
  
  Box box4(a_domain);
  box4.setSmall(IntVect(32,32));
  box4.setBig  (IntVect(40,40));
  boxVec[4] = box4;

  for(int iBox = 0;iBox <numBox;++iBox)
    {
      //pout()<<"iBox["<<iBox<<"] = " <<boxVec[iBox]<<endl;
    }
  
  return boxVec;
}
#endif

 Vector<Box> makeReducedBoxes(const Box &a_reducedDomain)
 {
  Vector<Box> reducedBoxVec;
  int numBox = 1;
  reducedBoxVec.resize(numBox);
  reducedBoxVec[0] = a_reducedDomain;
  return reducedBoxVec;
}

int main(int    a_argc,
         char** a_argv)
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
#endif

  // Begin forever present scoping trick
  {
    // test sumOp    
    int nGhost = 1;
    IntVect ghost = IntVect::Unit*nGhost;
    //IntVect nCell((D_DECL6(13,13,18,2,2,2));
   
    IntVect nCell;
    for (int idir = 0; idir < SpaceDim; ++idir)
      {
        nCell[idir] = 13 + 2*idir;
      }
    Real domainLength = nCell[0];
    Real dx = 1.0/(Real)domainLength;
    dx = 1.0;
  
    // ProblemDomain
    IntVect domainLo = IntVect::Zero;
    IntVect domainHi;
    for (int iDir = 0; iDir< SpaceDim;++iDir)
      {
        domainHi[iDir] = nCell[iDir]- 1;
      }
    Box domainBox(domainLo,domainHi);
    ProblemDomain domain(domainBox);

    //Vector<Box> boxVec = makeBoxes(domainBox);
    Vector<Box> boxVec = makeBoxes(domain);

    // check total number of faces
    int totalNumberFaces = serialCountFaces(boxVec,domainBox);
    //pout()<<"totalNumberFaces = "<<totalNumberFaces<<endl;
  
    int numBoxes = boxVec.size();
    Vector<int> procAssign;
    procAssign.resize(numBoxes);
    int numProcessor = numProc();
    LoadBalance(procAssign,boxVec,numProcessor);

    
    //pout()<<"procAssign = "<<procAssign<<endl;
  
    // repeat for the reduced box layout
    Box reducedDomainBox(domainBox);
    int sumDir = 0;

    CH_assert(sumDir < SpaceDim && sumDir >= 0);
  
    reducedDomainBox.setSmall(sumDir,0);
    reducedDomainBox.setBig  (sumDir,0);
    ProblemDomain reducedDomain                 (reducedDomainBox);
    Vector<Box> reducedBoxVec = makeReducedBoxes(reducedDomainBox);
  
    int reducedNumBoxes = reducedBoxVec.size();
    Vector<int> reducedProcAssign;
    reducedProcAssign.resize(reducedNumBoxes);
    LoadBalance(reducedProcAssign,reducedBoxVec,1);
      
    // dbl
    DisjointBoxLayout dbl       (boxVec            , procAssign       , domain       );
    DisjointBoxLayout reducedDbl(reducedBoxVec     , reducedProcAssign, reducedDomain);
        
    int nComp = 1;
    LevelData<FluxBox> faceData       (dbl       , nComp);
    LevelData<FluxBox> weights        (dbl       , nComp);
    LevelData<FluxBox> reducedFaceData(reducedDbl, nComp);
    
    // initialize data
    Real testVal = 1.0;

    initialize(faceData,dbl,testVal);

    // make weights and modify the faces to prevent double counting
    computeFaceReductionWeights(weights);
    useWeights(faceData,
               weights ,
               dbl     ); 
  
    // initialize reduced data
    Real zeroVal = 0.0;
    initialize(reducedFaceData,reducedDbl,zeroVal);
    
    // face summing operator
    //pout()<<"sumDir = "<<sumDir<<endl;
    FaceSumOp FaceSumOp(sumDir);
    FaceSumOp.setScale(dx);

    // make the copier that calls op
    ReductionCopier reduceCopier(dbl          , 
                                 reducedDbl   , 
                                 domain       ,
                                 IntVect::Zero, 
                                 sumDir       );

    //reduceCopier.print();
  
    //use the copier, which calls op
    faceData.copyTo(faceData.interval()       , 
                    reducedFaceData           , 
                    reducedFaceData.interval(),
                    reduceCopier              , 
                    FaceSumOp                    );

    //check answer
    checkAnswer(reducedFaceData,reducedDbl,sumDir,totalNumberFaces);
#if 0
    if (numProcessor == 1)
      {
        for (int idir = 0; idir< SpaceDim; ++idir)
          {
            DataIterator dit = reducedDbl.dataIterator();
            int boxId;
            for (dit.begin(); dit.ok(); ++dit)
              {
                FArrayBox& reduceFab = reducedFaceData[dit()][idir];
                boxId = dit().datInd();
                writeFab(reduceFab,
                         idir     ,
                         sumDir   ,
                         boxId   );
            
              }
    
          }
      }
#endif
  }// End scoping trick 
  
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  
  return 0;
}
  


#if 0    
// faceFab factory
int dir = 1;
FaceFabDataFactory factory(dir);
   
// construct vector of boxed
Box box1(IntVect::Zero,IntVect::Unit*3);
Box box2(IntVect(0,4),IntVect(3,6));
Vector <Box> boxVec;
boxVec.resize(2);
boxVec[0] = box1;
//pout() << "box1 = "<<box1<<endl;

boxVec[1] = box2;
//pout() << "box2 = "<<box2<<endl;
// assign proc ids
Vector<int> procId;
    
LoadBalance(procId,boxVec);

// construct disjoint box layout
DisjointBoxLayout dbl(boxVec,procId);

// components
int nComp = 1;

// ghost
int nGhost = 0;
IntVect ghost = IntVect::Unit*nGhost;
      
// define levelData of faceFab
LevelData<FArrayBox> ldFaceFab;
ldFaceFab.define(dbl,nComp,ghost,factory);
       
int sumDir = 0;
FaceSumOp FaceSumOp(sumDir);
FaceSumOp.setDir(sumDir);
#endif
