#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LOCALFUNCTIONS_H_
#define _LOCALFUNCTIONS_H_
#include  <iostream>
#include  <cstdio>
#include "LevelOp.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "PoissonBC.H"
#include "ParmParse.H"
#include "AMRSolver.H"
#include "PoissonOp.H"
#include "Vector.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"

typedef LevelData<FArrayBox> LDFB;
typedef DisjointBoxLayout    DBL;

// global IntVect ... ?
const IntVect unitVector(IntVect::Unit);

void newTagCells(const Vector<LDFB*>&   rhs,
                 const Vector<LDFB*>&   phi,
                 const int numLevels,
                 //const AMRSolver& amrSolver,
                 const Vector<DBL>& DBLVector,
                 const Vector<Box>& boundingBoxVector,
                 const Vector<Real>&   DxVector,
                 const Vector<int>&    refinementRatioVector,
                 const DomainGhostBC& domghostbc,
                 const int baseLevel,
                 const bool verbose,
                 const int maxIterations,
                 const int numberVcyclesBottom,
                 const Real threshold,
                 Vector<IntVectSet>&   tags);


Box createBox(const int *Ncells, const int dimensions);


void createHOG(Vector<LDFB*>& second, const int numLevels,
               const Vector<DBL>& DBLVector);
void multiplyHOGbyConstant(Vector<LDFB*>& hog,
                           const int numLevels,
                           const Real k1);
void plusEqualHOG(Vector<LDFB*>& hog1,
                  const int numLevels,
                  const Vector<LDFB*>& hog2);
void compareTwoHOGs(const Vector<LDFB*>& grid1, const Vector<LDFB*>& grid2,
                    const int numLevels, const bool returnRawDiff,
                    Vector<LDFB*>& absDifference);

void subtractTwoHOGs(const Vector<LDFB*>& grid1, const Vector<LDFB*>& grid2,
                     const int numLevels, Vector<LDFB*>& difference);

void copyFirstHOGintoSecond(const Vector<LDFB*>& first, Vector<LDFB*>& second,
                            const int numLevels);


void averageDownHOG(Vector<LDFB*>& grid,
                    const Vector<int>& refinementRatioVector,
                    const int numLevels);


void makeHDF5fromIVS(const IntVectSet& tags,
                            const Box& boundingBox,
                            const string& name);

void fillColorFAB(LDFB* colorLDFB,
                  const IntVectSet& coarseFineSkinSet,
                  const IntVectSet& fineCoarseSkinSet,
                  const IntVectSet& problemBoundarySkinSet);

void coarseFineBoundaryCount(const DisjointBoxLayout& dbl,
                             const Box& boundingBoxLevel0,
                             IntVectSet& coarseFineSkinSet,
                             int& numberCoarseFineBoundaryFaces);

void fineCoarseBoundaryCount(const DisjointBoxLayout& dbl,
                             const Box& boundingBoxLevel0,
                             IntVectSet& fineCoarseSkinSet);


void problemBoundaryCount(const DisjointBoxLayout& dbl,
                          const Box& boundingBoxLevel0,
                          IntVectSet& problemBoundarySkinSet,
                          int& numberProblemBoundaryFaces);

int cellCount(const DisjointBoxLayout& dbl);

void allBoundaryCount(const DisjointBoxLayout& dbl,
                      const Box& boundingBoxLevel0,
                      IntVectSet& coarseFineSkinSet,
                      IntVectSet& fineCoarseSkinSet,
                      IntVectSet& problemBoundarySkinSet,
                      int& numberCoarseFineBoundaryFaces,
                      int& numberProblemBoundaryFaces,
                      int& numberCells);


///
//   Prints inputs to standard out.  For use by me.
void prettyPrintLDFB(const LevelData<FArrayBox>*  memLDFB_Ptr,
                    const string& label);

void prettyPrintLDFB(const LevelData<FArrayBox>*  memLDFB_Ptr);



//  set domain boundary conditions from input file
int setDomainBC(DomainGhostBC& domghostbc,
                const int* ibclo,
                const int* ibchi,
                const bool verbose);

// out output one data set to specified filename
int outputHDF5(const Vector<LevelData<FArrayBox>* >& vectData,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<Box>& vectDomain,
               const Vector<int>& vectRatio,
               const int numlevels,
               const string& filename);

  //    Output phi and rhs to file
int outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
               const Vector<LevelData<FArrayBox>* >& vectRhs,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<Box>& domain,
               const Vector<int>& vectRatio,
               const int numlevels,
               const bool verbose);


//int regrid(Vector<LevelData<FArrayBox>*>& vectData,
//         Vector<DisjointBoxLayout>&     vectGrids,
//         const Vector<Box>&             vectDomain,
//         const Vector<int>&             vectRefRatio,
//         const int  blockFactor,
//         const int  maxBoxSize,
//         const int  bufferSize,
//         const Real fillRatio,
//         const Real refineThreshold,
//         const int maxLevel,
//         const int base_level,
//         int& numlevels);



#endif
