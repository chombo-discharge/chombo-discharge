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
#include  <iostream>
#include "parstream.H"
#include "CONSTANTS.H"
#include "FABView.H"
#include "DebugDump.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "ParmParse.H"
#include "BRMeshRefine.H" // for domainSplit
#include "LoadBalance.H"
#include "MultiBlockLevelCopierCenter.H"
#include "MultiBlockLevelCopierAverage.H"
#include "CylinderEquiangularCS.H"
#include "CubedSphere2DCS.H"
#include "CubedSphereShellCS.H"
#include "DoubleCartesianCS.H"
#include "DoubleCartesianRotateCS.H"
#include "CubedSphereSolidCS.H"
#include "CylindricalHarmonic.H"
#include "CylindricalHarmonicGrad.H"
#include "SphericalHarmonic.H"
#include "SphericalHarmonicGrad.H"
#include "Trig2.H"
#include "Trig2Grad.H"
#include "FourthOrderUtil.H"
// for parallel debugging
#include "CH_Attach.H"
using std::cerr;
using std::endl;
#include "CH_Timer.H"

#include "UsingNamespace.H"

/**
 * Test MultiBlockLevelCopier.
 * petermc, 13 Jun 2012
 */

enum MapCode
  {
    DOUBLECARTESIAN,
    DOUBLECARTESIANROTATE,
    CYLINDEREQUIANGULAR,
    CUBEDSPHERE2D,
    CUBEDSPHERESHELL,
    CUBEDSPHERESOLID
  };
const std::string mapStr[] =
  {
    "DoubleCartesianCS",
    "DoubleCartesianRotateCS",
    "CylinderEquiangularCS",
    "CubedSphere2DCS",
    "CubedSphereShellCS",
    "CubedSphereSolidCS"
  };

enum CenteringCode
  {
    CENTER,
    AVERAGE
  };
const std::string centeringStr[] =
  {
    "Center",
    "Average"
  };

// ---------------------------------------------------------
#ifdef CH_MPI
void reduceReal(Real&           a_val,
                const MPI_Op&   a_mpiOp)
{
  // Kludge to return NaN if a_val is NaN on any processor.
  unsigned short int isNaNlocal = (a_val == a_val) ? 0 : 1;
  unsigned short int isNaNglobal;
  int redMPI = MPI_Allreduce(&isNaNlocal, &isNaNglobal,
                             1, MPI_UNSIGNED_SHORT,
                             MPI_MAX, Chombo_MPI::comm);
  if (redMPI != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on reduceReal");
    }
  if (isNaNglobal == 1)
    {
      a_val = 0./0.;
      return;
    }

  Real recv;
  int resultMPI = MPI_Allreduce(&a_val, &recv,
                                1, MPI_CH_REAL,
                                a_mpiOp, Chombo_MPI::comm);

  if (resultMPI != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on reduceReal");
    }
  a_val = recv;
}
#endif

// ---------------------------------------------------------
Real vectorMax(const Vector<Real>& a_vec)
{
  Real vecMax = 0.;
  for (int i = 0; i < a_vec.size(); i++)
    {
      if ( (a_vec[i] != a_vec[i]) || // also known as NaN
           (a_vec[i] > vecMax) )
        {
          vecMax = a_vec[i];
        }
    }
  return vecMax;
}

// ---------------------------------------------------------
RealVect vectorVecMax(const Vector<RealVect>& a_vec)
{
  RealVect vecMax = RealVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int i = 0; i < a_vec.size(); i++)
        {
          if ( (a_vec[i][idir] != a_vec[i][idir]) || // also known as NaN
               (a_vec[i][idir] > vecMax[idir]) )
            {
              vecMax[idir] = a_vec[i][idir];
            }
        }
    }
  return vecMax;
}

// ---------------------------------------------------------
Vector<Real> blockMaxima(const BoxLayoutData<FArrayBox>& a_data,
                         const MultiBlockCoordSys* a_coordSysPtr)
{
  CH_assert(a_data.nComp() == 1);
  int nblocks = a_coordSysPtr->numBlocks();
  Vector<Real> blockMax(nblocks, 0.);
  const BoxLayout& grids = a_data.boxLayout();
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& bx = grids[dit];
      int blockNum = a_coordSysPtr->whichBlockOverlap(bx);
      const FArrayBox& dataFab = a_data[dit];
      // Norms all on bx:  ignore external stuff.
      Real dataFabMax = dataFab.norm(bx, 0);
      if (dataFabMax > blockMax[blockNum])
        blockMax[blockNum] = dataFabMax;
    }
# ifdef CH_MPI
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      reduceReal(blockMax[iblock], MPI_MAX);
    }
# endif
  return blockMax;
}

// ---------------------------------------------------------
Vector<RealVect> blockVecMaxima(const BoxLayoutData<FArrayBox>& a_data,
                                const MultiBlockCoordSys* a_coordSysPtr)
{
  CH_assert(a_data.nComp() == SpaceDim);
  int nblocks = a_coordSysPtr->numBlocks();
  Vector<RealVect> blockVecMax(nblocks, RealVect::Zero);
  const BoxLayout& grids = a_data.boxLayout();
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& bx = grids[dit];
      int blockNum = a_coordSysPtr->whichBlockOverlap(bx);
      const FArrayBox& dataFab = a_data[dit];
      // Max-norms all on bx:  ignore external stuff.
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real dataFabMax = dataFab.norm(bx, 0, idir);
          if (dataFabMax > blockVecMax[blockNum][idir])
            blockVecMax[blockNum][idir] = dataFabMax;
        }
    }
# ifdef CH_MPI
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          reduceReal(blockVecMax[iblock][idir], MPI_MAX);
        }
    }
# endif
  return blockVecMax;
}

template <class T>
const std::string paddedStr(T r, std::streamsize n)
{
  std::ostringstream os;
  os << std::setw(n) << r;
  return os.str();
}

int main(int argc, char* argv[])
{
  int status = 0; // number of errors detected.
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  //scoping trick
  {
    // test parameters
    const int testOrder = 4; // expected order of convergence

    // A test is considered as a failure
    //  if its convergence rate is smaller than below.
    const Real targetConvergeRate = testOrder*0.9;

    // Maximum dimension of a grid
    int maxGridSize = 32;
    // Minimum dimension of a grid
    int blockFactor = 4;

    // Do vector field too?
    bool doVector = true;

    // number of ghost cells being interpolated
    int numGhost = 4;

    Vector<CenteringCode> centerings;
    centerings.push_back(AVERAGE);
    centerings.push_back(CENTER);
    int nCenterings = centerings.size();

    Vector<MapCode> maps;
    if (SpaceDim >= 2)
      { // MultiBlockLevelCopier works only if SpaceDim >= 2.
        // Too many bad calls to mappedCoord with CUBEDSPHERESOLID.
        // if (SpaceDim == 3) maps.push_back(CUBEDSPHERESOLID);
        if (SpaceDim <= 3) maps.push_back(CUBEDSPHERESHELL);
        if (SpaceDim <= 3) maps.push_back(CUBEDSPHERE2D);
        if (SpaceDim <= 3) maps.push_back(CYLINDEREQUIANGULAR);
        maps.push_back(DOUBLECARTESIAN);
        maps.push_back(DOUBLECARTESIANROTATE);
      }
    int nMaps = maps.size();

    // The big nested loop for different cases.
    for (int indCen = 0; indCen < nCenterings; indCen++)
      {
        CenteringCode codeCen = centerings[indCen];
        RefCountedPtr<MultiBlockLevelCopier> mblcPtr;
        // MultiBlockLevelCopier* mblcPtr;
        switch (codeCen)
          {
          case CENTER :
            {
              mblcPtr = RefCountedPtr<MultiBlockLevelCopierCenter>(new MultiBlockLevelCopierCenter());
              // mblcPtr = new MultiBlockLevelCopierCenter();
              break;
            }
          case AVERAGE :
            {
              mblcPtr = RefCountedPtr<MultiBlockLevelCopierAverage>(new MultiBlockLevelCopierAverage());
              // mblcPtr = new MultiBlockLevelCopierAverage();
              break;
            }
          default : 
            {
              MayDay::Error("unknown centering");
            }
          }
        bool atCellCenters = (codeCen == CENTER);

        int extraCentering = 0;
        if (!atCellCenters)
          { // need an additional layer to use for averaging
            extraCentering++;
          }
        IntVect extraVect = extraCentering * IntVect::Unit;

        for (int indMap = 0; indMap < nMaps; indMap++)
          {
            MapCode codeMap = maps[indMap];

            RefCountedPtr<MultiBlockCoordSysFactory> coordSysFactPtr;

            IntVect basicDomainLo, basicDomainHi;
            // Once we get int domainLength, we will set:
            // levelDomainLo = domainLength * basicDomainLo;
            // levelDomainHi = domainLength * basicDomainHi - IntVect::Unit;

            RefCountedPtr<FunctionOnSpace> funPtr;
            RefCountedPtr<VectorFunctionOnSpace> vectorFunPtr;
            
            switch (codeMap)
              {
              case DOUBLECARTESIAN :
                {
                  RefCountedPtr<DoubleCartesianCSFactory>
                    doubleCartesianCSFactPtr =
                    RefCountedPtr<DoubleCartesianCSFactory>(new DoubleCartesianCSFactory); 
                  coordSysFactPtr = doubleCartesianCSFactPtr;

                  basicDomainLo = IntVect::Zero;
                  basicDomainHi = 3 * IntVect::Unit;

                  funPtr = RefCountedPtr<FunctionOnSpace>
                    (new Trig2());
                  if (doVector)
                    vectorFunPtr = RefCountedPtr<VectorFunctionOnSpace>
                      (new Trig2Grad());
                  break;
                }
              case DOUBLECARTESIANROTATE :
                {
                  RefCountedPtr<DoubleCartesianRotateCSFactory>
                    doubleCartesianRotateCSFactPtr =
                    RefCountedPtr<DoubleCartesianRotateCSFactory>(new DoubleCartesianRotateCSFactory); 
                  coordSysFactPtr = doubleCartesianRotateCSFactPtr;

                  basicDomainLo = IntVect::Zero;
                  basicDomainHi = 3 * IntVect::Unit;

                  funPtr = RefCountedPtr<FunctionOnSpace>
                    (new Trig2());
                  if (doVector)
                    vectorFunPtr = RefCountedPtr<VectorFunctionOnSpace>
                      (new Trig2Grad());
                  break;
                }
              case CYLINDEREQUIANGULAR :
                {
                  RealVect centerPoint = RealVect::Zero;
                  RealVect centralRectSize = RealVect::Unit;
                  Real outerRadius = 1.5;
                  RefCountedPtr<CylinderEquiangularCSFactory>
                    cylinderCSFactPtr =
                    RefCountedPtr<CylinderEquiangularCSFactory>(new CylinderEquiangularCSFactory);
                  cylinderCSFactPtr->setCenterPoint(centerPoint);
                  cylinderCSFactPtr->setCentralRectSize(centralRectSize);
                  cylinderCSFactPtr->setOuterRadius(outerRadius);
                  coordSysFactPtr = cylinderCSFactPtr;

                  basicDomainLo = IntVect(D_DECL(-2, -2, 0));
                  basicDomainHi = IntVect(D_DECL(3, 3, 1));

                  int cylinderN = 4;
                  int cylinderK = 1;
                  funPtr = RefCountedPtr<FunctionOnSpace>
                    (new CylindricalHarmonic(cylinderN, cylinderK));
                  if (doVector)
                    vectorFunPtr = RefCountedPtr<VectorFunctionOnSpace>
                      (new CylindricalHarmonicGrad(cylinderN, cylinderK));
                  break;
                }
              case CUBEDSPHERE2D :
                {
                  Real rMin = 0.03;
                  Real rMax = 1.;
                  RefCountedPtr<CubedSphere2DCSFactory>
                    cubedSphere2DCSFactPtr =
                    RefCountedPtr<CubedSphere2DCSFactory>(new CubedSphere2DCSFactory);
                  cubedSphere2DCSFactPtr->setRadii(rMin, rMax);
                  coordSysFactPtr = cubedSphere2DCSFactPtr;

                  basicDomainLo = IntVect::Zero;
                  basicDomainHi = IntVect(D_DECL(11, 1, 1));

                  int legendreL = 4;
                  int legendreM = 1;
                  funPtr = RefCountedPtr<FunctionOnSpace>
                    (new SphericalHarmonic(legendreL, legendreM));
                  if (doVector)
                    vectorFunPtr = RefCountedPtr<VectorFunctionOnSpace>
                      (new SphericalHarmonicGrad(legendreL, legendreM));
                  break;
                }
              case CUBEDSPHERESHELL :
                {
                  RefCountedPtr<CubedSphereShellCSFactory>
                    cubedSphereShellCSFactPtr =
                    RefCountedPtr<CubedSphereShellCSFactory>(new CubedSphereShellCSFactory);
                  cubedSphereShellCSFactPtr->setHeight(10.);
                  cubedSphereShellCSFactPtr->setRadius(6378.);
                  coordSysFactPtr = cubedSphereShellCSFactPtr;

                  basicDomainLo = IntVect::Zero;
                  basicDomainHi = IntVect(D_DECL(11, 1, 1));

                  int legendreL = 4;
                  int legendreM = 1;
                  funPtr = RefCountedPtr<FunctionOnSpace>
                    (new SphericalHarmonic(legendreL, legendreM));
                  if (doVector)
                    vectorFunPtr = RefCountedPtr<VectorFunctionOnSpace>
                      (new SphericalHarmonicGrad(legendreL, legendreM));
                  break;
                }
              case CUBEDSPHERESOLID :
                {
                  RefCountedPtr<CubedSphereSolidCSFactory>
                    cubedSphereSolidCSFactPtr =
                    RefCountedPtr<CubedSphereSolidCSFactory>(new CubedSphereSolidCSFactory);
                  cubedSphereSolidCSFactPtr->setR0(1.);
                  cubedSphereSolidCSFactPtr->setR1(3.);
                  coordSysFactPtr = cubedSphereSolidCSFactPtr;

                  basicDomainLo = -2 * IntVect::Unit;
                  basicDomainHi = 3 * IntVect::Unit;

                  int legendreL = 4;
                  int legendreM = 1;
                  // funPtr = RefCountedPtr<FunctionOnSpace>
                  //   (new SphericalRadius());
                  funPtr = RefCountedPtr<FunctionOnSpace>
                    (new SphericalHarmonic(legendreL, legendreM));
                  if (doVector)
                    vectorFunPtr = RefCountedPtr<VectorFunctionOnSpace>
                      (new SphericalHarmonicGrad(legendreL, legendreM));
                  break;
                }
              default :
                {
                  MayDay::Error("unknown mapping");
                }
              }

            Vector<Interval> fixedDimsAll;
            if ( (SpaceDim == 3) &&
                 ( (codeMap == CYLINDEREQUIANGULAR) ||
                   (codeMap == CUBEDSPHERE2D) ||
                   (codeMap == CUBEDSPHERESHELL) ) )
              { // This is the only coordsys we can try it on here.
                fixedDimsAll.push_back(Interval(SpaceDim-1, SpaceDim-1));
              }
            // Always include empty interval, meaning NO fixed dimensions.
            fixedDimsAll.push_back(Interval());
            
            for (int ifixed = 0; ifixed < fixedDimsAll.size(); ifixed++)
              {
                const Interval& fixedDims = fixedDimsAll[ifixed];
                int nfixed = fixedDims.size();
                Vector<int> fixedPt(nfixed);
                // Set the components in all fixed dimensions to 5.
                fixedPt.assign(5);

                IntVect ghostVect = numGhost * IntVect::Unit;
                // No ghost cells in fixed dimensions.
                for (int idir = fixedDims.begin(); idir <= fixedDims.end(); idir++)
                  {
                    ghostVect[idir] = 0;
                  }

                int domainLengthMin = 16;
                if  (codeMap == CYLINDEREQUIANGULAR)
                  { // 4 ghost cells takes too much of a bite out of length 16
                    domainLengthMin = 32;
                  }

#ifdef CH_USE_FLOAT
                // Single precision doesn't get us very far.
                int domainLengthMax = 32;
                if (codeMap == CYLINDEREQUIANGULAR)
                  { // In this CS, 32 to 64 not good enough in single precision.
                    domainLengthMax = -1;
                  }
                else if ( (codeMap == CUBEDSPHERE2D) ||
                          (codeMap == CUBEDSPHERESHELL) )
                  { // In this CS, we can get finer.
                    domainLengthMax = 64;
                  }
#endif

#ifdef CH_USE_DOUBLE
#if CH_SPACEDIM >= 3
                int domainLengthMax = 64;
                if (nfixed > 0) domainLengthMax = 128;
#else
                int domainLengthMax = 512;
#endif
#endif
                
                Vector<int> domainLengths;
                int len = domainLengthMin;
                while (len <= domainLengthMax)
                  {
                    domainLengths.push_back(len);
                    len *= 2;
                  }
                int nGrids = domainLengths.size();

                if (nGrids > 0)
                  {
                    pout() << endl
                           << "Testing MultiBlockLevelCopier"
                           << centeringStr[codeCen]
                           << "; DIM=" << SpaceDim
                           << "; " << mapStr[codeMap];
                    if (nfixed > 0)
                      {
                        pout() << " fixing " << fixedDims.begin()
                               << ":" << fixedDims.end();
                      }
                    pout() << endl;

                    pout() << "size  "
                           << "scalar     rate";
                    if (doVector)
                      {
                        for (int idir = 0; idir < SpaceDim; idir++)
                          {
                            pout() << "  vector " << idir << "   rate";
                          }
                      }
                    pout() << endl;
                  }
                
                bool passTest = true;
                Vector<Real> diffNormsCoarser;
                for (int iGrid = 0; iGrid < nGrids; iGrid++)
                  {
                    Vector<Real> diffNorms;
                    
                    int domainLength = domainLengths[iGrid];
                    
                    IntVect levelDomainLo =
                      domainLength * basicDomainLo;
                    IntVect levelDomainHi =
                      domainLength * basicDomainHi - IntVect::Unit;
                    Box levelDomainBox(levelDomainLo, levelDomainHi);
                    ProblemDomain levelDomain(levelDomainBox);
                    RealVect dx = (1.0 / Real(domainLength)) * RealVect::Unit;
                    MultiBlockCoordSys* coordSysPtr =
                      coordSysFactPtr->getCoordSys(levelDomain, dx);
                    
                    int realDim = coordSysPtr->realDim();
                    
                    const Vector<Box>& blockBoxes = coordSysPtr->mappingBlocks();
                    int nblocks = blockBoxes.size();
                    
                    Vector<Box> allBoxes;
                    for (int iblock = 0; iblock < nblocks; iblock++)
                      {
                        Vector<Box> thisBlockBoxes;
                        domainSplit(blockBoxes[iblock], thisBlockBoxes,
                                    maxGridSize, blockFactor);
                        allBoxes.append(thisBlockBoxes);
                      }
                    mortonOrdering(allBoxes);
                    
                    Vector<int> allProcs(allBoxes.size());
                    LoadBalance(allProcs, allBoxes);
                    DisjointBoxLayout srcLayout(allBoxes, allProcs);
                    
                    // Set dstLayout to be srcLayout expanded by ghostVect.
                    Vector<Box> allExpandedBoxes(allBoxes.size());
                    for (int i = 0; i < allBoxes.size(); i++)
                      {
                        Box bx(allBoxes[i]);
                        bx.grow(ghostVect);
                        int blockNum = coordSysPtr->whichBlockOverlap(bx);
                        coordSysPtr->keepInDomain(bx, blockNum);
                        allExpandedBoxes[i] = bx;
                      }
                    BoxLayout dstLayout(allExpandedBoxes, allProcs);
                    
                    // petermc, 19 Nov 2012: WAS extraCentering, not 0.
                    MultiBlockLevelGeom geom(coordSysPtr, srcLayout, 0,
                                             fixedDims, fixedPt);
                    mblcPtr->define(&geom, dstLayout, testOrder);
                    
                    ///////////////////////////////////////////
                    // SCALAR FIELD
                    ///////////////////////////////////////////
                    
                    /*
                      Define data holders for scalar.
                    */
                    int ncomp = 1;
                    Interval intvlCalc(0, ncomp-1);
                    Interval intvlExact(ncomp, 2*ncomp-1);
                    Interval intvlDiff(2*ncomp, 3*ncomp-1);
                    
                    LevelData<FArrayBox> src(srcLayout, ncomp);
                    
                    // Destination data holder.
                    BoxLayoutData<FArrayBox> dstAll(dstLayout, 3*ncomp);
                    BoxLayoutData<FArrayBox> dstCalc, dstExact, dstDiff;
                    AliasDataFactory<FArrayBox> factoryCalc(&dstAll, intvlCalc);
                    AliasDataFactory<FArrayBox> factoryExact(&dstAll, intvlExact);
                    AliasDataFactory<FArrayBox> factoryDiff(&dstAll, intvlDiff);
                    dstCalc.define(dstLayout, intvlCalc.size(), factoryCalc);
                    dstExact.define(dstLayout, intvlExact.size(), factoryExact);
                    dstDiff.define(dstLayout, intvlDiff.size(), factoryDiff);
                    
                    /*
                      Find cell centers.
                    */
                    // MultiBlockLevelGeom geom(coordSysPtr, srcLayout, extraCentering);
                    // const LevelData<FArrayBox>& cellCenters = geom.cellCenters();
                    LevelData<FArrayBox> cellCenters(srcLayout, realDim, extraVect);
                    
                    // Set cellCenters.
                    for (DataIterator dit = srcLayout.dataIterator(); dit.ok(); ++dit)
                      {
                        CH_TIME("computing cellCenters");
                        
                        FArrayBox& cellCentersFab = cellCenters[dit];
                        const Box& baseBox = srcLayout[dit];
                        int thisBlock = coordSysPtr->whichBlock(baseBox);
                        const NewCoordSys* thisCoordSysPtr =
                          coordSysPtr->getCoordSys(thisBlock);
                        Box cellCentersBox = grow(baseBox, extraVect);
                        FArrayBox xiCentersFab(cellCentersBox, SpaceDim);
                        thisCoordSysPtr->getCenterMappedCoordinates(xiCentersFab,
                                                                    cellCentersBox);
                        thisCoordSysPtr->realCoord(cellCentersFab,
                                                   xiCentersFab,
                                                   cellCentersBox);
                      }
                    
                    /*
                      Fill in data holders for scalar.
                    */
                    
                    // Set src to exact value of function on srcLayout.
                    for (DataIterator ditSrc = srcLayout.dataIterator(); ditSrc.ok(); ++ditSrc)
                      {
                        CH_TIME("computing srcFab");
                        const Box& bxBase = srcLayout[ditSrc];
                        const FArrayBox& cellCentersFab = cellCenters[ditSrc];
                        FArrayBox& srcFab = src[ditSrc];
                        if (atCellCenters)
                          {
                            funPtr->setFunctionFromPhysical(srcFab,
                                                            bxBase,
                                                            cellCentersFab);
                          }
                        else
                          { // Get average of function on each cell.
                            Box bx1 = grow(bxBase, 1);
                            FArrayBox srcMoreFab(bx1, ncomp);
                            funPtr->setFunctionFromPhysical(srcMoreFab,
                                                            bx1,
                                                            cellCentersFab);
                            fourthOrderAverageCell(srcMoreFab);
                            srcFab.copy(srcMoreFab, bxBase);
                          }
                      }
                    
                    src.copyTo(dstCalc);
                    mblcPtr->copyTo(src, dstCalc);
                    
                    /*
                      Got results in dstCalc.
                      Now find dstExact and dstDiff = dstCalc - dstExact.
                    */
                    const LayoutData<int>& blockNumber = geom.block();
                    for (DataIterator ditDst = dstLayout.dataIterator(); ditDst.ok(); ++ditDst)
                      {
                        const Box& bx = dstLayout[ditDst];
                        
                        const FArrayBox& calcFab = dstCalc[ditDst];
                        FArrayBox& exactFab = dstExact[ditDst];
                        FArrayBox& diffFab = dstDiff[ditDst];
                        
                        // Find cellCentersFab on bxExpanded.
                        Box bxExpanded = grow(bx, extraVect);
                        FArrayBox cellCentersFab(bxExpanded, realDim);
                        FArrayBox xiCentersFab(bxExpanded, SpaceDim);
                        int blockNum = coordSysPtr->whichBlockOverlap(bx);
                        const NewCoordSys* blockCoordSysPtr =
                          coordSysPtr->getCoordSys(blockNum);
                        blockCoordSysPtr->getCenterMappedCoordinates(xiCentersFab, bxExpanded);
                        blockCoordSysPtr->realCoord(cellCentersFab, xiCentersFab, bxExpanded);
                        
                        // Find exactFab.
                        if (atCellCenters)
                          {
                            funPtr->setFunctionFromPhysical(exactFab,
                                                            bxExpanded, // == bx
                                                            cellCentersFab);
                          }
                        else
                          { // Get average of function on each cell.
                            FArrayBox exactExpandedFab(bxExpanded, ncomp);
                            funPtr->setFunctionFromPhysical(exactExpandedFab,
                                                            bxExpanded, // == grow(bx, 1)
                                                            cellCentersFab);
                            fourthOrderAverageCell(exactExpandedFab);
                            exactFab.copy(exactExpandedFab, bx);
                          }
                        
                        diffFab.copy(calcFab);
                        diffFab.minus(exactFab);
                      }
                    
                    Vector<Real> diffBlockMax =
                      blockMaxima(dstDiff, coordSysPtr);
                    Vector<Real> exactBlockMax =
                      blockMaxima(dstExact, coordSysPtr);
                    Real diffAllMax = vectorMax(diffBlockMax);
                    Real exactAllMax = vectorMax(exactBlockMax);
                    diffNorms.push_back(diffAllMax / exactAllMax);
                    
                    if (doVector)
                      {
                        ///////////////////////////////////////////
                        // VECTOR FIELD
                        ///////////////////////////////////////////
                        
                        mblcPtr->defineVector();

                        /*
                          Define data holders for vector.
                        */

                        Interval intvlVectorCalc(0, SpaceDim-1);
                        Interval intvlVectorExact(SpaceDim, 2*SpaceDim-1);
                        Interval intvlVectorRealExact(SpaceDim, SpaceDim+realDim-1);
                        Interval intvlVectorDiff(2*SpaceDim, 3*SpaceDim-1);

                        // Recall realDim >= SpaceDim.
                        LevelData<FArrayBox> srcVectorReal(srcLayout, realDim);
                        LevelData<FArrayBox> srcVector;
                        Interval intvlVectorReal(0, realDim-1);
                        Interval intvlVector(0, SpaceDim-1);
                        aliasLevelData(srcVector, &srcVectorReal, intvlVector);

                        BoxLayoutData<FArrayBox> dstVectorAll(dstLayout, 3*SpaceDim);
                        BoxLayoutData<FArrayBox> dstVectorCalc, dstVectorRealExact, dstVectorExact, dstVectorDiff;
                        AliasDataFactory<FArrayBox> factoryVectorCalc(&dstVectorAll, intvlVectorCalc);
                        AliasDataFactory<FArrayBox> factoryVectorExact(&dstVectorAll, intvlVectorExact);
                        AliasDataFactory<FArrayBox> factoryVectorRealExact(&dstVectorAll, intvlVectorRealExact);
                        AliasDataFactory<FArrayBox> factoryVectorDiff(&dstVectorAll, intvlVectorDiff);
                        dstVectorCalc.define(dstLayout, intvlVectorCalc.size(), factoryVectorCalc);
                        dstVectorExact.define(dstLayout, intvlVectorExact.size(), factoryVectorExact);
                        dstVectorRealExact.define(dstLayout, intvlVectorRealExact.size(), factoryVectorRealExact);
                        dstVectorDiff.define(dstLayout, intvlVectorDiff.size(), factoryVectorDiff);

                        /*
                          Fill in srcVectorReal with vector function values at all cells,
                          with the basis in real space (realDim components).
                          Then convert to the basis in mapped space (SpaceDim components)
                          by filling in srcVector.
                        */
                        for (DataIterator ditSrc = srcLayout.dataIterator(); ditSrc.ok(); ++ditSrc)
                          {
                            CH_TIME("computing srcVectorFab");
                            const Box& bxBase = srcLayout[ditSrc];
                            const FArrayBox& cellCentersFab = cellCenters[ditSrc];
                            FArrayBox& srcVectorFab = srcVector[ditSrc];
                            FArrayBox& srcVectorRealFab = srcVectorReal[ditSrc];

                            // int blockNum = coordSysPtr->whichBlockOverlap(bxBase);
                            int blockNum = blockNumber[ditSrc];
                            const NewCoordSys* blockCoordSysPtr =
                              coordSysPtr->getCoordSys(blockNum);
                            if (atCellCenters)
                              { // In srcVectorRealFab, store function with basis in physical space (realDim components) at center of each cell.
                                vectorFunPtr->setVectorFunctionFromPhysical(srcVectorRealFab,
                                                                            bxBase,
                                                                            cellCentersFab);
                                // Convert srcVectorRealFab from basis in physical space (realDim components) to basis in mapped space (SpaceDim components) at cell centers.
                                blockCoordSysPtr->vectorTransformRealToMappedCenterFab(srcVectorRealFab);
                              }
                            else
                              { // In srcVectorRealMoreFab, store function with basis in physical space (realDim components) at center of each cell.
                                Box bx1 = grow(bxBase, 1);
                                FArrayBox srcVectorRealMoreFab(bx1, realDim);
                                vectorFunPtr->setVectorFunctionFromPhysical(srcVectorRealMoreFab,
                                                                            bx1,
                                                                            cellCentersFab);
                                // Convert srcVectorRealMoreFab from basis in physical space (realDim components) to basis in mapped space (SpaceDim components) at cell centers.
                                blockCoordSysPtr->vectorTransformRealToMappedCenterFab(srcVectorRealMoreFab);
                                // Set srcVectorFab to cell-averaged function in basis in mapped space (SpaceDim components).
                                FArrayBox srcVectorMoreFab(intvlVector, srcVectorRealMoreFab);
                                fourthOrderAverageCell(srcVectorMoreFab);
                                srcVectorFab.copy(srcVectorMoreFab, bxBase);
                              }
                          }
                    
                        srcVector.copyTo(dstVectorCalc);
                        mblcPtr->copyToVector(srcVector, dstVectorCalc);

                        /*
                          Got results in dstVectorCalc.
                          Now find dstVectorExact and dstVectorDiff = dstVectorCalc - dstVectorExact.
                        */
                        for (DataIterator ditDst = dstLayout.dataIterator(); ditDst.ok(); ++ditDst)
                          {
                            const Box& bx = dstLayout[ditDst];

                            const FArrayBox& calcFab = dstVectorCalc[ditDst];
                            FArrayBox& exactFab = dstVectorExact[ditDst];
                            FArrayBox& exactRealFab = dstVectorRealExact[ditDst];
                            FArrayBox& diffFab = dstVectorDiff[ditDst];

                            // Find cellCentersFab on bxExpanded.
                            Box bxExpanded = grow(bx, extraVect);
                            FArrayBox cellCentersFab(bxExpanded, realDim);
                            FArrayBox xiCentersFab(bxExpanded, SpaceDim);
                            int blockNum = coordSysPtr->whichBlockOverlap(bx);
                            const NewCoordSys* blockCoordSysPtr =
                              coordSysPtr->getCoordSys(blockNum);
                            blockCoordSysPtr->getCenterMappedCoordinates(xiCentersFab, bxExpanded);
                            blockCoordSysPtr->realCoord(cellCentersFab, xiCentersFab, bxExpanded);
                            if (atCellCenters)
                              { // In exactRealFab, store function with basis in physical space (realDim components) at center of each cell.
                                vectorFunPtr->setVectorFunctionFromPhysical(exactRealFab,
                                                                            bxExpanded, // == bx
                                                                            cellCentersFab);
                                // Convert exactRealFab from basis in physical space (realDim components) to basis in mapped space (SpaceDim components) at cell centers.
                                blockCoordSysPtr->vectorTransformRealToMappedCenterFab(exactRealFab);
                              }
                            else
                              { // In exactRealMoreFab, store function with basis in physical space (realDim components) at center of each cell.
                                FArrayBox exactRealMoreFab(bxExpanded, realDim);
                                vectorFunPtr->setVectorFunctionFromPhysical(exactRealMoreFab,
                                                                            bxExpanded, // == grow(bx, 1)
                                                                            cellCentersFab);
                                // Convert exactRealMoreFab from basis in physical space (realDim components) to basis in mapped space (SpaceDim components) at cell centers.
                                blockCoordSysPtr->vectorTransformRealToMappedCenterFab(exactRealMoreFab);
                            
                                // Set exactFab to cell-averaged function in basis in mapped space (SpaceDim components).
                                FArrayBox exactMoreFab(intvlVector, exactRealMoreFab);
                                fourthOrderAverageCell(exactMoreFab);
                                exactFab.copy(exactMoreFab, bx);
                              }
                            diffFab.copy(calcFab);
                            diffFab.minus(exactFab);
                            // dummy statement in order to get around gdb bug
                            int dummy_unused = 0; dummy_unused = 0;
                          }

                        Vector<RealVect> diffVectorBlockMax =
                          blockVecMaxima(dstVectorDiff, coordSysPtr);
                        Vector<RealVect> exactVectorBlockMax =
                          blockVecMaxima(dstVectorExact, coordSysPtr);
                        RealVect diffVectorAllMax = vectorVecMax(diffVectorBlockMax);
                        RealVect exactVectorAllMax = vectorVecMax(exactVectorBlockMax);

                        for (int idir = 0; idir < SpaceDim; idir++)
                          {
                            Real diffRel = diffVectorAllMax[idir] /
                              exactVectorAllMax[idir];
                            diffNorms.push_back(diffRel);
                          }

                      } // end if (doVector)

                    pout() << setw(4) << domainLength;
                    int nNorms = diffNorms.size();
                    for (int iNorm = 0; iNorm < nNorms; iNorm++)
                      {
                        pout() << "  " << scientific << setprecision(4)
                               << diffNorms[iNorm];
                        if (iGrid > 0)
                          {
                            Real ratio = diffNormsCoarser[iNorm] /
                              diffNorms[iNorm];
                            Real rate = log(ratio) / log(2.0);
                            pout () << " " << fixed << setprecision(2) << rate;
                            if (rate < targetConvergeRate)
                              {
                                passTest = false;
                                status += 1;
                              }
                          }
                        else
                          {
                            pout() << "     ";
                          }
                      }
                    pout() << endl;

                    diffNormsCoarser = Vector<Real>(diffNorms);
                          
                    delete coordSysPtr;
                  } // end loop over grid sizes
                if (passTest)
                  {
                    pout() << "Test passed." << endl;
                  }
                else
                  {
                    pout() << "Test failed." << endl;
                  }
              } // end loop over which dimensions are fixed
          } // end mapping type
        // delete mblcPtr;
      } // end centering type
    if (status==0)
      {
        pout() <<  "All tests passed!\n";
      }
    else
      {
        pout() <<  status << " tests failed!\n";
      }
  } // end scoping trick
#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return status;
}
