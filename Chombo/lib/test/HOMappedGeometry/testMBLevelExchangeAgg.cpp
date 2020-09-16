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
#include "MultiBlockLevelExchangeCenter.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "CylinderEquiangularCS.H"
#include "CubedSphere2DCS.H"
#include "CubedSphereShellCS.H"
#include "DoubleCartesianCS.H"
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
 * Test MultiBlockLevelExchange.
 * petermc, 29 May 2012
 */

enum MapCode
  {
    DOUBLECARTESIAN,
    CYLINDEREQUIANGULAR,
    CUBEDSPHERE2D,
    CUBEDSPHERESHELL
  };
const std::string mapStr[] =
  {
    "DoubleCartesianCS",
    "CylinderEquiangularCS",
    "CubedSphere2DCS",
    "CubedSphereShellCS"
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
Vector<Real> blockMaxima(const LevelData<FArrayBox>& a_data,
                         const LayoutData<int>& a_blockNumber,
                         const MultiBlockCoordSys* a_coordSysPtr,
                         int a_nblocks,
                         const IntVect& a_ghostVect)
{
  CH_assert(a_data.nComp() == 1);
  Vector<Real> blockMax(a_nblocks, 0.);
  const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      int blockNum = a_blockNumber[dit];

      const FArrayBox& dataFab = a_data[dit];
      const Box& bxBase = grids[dit];
      // bxWithin:  ghosted box excluding cells outside the whole domain
      Box bxWithin = grow(bxBase, a_ghostVect);
      a_coordSysPtr->keepInDomain(bxWithin, blockNum);
      // bxInterior:  base box excluding cells next to block boundary
      // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

      // Norms all on bxWithin:  ignore external stuff.
      Real dataFabMax = dataFab.norm(bxWithin, 0);
      if ( (dataFabMax != dataFabMax) || // also known as NaN
           (dataFabMax > blockMax[blockNum]) )
        {
          blockMax[blockNum] = dataFabMax;
        }
    }
# ifdef CH_MPI
  for (int iblock = 0; iblock < a_nblocks; iblock++)
    {
      reduceReal(blockMax[iblock], MPI_MAX);
    }
# endif
  return blockMax;
}

// ---------------------------------------------------------
Vector<RealVect> blockVecMaxima(const LevelData<FArrayBox>& a_data,
                                const LayoutData<int>& a_blockNumber,
                                const MultiBlockCoordSys* a_coordSysPtr,
                                int a_nblocks,
                                const IntVect& a_ghostVect)
{
  CH_assert(a_data.nComp() == SpaceDim);
  Vector<RealVect> blockVecMax(a_nblocks, RealVect::Zero);
  const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      int blockNum = a_blockNumber[dit];

      const FArrayBox& dataFab = a_data[dit];
      const Box& bxBase = grids[dit];
      // bxWithin:  ghosted box excluding cells outside the whole domain
      Box bxWithin = grow(bxBase, a_ghostVect);
      a_coordSysPtr->keepInDomain(bxWithin, blockNum);
      // bxInterior:  base box excluding cells next to block boundary
      // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

      // Max-norms all on bxWithin:  ignore external stuff.
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real dataFabMax =
            dataFab.norm(bxWithin, 0, idir);
          if ( (dataFabMax != dataFabMax) || // also known as NaN
               (dataFabMax > blockVecMax[blockNum][idir]) )
            {
              blockVecMax[blockNum][idir] = dataFabMax;
            }
        }
    }
# ifdef CH_MPI
  for (int iblock = 0; iblock < a_nblocks; iblock++)
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
    centerings.push_back(CENTER);
    centerings.push_back(AVERAGE);
    int nCenterings = centerings.size();

    Vector<MapCode> maps;
    if (SpaceDim >= 2)
      { // MultiBlockLevelExchange works only if SpaceDim >= 2.
        if (SpaceDim <= 3) maps.push_back(CUBEDSPHERESHELL);
        if (SpaceDim <= 3) maps.push_back(CUBEDSPHERE2D);
        if (SpaceDim <= 3) maps.push_back(CYLINDEREQUIANGULAR);
        maps.push_back(DOUBLECARTESIAN);
      }
    int nMaps = maps.size();

    // The big nested loop for different cases.
    for (int indCen = 0; indCen < nCenterings; indCen++)
      {
        CenteringCode codeCen = centerings[indCen];
        RefCountedPtr<MultiBlockLevelExchange> mblexPtr;
        // MultiBlockLevelExchange* mblexPtr;
        switch (codeCen)
          {
          case CENTER :
            {
              mblexPtr = RefCountedPtr<MultiBlockLevelExchangeCenter>(new MultiBlockLevelExchangeCenter());
              // mblexPtr = new MultiBlockLevelExchangeCenter();
              break;
            }
          case AVERAGE :
            {
              mblexPtr = RefCountedPtr<MultiBlockLevelExchangeAverage>(new MultiBlockLevelExchangeAverage());
              // mblexPtr = new MultiBlockLevelExchangeAverage();
              break;
            }
          default : 
            {
              MayDay::Error("unknown centering");
            }
          }
        bool atCellCenters = (codeCen == CENTER);

        // number of ghost cells used for exact evaluation
        int numGhostMore = numGhost;
        if (!atCellCenters)
          { // need an additional layer to use for averaging
            numGhostMore++;
          }

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
              { // Fixing dimension 2, and interpolating in dimensions 0:1.
                fixedDimsAll.push_back(Interval(SpaceDim-1, SpaceDim-1));
              }
            // Empty interval: no fixed dimensions.
            fixedDimsAll.push_back(Interval());
            
            for (int ifixed = 0; ifixed < fixedDimsAll.size(); ifixed++)
              {
                const Interval& fixedDims = fixedDimsAll[ifixed];
                int nfixed = fixedDims.size();
                Vector<int> fixedPt(nfixed);
                // Set the components in all fixed dimensions to 3.
                fixedPt.assign(3);

                IntVect interpUnit = IntVect::Unit;
                for (int idir = fixedDims.begin(); idir <= fixedDims.end(); idir++)
                  {
                    interpUnit[idir] = 0;
                  }

                // interpUnit instead of IntVect::Unit,
                // because no ghost cells in fixed dimensions.
                IntVect ghostVect = numGhost * interpUnit;

                int domainLengthMin = 16;
                if (codeMap == CYLINDEREQUIANGULAR)
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
                           << "Testing MultiBlockLevelExchange"
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

                    DisjointBoxLayout grids(allBoxes, allProcs);

                    // petermc, 19 Nov 2012: WAS numGhostMore, not numGhost.
                    MultiBlockLevelGeom geom(coordSysPtr, grids, numGhost,
                                             fixedDims, fixedPt);
                    // mblexPtr->setGetConditionNumber(true);
                    mblexPtr->define(&geom, numGhost, testOrder);

                    ///////////////////////////////////////////
                    // SCALAR FIELD
                    ///////////////////////////////////////////

                    /*
                      Define data holders for scalar.
                    */
                    int ncomp = 1;
                    Interval intvlExact(0, ncomp-1);
                    Interval intvlCalc(ncomp, 2*ncomp-1);
                    Interval intvlDiff(2*ncomp, 3*ncomp-1);
                    LevelData<FArrayBox> allFun;
                    LevelData<FArrayBox> exactFun, calcFun, diffFun;
                    { CH_TIME("allocate allFun");
                      allFun.define(grids, 3*ncomp, ghostVect);

                      aliasLevelData(exactFun, &allFun, intvlExact);
                      aliasLevelData(calcFun, &allFun, intvlCalc);
                      aliasLevelData(diffFun, &allFun, intvlDiff);
                    }
                
                    /*
                      Find cell centers.
                    */
                    // const LevelData<FArrayBox>& cellCenters = geom.cellCenters();
                    IntVect ghostMoreVect = numGhostMore * IntVect::Unit;
                    LevelData<FArrayBox> cellCenters(grids, realDim, ghostMoreVect);
                
                    // Set cellCenters.
                    for (DataIterator dit = exactFun.dataIterator(); dit.ok(); ++dit)
                      {
                        CH_TIME("computing cellCenters");

                        FArrayBox& cellCentersFab = cellCenters[dit];
                        const Box& baseBox = grids[dit];
                        int thisBlock = coordSysPtr->whichBlock(baseBox);
                        const NewCoordSys* thisCoordSysPtr =
                          coordSysPtr->getCoordSys(thisBlock);
                        Box cellCentersBox = grow(baseBox, ghostMoreVect);
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

                    // Set exact value of function.
                    for (DataIterator dit = exactFun.dataIterator(); dit.ok(); ++dit)
                      {
                        CH_TIME("computing exactFun on FAB");
                        FArrayBox& exactFunFab = exactFun[dit];
                        const Box& finalBox = exactFunFab.box();
                        const FArrayBox& cellCentersFab = cellCenters[dit];
                        if (atCellCenters)
                          {
                            funPtr->setFunctionFromPhysical(exactFunFab,
                                                            finalBox,
                                                            cellCentersFab);
                          }
                        else
                          { // Get average of function on each cell.
                            Box expandedBox = grow(finalBox, 1);
                            FArrayBox exactFunMoreFab(expandedBox, ncomp);
                            funPtr->setFunctionFromPhysical(exactFunMoreFab,
                                                            expandedBox,
                                                            cellCentersFab);
                            fourthOrderAverageCell(exactFunMoreFab);
                            exactFunFab.copy(exactFunMoreFab, finalBox);
                          }
                      }

                    // Fill in calcFun with function values at valid cells.
                    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
                      { CH_TIME("copying exactFun to calcFun");
                        const Box& bxBase = grids[dit];
                        FArrayBox& calcFunFab = calcFun[dit];
                        const FArrayBox& exactFunFab = exactFun[dit];
                        calcFunFab.copy(exactFunFab, bxBase);
                      }
                    calcFun.exchange();

                    mblexPtr->interpGhosts(calcFun);

                    /*
                      Got results.  Now take differences,
                      diffFun = calcFun - exactFun.
                    */

                    const LayoutData<int>& blockNumber = geom.block();
                    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
                      { CH_TIME("calculating diffFun");
                        FArrayBox& diffFunFab = diffFun[dit];
                        const FArrayBox& exactFunFab = exactFun[dit];
                        const FArrayBox& calcFunFab = calcFun[dit];

                        int blockNum = blockNumber[dit];
                        const Box& bxBase = grids[dit];
                        // bxWithin:  ghosted box excluding cells outside the whole domain
                        Box bxWithin = grow(bxBase, ghostVect);
                        coordSysPtr->keepInDomain(bxWithin, blockNum);
                        // bxInterior:  base box excluding cells next to block boundary
                        // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

                        diffFunFab.setVal(0.);
                        diffFunFab.copy(exactFunFab, bxWithin);
                        diffFunFab.minus(calcFunFab, bxWithin, 0, 0);
                        // dummy statement in order to get around gdb bug
                        int dummy_unused = 0; dummy_unused = 0;
                      }
                
                    Vector<Real> diffFunBlockMax =
                      blockMaxima(diffFun, blockNumber, coordSysPtr, nblocks, ghostVect);
                    Vector<Real> exactFunBlockMax =
                      blockMaxima(exactFun, blockNumber, coordSysPtr, nblocks, IntVect::Zero);
                    Real diffFunAllMax = vectorMax(diffFunBlockMax);
                    Real exactFunAllMax = vectorMax(exactFunBlockMax);
                    diffNorms.push_back(diffFunAllMax / exactFunAllMax);

                    if (doVector)
                      {
                        ///////////////////////////////////////////
                        // VECTOR FIELD
                        ///////////////////////////////////////////

                        mblexPtr->defineVector();

                        /*
                          Define data holders for vector.
                        */
                        Interval intvlVectorExact(0, SpaceDim-1);
                        Interval intvlVectorCalc(SpaceDim, 2*SpaceDim-1);
                        Interval intvlVectorDiff(2*SpaceDim, 3*SpaceDim-1);
                        LevelData<FArrayBox> allVectorFun;

                        // In case realDim > SpaceDim, exactVectorRealFun takes more
                        // components than exactVectorFun.
                        Interval intvlVectorRealExact(0, realDim-1);
                        LevelData<FArrayBox> exactVectorRealFun;
                        LevelData<FArrayBox> exactVectorFun, calcVectorFun, diffVectorFun;

                        { CH_TIME("allocate allVectorFun");
                      
                          allVectorFun.define(grids, 3*SpaceDim, ghostVect);
                    
                          aliasLevelData(exactVectorFun, &allVectorFun, intvlVectorExact);
                          aliasLevelData(calcVectorFun, &allVectorFun, intvlVectorCalc);
                          aliasLevelData(diffVectorFun, &allVectorFun, intvlVectorDiff);

                          aliasLevelData(exactVectorRealFun, &allVectorFun, intvlVectorRealExact);
                        }

                        /*
                          Fill in exactVectorRealFun with vector function values
                          at all cells,
                          with the basis in real space (realDim components).
                          Then convert to
                          the basis in mapped space (SpaceDim components)
                          by filling in exactVectorFun.
                        */
                        for (DataIterator dit = exactVectorFun.dataIterator(); dit.ok(); ++dit)
                          {
                            CH_TIME("computing exactVectorFun on FAB");
                            FArrayBox& exactVectorRealFunFab = exactVectorRealFun[dit];
                            FArrayBox& exactVectorFunFab = exactVectorFun[dit];
                            const Box& finalBox = exactVectorFunFab.box();
                            const FArrayBox& cellCentersFab = cellCenters[dit];
                            int blockNum = blockNumber[dit];
                            const NewCoordSys* blockCoordSysPtr =
                              coordSysPtr->getCoordSys(blockNum);
                            if (atCellCenters)
                              {
                                // In exactVectorRealFunFab, store function with
                                // basis in physical space (realDim components)
                                // at center of each cell.
                                vectorFunPtr->setVectorFunctionFromPhysical(exactVectorRealFunFab,
                                                                            finalBox,
                                                                            cellCentersFab);
                                // Convert exactVectorRealFunFab
                                // from basis in physical space (realDim components)
                                // to basis in mapped space (SpaceDim components)
                                // at cell centers.
                                blockCoordSysPtr->vectorTransformRealToMappedCenterFab(exactVectorRealFunFab);
                              }
                            else
                              { // In exactVectorRealFunMoreFab, store function
                                // with basis in physical space (realDim components)
                                // at center of each cell.
                                Box expandedBox = grow(finalBox, 1);
                                FArrayBox exactVectorRealFunMoreFab(expandedBox, realDim);
                                vectorFunPtr->setVectorFunctionFromPhysical(exactVectorRealFunMoreFab,
                                                                        
                                                                            expandedBox,
                                                                            cellCentersFab);
                                // Convert exactVectorRealFunMoreFab
                                // from basis in physical space (realDim components)
                                // to basis in mapped space (SpaceDim components)
                                // at cell centers.
                                blockCoordSysPtr->vectorTransformRealToMappedCenterFab(exactVectorRealFunMoreFab);
                            
                                FArrayBox exactVectorFunMoreFab(intvlVectorExact,
                                                                exactVectorRealFunMoreFab);
                                // Convert from cell-centered exactVectorFunMoreFab
                                // to cell-averaged exactVectorFunFab.
                                // Both are in basis
                                // in mapped space (SpaceDim components).
                                fourthOrderAverageCell(exactVectorFunMoreFab);
                                exactVectorFunFab.copy(exactVectorFunMoreFab,
                                                       finalBox);
                              }
                          }
                    
                        /*
                          Fill in calcVectorFun with function values at all cells,
                          on basis in mapped space (SpaceDim components).
                        */
                        for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
                          { CH_TIME("copying exactVectorFun to calcVectorFun");
                            const Box& bxBase = grids[dit];
                            FArrayBox& calcVectorFunFab = calcVectorFun[dit];
                            const FArrayBox& exactVectorFunFab = exactVectorFun[dit];
                            calcVectorFunFab.copy(exactVectorFunFab, bxBase);
                          }
                    
                        // Fill in within-block ghost cells of calcVectorFun.
                        calcVectorFun.exchange();

                        /*
                          Fill in extra-block ghost cells of calcVectorFun,
                          on basis in mapped space (SpaceDim components).
                        */
                        mblexPtr->interpGhostsVector(calcVectorFun);

                        /*
                          Got results.  Now take differences,
                          diffVectorFun = calcVectorFun - exactVectorFun,
                          on basis in mapped space (SpaceDim components).
                        */
                        for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
                          { CH_TIME("calculating diffVectorFun");
                            FArrayBox& diffVectorFunFab = diffVectorFun[dit];
                            const FArrayBox& exactVectorFunFab = exactVectorFun[dit];
                            const FArrayBox& calcVectorFunFab = calcVectorFun[dit];

                            const Box& bxBase = grids[dit];
                            int blockNum = blockNumber[dit];
                            // bxWithin:  ghosted box excluding cells outside the whole domain
                            Box bxWithin = grow(bxBase, ghostVect);
                            coordSysPtr->keepInDomain(bxWithin, blockNum);
                            // bxInterior:  base box excluding cells
                            // next to block boundary
                            // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter
                            diffVectorFunFab.setVal(0.);
                            diffVectorFunFab.copy(exactVectorFunFab, bxWithin);
                            diffVectorFunFab.minus(calcVectorFunFab, bxWithin, 0, 0, SpaceDim);
                          }

                        Vector<RealVect> diffVectorFunBlockMax =
                          blockVecMaxima(diffVectorFun, blockNumber, coordSysPtr, nblocks, ghostVect);
                        Vector<RealVect> exactVectorFunBlockMax =
                          blockVecMaxima(exactVectorFun, blockNumber, coordSysPtr, nblocks, IntVect::Zero);
                        RealVect diffVectorFunAllMax =
                          vectorVecMax(diffVectorFunBlockMax);
                        RealVect exactVectorFunAllMax =
                          vectorVecMax(exactVectorFunBlockMax);
                        for (int idir = 0; idir < SpaceDim; idir++)
                          {
                            Real diffRel = diffVectorFunAllMax[idir] /
                              exactVectorFunAllMax[idir];
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
              } // end loop over which dimensions are fixed
          } // end mapping type
        // delete mblexPtr;
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
