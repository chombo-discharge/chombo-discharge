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
#include "CubedSphereShellCS.H"
#include "CubedSphereShellPanelCS.H"
#include "ClimateIBCF_F.H"
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
      }
    int nMaps = maps.size();

    // Fixed length in radial dimension.
    int fixedLengthRadial = 32;
    Real magnitude = 2.0 * M_PI;
    Real shift = 0.;
    Real angle = 0.;

    int nVelDims = 2;

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

            switch (codeMap)
              {
              case CUBEDSPHERESHELL :
                {
                  RefCountedPtr<CubedSphereShellCSFactory>
                    cubedSphereShellCSFactPtr =
                    RefCountedPtr<CubedSphereShellCSFactory>(new CubedSphereShellCSFactory);
                  cubedSphereShellCSFactPtr->setHeight(1.e4);
                  cubedSphereShellCSFactPtr->setRadius(6.37122e6);
                  coordSysFactPtr = cubedSphereShellCSFactPtr;

                  basicDomainLo = IntVect::Zero;
                  basicDomainHi = IntVect(D_DECL(11, 1, 1));

                  break;
                }
              default :
                {
                  MayDay::Error("unknown mapping");
                }
              }

            Vector<Interval> fixedDimsAll;
            // Fixing dimension 2, and interpolating in dimensions 0:1.
            fixedDimsAll.push_back(Interval(SpaceDim-1, SpaceDim-1));
            // Empty interval: no fixed dimensions.
            // fixedDimsAll.push_back(Interval());
            
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

#ifdef CH_USE_FLOAT
                // Single precision doesn't get us very far.
                domainLengthMax = 64;
#endif

#ifdef CH_USE_DOUBLE
#if CH_SPACEDIM >= 3
                int domainLengthMax = 64;
                if (nfixed > 0) domainLengthMax = 256;
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

                    pout() << "size";
                    if (doVector)
                      {
                        for (int idir = 0; idir < nVelDims; idir++)
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
                    IntVect domainLengthVect = domainLength * IntVect::Unit;
                    if (SpaceDim >= 3) domainLengthVect[2] = fixedLengthRadial;

                    IntVect levelDomainLo =
                      domainLengthVect * basicDomainLo;
                    IntVect levelDomainHi =
                      domainLengthVect * basicDomainHi - IntVect::Unit;
                    Box levelDomainBox(levelDomainLo, levelDomainHi);
                    ProblemDomain levelDomain(levelDomainBox);
                    RealVect dx = (1.0 / Real(domainLength)) * RealVect::Unit;
                    if (SpaceDim >= 3) dx[2] = 1.0 / Real(fixedLengthRadial);

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

                    IntVect ghostMoreVect = numGhostMore * IntVect::Unit;
                    LevelData<FArrayBox> cellMappedCenters(grids, SpaceDim, ghostMoreVect);
                
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
                        const LayoutData<int>& blockNumber = geom.block();
                        for (DataIterator dit = exactVectorFun.dataIterator(); dit.ok(); ++dit)
                          {
                            CH_TIME("computing exactVectorFun on FAB");
                            FArrayBox& exactVectorRealFunFab = exactVectorRealFun[dit];
                            FArrayBox& exactVectorFunFab = exactVectorFun[dit];
                            const Box& finalBox = exactVectorFunFab.box();
                            int blockNum = blockNumber[dit];
                            const CubedSphereShellPanelCS* blockCoordSysPtr =
                              dynamic_cast<const CubedSphereShellPanelCS*>
                              (coordSysPtr->getCoordSys(blockNum));
                            Box calcBox(finalBox);
                            if (!atCellCenters) calcBox.grow(1);
                            FArrayBox xiFab(calcBox, SpaceDim);

                            FArrayBox rllFab(calcBox, SpaceDim);
                            blockCoordSysPtr->getCellMappedCoordinates(xiFab, calcBox);
                            blockCoordSysPtr->fabTransformEquiangularToLonLat(xiFab, rllFab);
                                                                      
                            FArrayBox vecRLLFab(calcBox, SpaceDim);
                            vecRLLFab.setVal(0.);
                            FORT_CSPVECLONLATSOLIDBODYFLOW(CHF_FRA(vecRLLFab),
                                                           CHF_CONST_FRA(rllFab),
                                                           CHF_CONST_REAL(magnitude),
                                                           CHF_CONST_REAL(angle),
                                                           CHF_CONST_REAL(shift));
                            if (atCellCenters)
                              {
                                // exactVectorRealFunFab: vector in equiangular basis
                                blockCoordSysPtr->fabVectorTransformLatLonToEquiangular(xiFab, vecRLLFab, exactVectorRealFunFab);
                              }
                            else
                              { // In exactVectorRealFunMoreFab, store function
                                // with basis in physical space (realDim components)
                                // at center of each cell.
                                FArrayBox exactVectorRealFunMoreFab(calcBox, realDim);

                                // exactVectorRealFunFab: vector in equiangular basis
                                blockCoordSysPtr->fabVectorTransformLatLonToEquiangular(xiFab, vecRLLFab, exactVectorRealFunMoreFab);

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
                    int nNorms = nVelDims; // diffNorms.size();
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
