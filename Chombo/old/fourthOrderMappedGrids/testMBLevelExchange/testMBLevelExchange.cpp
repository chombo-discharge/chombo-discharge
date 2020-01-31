#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// testMBLevelExchange.cpp
// petermc, 23 Apr 2010

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
// #include "CylinderSpokes.H"
// #include "CylinderEqualAngles.H"
// #include "CylinderTransition.H"
// #include "Sphere.H"
#include "CylinderEquiangularCS.H"
#include "CubedSphere2DCS.H"
#include "DoubleCartesianCS.H"
#include "CylindricalHarmonic.H"
#include "CylindricalHarmonicGrad.H"
#include "SphericalHarmonic.H"
#include "SphericalHarmonicGrad.H"
#include "Trig2.H"
#include "Trig2Grad.H"
#include "FourthOrderUtil.H"
#include "GaussianAdvectFun.H"
// for parallel debugging
#include "CH_Attach.H"
using std::cerr;
using std::endl;
#include "CH_Timer.H"

#include "memusage.H"

enum MapCode
{
  CYLINDERSPOKES,
  CYLINDEREQUIANGULAR,
  CYLINDERTRANSITION,
  CUBEDSPHERE2D,
  DOUBLECARTESIAN
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
      if (a_vec[i] > vecMax)
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
          if (a_vec[i][idir] > vecMax[idir])
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
                         int a_numGhost)
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
      Box bxWithin = grow(bxBase, a_numGhost);
      a_coordSysPtr->keepInDomain(bxWithin, blockNum);
      // bxInterior:  base box excluding cells next to block boundary
      // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

      // Norms all on bxWithin:  ignore external stuff.
      Real dataFabMax = dataFab.norm(bxWithin, 0);
      if (dataFabMax > blockMax[blockNum])
        blockMax[blockNum] = dataFabMax;
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
                                int a_numGhost)
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
      Box bxWithin = grow(bxBase, a_numGhost);
      a_coordSysPtr->keepInDomain(bxWithin, blockNum);
      // bxInterior:  base box excluding cells next to block boundary
      // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

      // Max-norms all on bxWithin:  ignore external stuff.
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real dataFabMax =
            dataFab.norm(bxWithin, 0, idir);
          if (dataFabMax > blockVecMax[blockNum][idir])
            blockVecMax[blockNum][idir] = dataFabMax;
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


// -------------------------------------------------------------------------
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //scoping trick
  {
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

    print_memory_line("before anything");
    get_memory_usage_from_OS();

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#ifndef NDEBUG
    // for parallel debugging
    // registerDebugger();
    setChomboMPIErrorHandler();
#endif
#endif

    if (argc < 2)
      {
        pout() << "Usage:  " << argv[0] << " <input_file_name> " << endl;
        CH_TIMER_REPORT();
#ifdef CH_MPI
        MPI_Finalize();
#endif
        exit(0);
      }

    {
      pout() << "Beginning testMBLevelExchange.  mem= "
             << get_memory_usage_from_OS() << " MB" << endl;
    }

    /* ==============================================================

    Read from input file.

    ============================================================== */

    FunctionOnSpace* funPtr;
    VectorFunctionOnSpace* vectorFunPtr;
    //    MultiBlock* mbclassPtr;
    int verbose = 1;
    //    int radius = 1; // not used
    // int degree = 2; // not used
    int domainLength = 16;
    int numGhost = 2;
    int order = 2;
    // MUST SET ghostFactor.  Surprised that 1 works.
    int ghostFactor = 1;
    bool atCellCenters = false;
    bool getConditionNumber = false;
    bool writeOutput = false;
    bool getDistances = false;
    bool newtonCotes = false;
    bool gaussQuadrature = true;
    int quadraturePoints = 4;
    MapCode funCode;
    RealVect centerPoint = RealVect::Zero;
    RealVect centralRectSize = RealVect::Unit;
    Real outerRadius = 1.5;
    // for sphere
    int legendreL = 4;
    int legendreM = 1;
    // for cylinder
    int useGaussianFunction = 0;
    int cylinderN = 4;
    int cylinderK = 1;
    // for gaussian
    Real gaussianRadius = 0.15;
    RealVect gaussianCenter = RealVect(D_DECL(1.5, 1.7, 0.5));
    Real gaussianOmega = 1.0;
    RealVect gaussianRotationCenter = RealVect::Zero;
    bool showStencils = false;

    // Maximum dimension of a grid
    int maxGridSize = 32;
    // Minimum dimension of a grid
    int blockFactor = 4;
    // Do vector field too?
    bool doVector = true;
    std::string fileName = "out.hdf5";
    std::string fileVectorName = "outVector.hdf5";
    {
      char* in_file = argv[1];
      ParmParse ppOriginal(argc-2, argv+2, NULL, in_file);
      ParmParse pp("main");

      string funcString;
      int gotFunction = pp.query("function", funcString);
      if (gotFunction != 0)
        {
          if (funcString == "cylinderspokes")
            funCode = CYLINDERSPOKES;
          else if (funcString == "cylinderequiangular")
            funCode = CYLINDEREQUIANGULAR;
          else if (funcString == "cylindertransition")
            funCode = CYLINDERTRANSITION;
          else if (funcString == "cubedsphere2d")
            funCode = CUBEDSPHERE2D;
          else if (funcString == "doublecartesian")
            funCode = DOUBLECARTESIAN;
        }
      pout() << funcString;

      pp.query("verbose", verbose);
      // pp.query("radius", radius);
      // pp.query("degree", degree);
      pp.query("domain_length", domainLength);
      pp.query("ghost", numGhost);
      pp.query("ghost_factor", ghostFactor);
      pp.query("order", order);

      switch (funCode)
        {
        case CUBEDSPHERE2D:
          {
            ParmParse ppSphere("sphere");
            ppSphere.query("l", legendreL);
            ppSphere.query("m", legendreM);
            pout() << " l = " << legendreL << " m = " << legendreM;
            funPtr = new SphericalHarmonic(legendreL, legendreM);
            vectorFunPtr = new SphericalHarmonicGrad(legendreL, legendreM);
            break;
          }
        case CYLINDERSPOKES:
        case CYLINDEREQUIANGULAR:
        case CYLINDERTRANSITION:
          {
            ParmParse ppCylinder("cylinder");
            ppCylinder.query("n", cylinderN);
            ppCylinder.query("k", cylinderK);
            pout() << " n = " << cylinderN << " k = " << cylinderK;
            ppCylinder.query("useGaussianFunction", useGaussianFunction);
            if (useGaussianFunction == 1)
              {
                ppCylinder.query("radius", gaussianRadius);
                ppCylinder.query("omega", gaussianOmega);
                Vector<Real> rotCtra(SpaceDim);
                ppCylinder.queryarr("rotationCenter", rotCtra, 0, SpaceDim);
                D_TERM(gaussianRotationCenter[0]=rotCtra[0];,
                       gaussianRotationCenter[1]=rotCtra[1];,
                       gaussianRotationCenter[2]=rotCtra[2];)
                  std::vector<Real> centera(SpaceDim,0.0);
                ppCylinder.queryarr("center",centera,0,SpaceDim);
                for ( int d=0 ; d<SpaceDim ; ++d ) gaussianCenter[d] = centera[d] ;
                Real evalTime = 0.;
                ppCylinder.query("time", evalTime);
                // x0 is location of lower left corner
                RealVect x0 = RealVect::Zero;
                funPtr = new GaussianAdvectFun(gaussianRadius,
                                               gaussianCenter,
                                               x0,
                                               evalTime,
                                               gaussianRotationCenter,
                                               gaussianOmega);
              }
            else
              {
                funPtr = new CylindricalHarmonic(cylinderN, cylinderK);
              }
            vectorFunPtr = new CylindricalHarmonicGrad(cylinderN, cylinderK);
            break;
          }
        case DOUBLECARTESIAN:
          {
            funPtr = new Trig2();
            vectorFunPtr = new Trig2Grad();
            break;
          }
        }
      pout() << endl;

      {
        int getDistancesInt = 0;
        pp.query("get_distances", getDistancesInt);
        getDistances = (getDistancesInt == 1);
      }

      {
        int getConditionNumberInt = 0;
        pp.query("get_cond", getConditionNumberInt);
        getConditionNumber = (getConditionNumberInt == 1);
      }

      {
        int doVectorInt = 1;
        pp.query("do_vector", doVectorInt);
        doVector = (doVectorInt == 1);
      }

      { // atCellCenters:  whether to use cell-centered values (as
        // opposed to cell-averaged values)
        int atCellCentersInt = 0;
        pp.query("cell_centers", atCellCentersInt);
        atCellCenters = (atCellCentersInt == 1);
      }

      {
        int showStencilsInt = 0;
        pp.query("show_stencils", showStencilsInt);
        showStencils = (showStencilsInt == 1);
      }

      {
        int writeOutputInt = 0;
        pp.query("write_output", writeOutputInt);
        writeOutput = (writeOutputInt == 1);
      }
      if (writeOutput)
        {
          pp.query("file_name", fileName);
          pp.query("file_vector_name", fileVectorName);
        }

      string quadString;
      int gotMethod = pp.query("quadrature", quadString);
      if (gotMethod != 0)
        {
          if (quadString == "gaussian")
            {
              gaussQuadrature = true;
              newtonCotes = false;
            }
          else if (quadString == "newton_cotes")
            {
              newtonCotes = true;
              gaussQuadrature = false;
            }
        }

      pp.query("quadrature_points", quadraturePoints);

      vector<Real> centerPointVect(SpaceDim, 0.);
      if (pp.queryarr("center_point", centerPointVect, 0, SpaceDim))
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            centerPoint[idir] = centerPointVect[idir];
        }

      vector<Real> centralRectSizeVect(SpaceDim, 0.);
      if (pp.queryarr("central_rect_size", centralRectSizeVect, 0, SpaceDim))
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            centralRectSize[idir] = centralRectSizeVect[idir];
        }

      pp.query("outer_radius", outerRadius);

      pp.query("max_grid_size", maxGridSize);
      pp.query("block_factor", blockFactor);
    }

    // one of them, but not both
    CH_assert(newtonCotes || gaussQuadrature);
    CH_assert(! (newtonCotes && gaussQuadrature));

    // Box domain(IntVect::Zero, (domainLength - 1) * IntVect::Unit);
    int realDim = SpaceDim;

    MultiBlockCoordSysFactory* coordSysFactPtr;
    //    SetMappedDomain* smdPtr;
    IntVect levelDomainLo, levelDomainHi;
    switch (funCode)
      {
      case CYLINDERSPOKES:
        {
          // smdPtr = new CylindricalSpokesDomain(domain, RealVect::Zero,
          //                                      bxWidth,
          //                                      outerRadius);
          //           ghostFactor = radius + 2;
          MayDay::Error("CylinderSpokes not implemented yet.");
          break;
        }
      case CYLINDEREQUIANGULAR:
        {
          CylinderEquiangularCSFactory* cylinderCSFactPtr =
            new CylinderEquiangularCSFactory;
          cylinderCSFactPtr->setCenterPoint(centerPoint);
          cylinderCSFactPtr->setCentralRectSize(centralRectSize);
          cylinderCSFactPtr->setOuterRadius(outerRadius);
          coordSysFactPtr = cylinderCSFactPtr;
          levelDomainLo = IntVect(D_DECL(-2*domainLength,
                                         -2*domainLength,
                                         0));
          levelDomainHi = IntVect(D_DECL(3*domainLength-1,
                                         3*domainLength-1,
                                         domainLength-1));
          break;
        }
      case CYLINDERTRANSITION:
        {
          // smdPtr = new CylindricalTransitionDomain(domain, RealVect::Zero,
          //                                          bxWidth, outerRadius);
          // ghostFactor = radius + 1;
          MayDay::Error("CylinderTransition not implemented yet.");
          break;
        }
      case CUBEDSPHERE2D:
        {
          // smdPtr = new SphericalDomain(domain, RealVect::Zero,
          //                              bxWidth, outerRadius);
          // ghostFactor = radius + 2;
          CubedSphere2DCSFactory* cubedSphere2DCSFactPtr =
            new CubedSphere2DCSFactory;
          // cubedSphere2DCSFactPtr->setCenterPoint(centerPoint);
          // cubedSphere2DCSFactPtr->setCentralRectSize(centralRectSize);
          // cubedSphere2DCSFactPtr->setOuterRadius(outerRadius);
          coordSysFactPtr = cubedSphere2DCSFactPtr;
          levelDomainLo = IntVect::Zero;
          levelDomainHi = IntVect(D_DECL(11*domainLength-1,
                                         domainLength-1,
                                         0)); // third dimension unused
          realDim = 3; // but SpaceDim == 2
          break;
        }
      case DOUBLECARTESIAN:
        {
          DoubleCartesianCSFactory* doubleCartesianCSFactPtr =
            new DoubleCartesianCSFactory;
          coordSysFactPtr = doubleCartesianCSFactPtr;
          levelDomainLo = IntVect::Zero;
          levelDomainHi = (3*domainLength - 1) * IntVect::Unit;
          break;
        }
      }
    Box levelDomainBox(levelDomainLo, levelDomainHi);
    ProblemDomain levelDomain(levelDomainBox);
    RealVect dx = (1.0 / Real(domainLength)) * RealVect::Unit;
    MultiBlockCoordSys* coordSysPtr =
      coordSysFactPtr->getCoordSys(levelDomain, dx);

    if (useGaussianFunction == 1)
      {
        GaussianAdvectFun* gaussianFunPtr = (GaussianAdvectFun*) funPtr;
        gaussianFunPtr->setCoordSys(coordSysPtr);
      }

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
    Vector<int> allProcs(allBoxes.size());
    LoadBalance(allProcs, allBoxes);
    DisjointBoxLayout grids(allBoxes, allProcs);

    int numGhostMore = (ghostFactor == 1 && !atCellCenters) ?
      (numGhost + 1) : (ghostFactor * numGhost);
    IntVect ghostVect = numGhost * IntVect::Unit;
    // IntVect ghostMoreVect = numGhostMore * IntVect::Unit;

    pout() << "ghosts " << numGhost
           << "  ghostFactor " << ghostFactor
           << "  order " << order << endl;
//     MultiBlockLevelData mbdata;
//     mbdata.setVerbose(verbose);
//     mbdata.setUseAverage(useAverage);
//     mbdata.setAverageGhost(averageGhost);
//     mbdata.setOrder(order);
//     mbdata.setGhostFactor(ghostFactor);
//     mbdata.define(md, grids, ghost, radius, degree);
//     if (gaussQuadrature)
//       mbdata.setGaussQuadrature(quadraturePoints);
//     else if (newtonCotes)
//       mbdata.setNewtonCotesQuadrature(quadraturePoints);

    /*
      Define MultiBlockLevelGeom and MultiBlockLevelExchange objects.
    */
    MultiBlockLevelGeom geom(coordSysPtr, grids, numGhostMore);
    MultiBlockLevelExchange* mblexPtr;
    if (atCellCenters)
      { // cell-centered
        mblexPtr = new MultiBlockLevelExchangeCenter();
      }
    else
      { // cell-averaged
        mblexPtr = new MultiBlockLevelExchangeAverage();
      }
    mblexPtr->setGetConditionNumber(getConditionNumber);
    mblexPtr->define(&geom, numGhost, order);
    const LayoutData< IVSFAB<Real>* >& conditionNumber =
      mblexPtr->conditionNumber();

    if (showStencils)
      {
        const LayoutData< IntVectSet >& ghostCells =
          mblexPtr->ghostCells();
        const LayoutData< IVSFAB< Vector<IntVect>* >* >& allStencils =
          mblexPtr->allStencils();
        const LayoutData< IVSFAB< Vector<Real>* >* >& allWeights =
          mblexPtr->allWeights();
        for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
          {
            const IntVectSet& ghostCellsIVS = ghostCells[dit];
            const IVSFAB< Vector<IntVect>* >& patchStencils = *allStencils[dit];
            const IVSFAB< Vector<Real>* >& patchWeights = *allWeights[dit];
            for (IVSIterator ivsit(ghostCellsIVS); ivsit.ok(); ++ivsit)
              {
                IntVect thisGhostCell = ivsit();
                const Vector<IntVect>& thisStencil =
                  *patchStencils(thisGhostCell, 0);
                const Vector<Real>& thisWeights =
                  *patchWeights(thisGhostCell, 0);
                int len = thisStencil.size();
                CH_assert(thisWeights.size() == len);
                pout() << "GHOST CELL " << thisGhostCell
                       << "  stencil size " << len << endl;
                Real sumWeights = 0.;
                for (int i = 0; i < len; i++)
                  {
                    //                    pout() << "  " << thisStencil[i]
                    //                           << "  " << thisWeights[i]
                    //                           << endl;
                    pout() << " " << thisStencil[i];
                    sumWeights += thisWeights[i];
                  }
                //                pout() << "  sum discrepancy " << (sumWeights-1.) << endl;
                pout() << endl;
              }
          }
      }

    ///////////////////////////////////////////////////////////////////
    // SCALAR FIELD
    ///////////////////////////////////////////////////////////////////

    /*
      Define data holders for scalar.
    */
    int ncomp = 1;
    LevelData<FArrayBox> allFun(grids, 3*ncomp, ghostVect);
    Interval intvlExact(0, ncomp-1);
    Interval intvlCalc(ncomp, 2*ncomp-1);
    Interval intvlDiff(2*ncomp, 3*ncomp-1);

    LevelData<FArrayBox> exactFun, calcFun, diffFun;
    aliasLevelData(exactFun, &allFun, intvlExact);
    aliasLevelData(calcFun, &allFun, intvlCalc);
    aliasLevelData(diffFun, &allFun, intvlDiff);

    /*
      Fill in data holders for scalar.
    */
    const LevelData<FArrayBox>& cellCenters = geom.cellCenters();
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
      {
        const Box& bxBase = grids[dit];
        FArrayBox& calcFunFab = calcFun[dit];
        const FArrayBox& exactFunFab = exactFun[dit];
        calcFunFab.copy(exactFunFab, bxBase);
      }
    calcFun.exchange();

    mblexPtr->interpGhosts(calcFun);

    /*
      Got results.  Now take differences.
    */

    // diffFun = calcFun - exactFun

    const LayoutData<int>& blockNumber = geom.block();
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& diffFunFab = diffFun[dit];
        const FArrayBox& exactFunFab = exactFun[dit];
        const FArrayBox& calcFunFab = calcFun[dit];

        int blockNum = blockNumber[dit];
        const Box& bxBase = grids[dit];
        // bxWithin:  ghosted box excluding cells outside the whole domain
        Box bxWithin = grow(bxBase, numGhost);
        coordSysPtr->keepInDomain(bxWithin, blockNum);
        // bxInterior:  base box excluding cells next to block boundary
        // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

        diffFunFab.setVal(0.);
        diffFunFab.copy(exactFunFab, bxWithin);
        diffFunFab.minus(calcFunFab, bxWithin, 0, 0);
      }

    Vector<Real> diffFunBlockMax =
      blockMaxima(diffFun, blockNumber, coordSysPtr, nblocks, numGhost);
    Real diffFunAllMax = vectorMax(diffFunBlockMax);
    Vector<Real> exactFunBlockMax =
      blockMaxima(exactFun, blockNumber, coordSysPtr, nblocks, numGhost);
    Real exactFunAllMax = vectorMax(exactFunBlockMax);

    Vector<Real> condBlockMax(nblocks, 0.);
    Real condAllMax;
    if (getConditionNumber)
      {
        for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
          {
            int blockNum = blockNumber[dit];

            const Box& bxBase = grids[dit];
            // bxWithin:  ghosted box excluding cells outside the whole domain
            Box bxWithin = grow(bxBase, numGhost);
            coordSysPtr->keepInDomain(bxWithin, blockNum);

            // bxInterior:  base box excluding cells next to block boundary
            // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

            const IVSFAB<Real>& conditionNumberIVSFAB = *conditionNumber[dit];
            const IntVectSet& ghostCellsIVS = conditionNumberIVSFAB.getIVS();
            Real condFabMax = 0.;
            for (IVSIterator ivsit(ghostCellsIVS); ivsit.ok(); ++ivsit)
              {
                IntVect thisGhostCell = ivsit();
                Real thisConditionNumber = conditionNumberIVSFAB(thisGhostCell, 0);
                if (thisConditionNumber > condFabMax)
                  condFabMax = thisConditionNumber;
              }
            if (condFabMax > condBlockMax[blockNum])
              condBlockMax[blockNum] = condFabMax;
          }
# ifdef CH_MPI
        for (int iblock = 0; iblock < nblocks; iblock++)
          {
            reduceReal(condBlockMax[iblock], MPI_MAX);
          }
# endif
        condAllMax = vectorMax(condBlockMax);
      }

    //    Real funMax = funPtr->setFunctionMax(bxWidth, outerRadius);
    //    pout() << "funMax = " << funMax << endl;
    for (int iblock = 0; iblock < nblocks; iblock++)
      {
        pout() << "block " << iblock
               << setprecision(4)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific)
               << " diff <= " << diffFunBlockMax[iblock]
               << " ; exact <= " << exactFunBlockMax[iblock];
        if (getConditionNumber)
          {
            pout() << " ; cond <= " << condBlockMax[iblock];
          }
        pout() << endl;
      }
    pout() << "ALL length " << domainLength
           << " ghost " << numGhost
           << " diff <= " << diffFunAllMax
           << " exact <= " << exactFunAllMax;
    if (getConditionNumber)
      {
        pout() << " cond <= " << condAllMax;
      }
    pout() << endl;

    pout() << "G = " << numGhost;
    pout() << " , P = " << (order-1);
    pout() << " , h = 1/" << domainLength;
    pout() << " , max(diff)/max(exact) = " << (diffFunAllMax / exactFunAllMax);
    pout() << endl;

    /*
      File output section
    */

    if (writeOutput)
      {
        Vector<DisjointBoxLayout> vecLayouts(1);
        vecLayouts[0] = grids;
        Vector<LevelData<FArrayBox>* > vecData(1);
        vecData[0] = &allFun;
        Vector<string> vecNames(3);
        vecNames[0] = "exact";
        vecNames[1] = "calc";
        vecNames[2] = "diff";
        Real dx0 = dx[0];
        Real dt = 1.;
        Real time = 0.;
        Vector<int> vecRatio(1);
        vecRatio[0] = 1;
        int numLevels = 1;
        WriteAMRHierarchyHDF5(fileName,
                              vecLayouts,
                              vecData,
                              vecNames,
                              levelDomainBox,
                              dx0, dt, time,
                              vecRatio,
                              numLevels);
      }

    if (doVector)
      {

    ///////////////////////////////////////////////////////////////////
    // VECTOR FIELD
    ///////////////////////////////////////////////////////////////////

    pout() << "VECTOR FIELD" << endl;
    mblexPtr->defineVector();

    /*
      Define data holders for vector.
    */
    LevelData<FArrayBox> allVectorFun(grids, 3*SpaceDim, ghostVect);
    Interval intvlVectorExact(0, SpaceDim-1);
    Interval intvlVectorCalc(SpaceDim, 2*SpaceDim-1);
    Interval intvlVectorDiff(2*SpaceDim, 3*SpaceDim-1);

    LevelData<FArrayBox> exactVectorFun, calcVectorFun, diffVectorFun;
    aliasLevelData(exactVectorFun, &allVectorFun, intvlVectorExact);
    aliasLevelData(calcVectorFun, &allVectorFun, intvlVectorCalc);
    aliasLevelData(diffVectorFun, &allVectorFun, intvlVectorDiff);

    // In case realDim > SpaceDim, exactVectorRealFun takes more
    // components than exactVectorFun.
    LevelData<FArrayBox> exactVectorRealFun;
    Interval intvlVectorRealExact(0, realDim-1);
    aliasLevelData(exactVectorRealFun, &allVectorFun, intvlVectorRealExact);

    /*
      Fill in exactVectorRealFun with vector function values at all cells,
      with the basis in real space (realDim components).
      Then convert to the basis in mapped space (SpaceDim components)
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
          { // In exactVectorRealFunFab, store function with basis in physical space (realDim components) at center of each cell.
            vectorFunPtr->setVectorFunctionFromPhysical(exactVectorRealFunFab,
                                                        finalBox,
                                                        cellCentersFab);
            // Convert exactVectorRealFunFab from basis in physical space (realDim components) to basis in mapped space (SpaceDim components) at cell centers.
            blockCoordSysPtr->vectorTransformRealToMappedCenterFab(exactVectorRealFunFab);
          }
        else
          { // In exactVectorRealFunMoreFab, store function with basis in physical space (realDim components) at center of each cell.
            Box expandedBox = grow(finalBox, 1);
            FArrayBox exactVectorRealFunMoreFab(expandedBox, realDim);
            vectorFunPtr->setVectorFunctionFromPhysical(exactVectorRealFunMoreFab,
                                                        expandedBox,
                                                        cellCentersFab);
            // Convert exactVectorRealFunMoreFab from basis in physical space (realDim components) to basis in mapped space (SpaceDim components) at cell centers.
            blockCoordSysPtr->vectorTransformRealToMappedCenterFab(exactVectorRealFunMoreFab);

            FArrayBox exactVectorFunMoreFab(intvlVectorExact,
                                            exactVectorRealFunMoreFab);
            // Convert from cell-centered exactVectorFunMoreFab to cell-averaged exactVectorFunFab.  Both are in basis in mapped space (SpaceDim components).
            fourthOrderAverageCell(exactVectorFunMoreFab);
            exactVectorFunFab.copy(exactVectorFunMoreFab, finalBox);
          }
      }

    /*
      Fill in calcVectorFun with function values at all cells,
      on basis in mapped space (SpaceDim components).
    */
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& calcVectorFunFab = calcVectorFun[dit];
        const FArrayBox& exactVectorFunFab = exactVectorFun[dit];
        // Earlier, I had this copy only on bxBase = grids[dit],
        // but changed because mblexPtr->interpGhostsVector(calcVectorFun)
        // calculates gradients, and to find gradients we need data.
        calcVectorFunFab.copy(exactVectorFunFab);
      }

    // Fill in within-block ghost cells of calcVectorFun.
    calcVectorFun.exchange();

    /*
      Fill in extra-block ghost cells of calcVectorFun,
      on basis in mapped space (SpaceDim components).
     */
    mblexPtr->interpGhostsVector(calcVectorFun);

    /*
      Got results.  Now take differences.
    */

    // diffVectorFun = calcVectorFun - exactVectorFun,
    // on basis in mapped space (SpaceDim cponents).
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& diffVectorFunFab = diffVectorFun[dit];
        const FArrayBox& exactVectorFunFab = exactVectorFun[dit];
        const FArrayBox& calcVectorFunFab = calcVectorFun[dit];

        const Box& bxBase = grids[dit];
        int blockNum = blockNumber[dit];
        // bxWithin:  ghosted box excluding cells outside the whole domain
        Box bxWithin = grow(bxBase, numGhost);
        coordSysPtr->keepInDomain(bxWithin, blockNum);
        // bxInterior:  base box excluding cells next to block boundary
        // Box bxInterior = bxBase & grow(mdDomain, -1); // for diffCenter

        diffVectorFunFab.setVal(0.);
        diffVectorFunFab.copy(exactVectorFunFab, bxWithin);
        diffVectorFunFab.minus(calcVectorFunFab, bxWithin, 0, 0, SpaceDim);
      }

    Vector<RealVect> exactVectorFunBlockMax =
      blockVecMaxima(exactVectorFun, blockNumber, coordSysPtr, nblocks, numGhost);
    RealVect exactVectorFunAllMax = vectorVecMax(exactVectorFunBlockMax);

    Vector<RealVect> diffVectorFunBlockMax =
      blockVecMaxima(diffVectorFun, blockNumber, coordSysPtr, nblocks, numGhost);
    RealVect diffVectorFunAllMax = vectorVecMax(diffVectorFunBlockMax);

    for (int iblock = 0; iblock < nblocks; iblock++)
      {
        pout() << "block " << iblock
               << setprecision(4)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific)
               << " diff <=";
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            pout() << " " << diffVectorFunBlockMax[iblock][idir];
          }
        pout() << " ; exact <=";
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            pout() << " " << exactVectorFunBlockMax[iblock][idir];
          }
        pout() << endl;
      }
    pout() << "ALL length " << domainLength
           << " ghost " << numGhost
           << " diffVector <=";
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        pout() << " " << diffVectorFunAllMax[idir];
      }
    pout() << " exactVector <=";
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        pout() << " " << exactVectorFunAllMax[idir];
      }
    pout() << endl;

    pout() << "G = " << numGhost;
    pout() << " , P = " << (order-1);
    pout() << " , h = 1/" << domainLength;
    pout() << " , max(diffVector)/max(exactVector) =";
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        pout() << " "
               << (diffVectorFunAllMax[idir] / exactVectorFunAllMax[idir]);
      }
    pout() << endl;

    /*
      File output section
    */

    if (writeOutput)
      {
        Vector<DisjointBoxLayout> vecLayouts(1);
        vecLayouts[0] = grids;
        Vector<LevelData<FArrayBox>* > vecData(1);
        vecData[0] = &allVectorFun;
        Vector<string> vecNames;
        D_EXPR(vecNames.push_back("exact0"),
               vecNames.push_back("exact1"),
               vecNames.push_back("exact2"));
        D_EXPR(vecNames.push_back("calc0"),
               vecNames.push_back("calc1"),
               vecNames.push_back("calc2"));
        D_EXPR(vecNames.push_back("diff0"),
               vecNames.push_back("diff1"),
               vecNames.push_back("diff2"));
        Real dx0 = dx[0];
        Real dt = 1.;
        Real time = 0.;
        Vector<int> vecRatio(1);
        vecRatio[0] = 1;
        int numLevels = 1;
        WriteAMRHierarchyHDF5(fileVectorName,
                              vecLayouts,
                              vecData,
                              vecNames,
                              levelDomainBox,
                              dx0, dt, time,
                              vecRatio,
                              numLevels);
      }
      } // if (doVector)

    delete mblexPtr;
    delete funPtr;
    delete vectorFunPtr;
    delete coordSysFactPtr;
    delete coordSysPtr;

    print_memory_line("after everything");
    get_memory_usage_from_OS();

#if defined(TIMER) && defined(CH_MPI)

    // Gather peak memory from procs and write a single line to screen.  that's all.
    Real end_memory = get_memory_usage_from_OS();
    Real avg_memory, min_memory, max_memory;
    gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}
