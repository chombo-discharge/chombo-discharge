#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// testMBLevelCopyTo.cpp
// petermc, 28 Sep 2011

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
      pout() << "Beginning testMBLevelCopyTo.  mem= "
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
    //    int ghostFactor = 1;
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
      //      pp.query("ghost_factor", ghostFactor);
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
    DisjointBoxLayout srcLayout(allBoxes, allProcs);

    // Set dstLayout to be srcLayout expanded by numGhost.
    Vector<Box> allExpandedBoxes(allBoxes.size());
    for (int i = 0; i < allBoxes.size(); i++)
      {
        Box bx(allBoxes[i]);
        bx.grow(numGhost);
        int blockNum = coordSysPtr->whichBlockOverlap(bx);
        coordSysPtr->keepInDomain(bx, blockNum);
        allExpandedBoxes[i] = bx;
      }
    BoxLayout dstLayout(allExpandedBoxes, allProcs);

    //    int numGhostMore = (ghostFactor == 1 && !atCellCenters) ?
    //      (numGhost + 1) : (ghostFactor * numGhost);
    //    IntVect ghostVect = numGhost * IntVect::Unit;
    // IntVect ghostMoreVect = numGhostMore * IntVect::Unit;

    pout() << "ghosts " << numGhost
      // << "  ghostFactor " << ghostFactor
           << "  order " << order << endl;
//     MultiBlockLevelData mbdata;
//     mbdata.setVerbose(verbose);
//     mbdata.setUseAverage(useAverage);
//     mbdata.setAverageGhost(averageGhost);
//     mbdata.setOrder(order);
//     mbdata.setGhostFactor(ghostFactor);
//     mbdata.define(md, srcLayout, ghost, radius, degree);
//     if (gaussQuadrature)
//       mbdata.setGaussQuadrature(quadraturePoints);
//     else if (newtonCotes)
//       mbdata.setNewtonCotesQuadrature(quadraturePoints);

    /*
      Define MultiBlockLevelGeom and MultiBlockLevelCopier.
    */
    int extraCentering = (atCellCenters) ? 0 : 1;
    MultiBlockLevelGeom geom(coordSysPtr, srcLayout, extraCentering);
    const LayoutData<int>& blockNumber = geom.block(); // on srcLayout
    MultiBlockLevelCopier* mblcPtr;
    if (atCellCenters)
      { // cell-centered
        mblcPtr = new MultiBlockLevelCopierCenter();
      }
    else
      { // cell-averaged
        mblcPtr = new MultiBlockLevelCopierAverage();
      }
    mblcPtr->setGetConditionNumber(getConditionNumber);
    // mblcPtr->define(&geom, numGhost, order);
    mblcPtr->define(&geom, dstLayout, order);
    const LayoutData< IVSFAB<Real>* >& conditionNumber =
      mblcPtr->conditionNumber();

    if (showStencils)
      {
        const LayoutData< IntVectSet >& extraCells =
          mblcPtr->extraCells();
        const LayoutData< IVSFAB< Vector<IntVect>* >* >& allStencils =
          mblcPtr->allStencils();
        const LayoutData< IVSFAB< Vector<Real>* >* >& allWeights =
          mblcPtr->allWeights();
        for (DataIterator ditDst = dstLayout.dataIterator(); ditDst.ok(); ++ditDst)
          {
            const IntVectSet& extraCellsIVS = extraCells[ditDst];
            const IVSFAB< Vector<IntVect>* >& patchStencils = *allStencils[ditDst];
            const IVSFAB< Vector<Real>* >& patchWeights = *allWeights[ditDst];
            for (IVSIterator ivsit(extraCellsIVS); ivsit.ok(); ++ivsit)
              {
                IntVect thisExtraCell = ivsit();
                const Vector<IntVect>& thisStencil =
                  *patchStencils(thisExtraCell, 0);
                const Vector<Real>& thisWeights =
                  *patchWeights(thisExtraCell, 0);
                int len = thisStencil.size();
                CH_assert(thisWeights.size() == len);
                pout() << "EXTRA CELL " << thisExtraCell
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
    Interval intvlCalc(0, ncomp-1);
    Interval intvlExact(ncomp, 2*ncomp-1);
    Interval intvlDiff(2*ncomp, 3*ncomp-1);

    // LevelData<FArrayBox> allFun(srcLayout, 3*ncomp, ghostVect);
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
      Fill in data holders for scalar.
    */

    // Set src to exact value of function on srcLayout.
    const LevelData<FArrayBox>& cellCenters = geom.cellCenters();
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
    for (DataIterator ditDst = dstLayout.dataIterator(); ditDst.ok(); ++ditDst)
      {
        const Box& bx = dstLayout[ditDst];

        const FArrayBox& calcFab = dstCalc[ditDst];
        FArrayBox& exactFab = dstExact[ditDst];
        FArrayBox& diffFab = dstDiff[ditDst];

        // Find cellCentersFab on bxExpanded.
        Box bxExpanded = grow(bx, extraCentering);
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

    Vector<Real> exactBlockMax = blockMaxima(dstExact, coordSysPtr);
    Real exactAllMax = vectorMax(exactBlockMax);

    Vector<Real> diffBlockMax = blockMaxima(dstDiff, coordSysPtr);
    Real diffAllMax = vectorMax(diffBlockMax);

    Vector<Real> condBlockMax(nblocks, 0.);
    Real condAllMax;
    if (getConditionNumber)
      {
        for (DataIterator ditDst = dstLayout.dataIterator(); ditDst.ok(); ++ditDst)
          {
            const Box& bx = dstLayout[ditDst];
            int blockNum = coordSysPtr->whichBlockOverlap(bx);
            const IVSFAB<Real>& conditionNumberIVSFAB = *conditionNumber[ditDst];
            const IntVectSet& extraCellsIVS = conditionNumberIVSFAB.getIVS();
            Real condFabMax = 0.;
            for (IVSIterator ivsit(extraCellsIVS); ivsit.ok(); ++ivsit)
              {
                IntVect thisExtraCell = ivsit();
                Real thisConditionNumber = conditionNumberIVSFAB(thisExtraCell, 0);
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
               << " diff <= " << diffBlockMax[iblock]
               << " ; exact <= " << exactBlockMax[iblock];
        if (getConditionNumber)
          {
            pout() << " ; cond <= " << condBlockMax[iblock];
          }
        pout() << endl;
      }
    pout() << "ALL length " << domainLength
           << " ghost " << numGhost
           << " diff <= " << diffAllMax
           << " exact <= " << exactAllMax;
    if (getConditionNumber)
      {
        pout() << " cond <= " << condAllMax;
      }
    pout() << endl;

    pout() << "G = " << numGhost;
    pout() << " , P = " << (order-1);
    pout() << " , h = 1/" << domainLength;
    pout() << " , max(diff)/max(exact) = " << (diffAllMax / exactAllMax);
    pout() << endl;

    /*
      File output section
    */


    /*
    if (writeOutput)
      {
        Vector<DisjointBoxLayout> vecLayouts(1);
        vecLayouts[0] = srcLayout;
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
    */

    if (doVector)
      {

    ///////////////////////////////////////////////////////////////////
    // VECTOR FIELD
    ///////////////////////////////////////////////////////////////////

    pout() << "VECTOR FIELD" << endl;
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
        Box bxExpanded = grow(bx, extraCentering);
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

    Vector<RealVect> exactVectorBlockMax =
      blockVecMaxima(dstVectorExact, coordSysPtr);
    RealVect exactVectorAllMax = vectorVecMax(exactVectorBlockMax);

    Vector<RealVect> diffVectorBlockMax =
      blockVecMaxima(dstVectorDiff, coordSysPtr);
    RealVect diffVectorAllMax = vectorVecMax(diffVectorBlockMax);

    for (int iblock = 0; iblock < nblocks; iblock++)
      {
        pout() << "block " << iblock
               << setprecision(4)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific)
               << " diff <=";
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            pout() << " " << diffVectorBlockMax[iblock][idir];
          }
        pout() << " ; exact <=";
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            pout() << " " << exactVectorBlockMax[iblock][idir];
          }
        pout() << endl;
      }
    pout() << "ALL length " << domainLength
           << " ghost " << numGhost
           << " diffVector <=";
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        pout() << " " << diffVectorAllMax[idir];
      }
    pout() << " exactVector <=";
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        pout() << " " << exactVectorAllMax[idir];
      }
    pout() << endl;

    pout() << "G = " << numGhost;
    pout() << " , P = " << (order-1);
    pout() << " , h = 1/" << domainLength;
    pout() << " , max(diffVector)/max(exactVector) =";
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        pout() << " "
               << (diffVectorAllMax[idir] / exactVectorAllMax[idir]);
      }
    pout() << endl;

    // dummy statement in order to get around gdb bug
    int dummy_unused = 0; dummy_unused = 0;
    /*
      File output section
    */

    /*
    if (writeOutput)
      {
        Vector<DisjointBoxLayout> vecLayouts(1);
        vecLayouts[0] = srcLayout;
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
    */
      } // if (doVector)

    delete mblcPtr;
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
