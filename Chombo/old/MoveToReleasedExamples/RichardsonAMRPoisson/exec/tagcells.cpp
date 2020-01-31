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
using std::pow;

#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "FineInterp.H"
#include "AMRSolver.H"
#include "LevelOp.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "PoissonBC.H"
#include "PoissonOp.H"
#include "Vector.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "IntVectSet.H"
#include "LoadBalance.H"

#include "localfunctions.H"
#include "fortranfunctions.H"
#include "DebugOut.H"

//typedef LevelData<FArrayBox> LDFB;
//typedef DisjointBoxLayout    DBL;
//const IntVect unitVector(IntVect::Unit);



void scaleDownBoundaryCells(Vector<LDFB*>& hog,
                            const Vector<DBL>& DBLVector,
                            const Vector<Box>& boundingBox,
                            const Vector<Real>& DxVector,
                            const Vector<int>& refinementRatioVector,
                            const int baseLevel,
                            const int numLevels,
                            const int component);

void newTagCells(const Vector<LDFB*>&   rhs,
                 const Vector<LDFB*>&   const_phi,
                 const int numLevels,
                 //const AMRSolver& amrSolver_,  // doesnt work
                 const Vector<DBL>& DBLVector,
                 const Vector<Box>& boundingBoxVector,
                 const Vector<Real>&    DxVector,
                 const Vector<int>& refinementRatioVector,
                 const DomainGhostBC& domghostbc,
                 const int baseLevel,
                 const bool verbose,
                 const int maxIterations,
                 const int numberVcyclesBottom,
                 const Real threshold,
                 Vector<IntVectSet>&   tags)
{
  // Something below wanted a non-const phi, but i am not trying to
  // change it and i really want to pass it in as const.  So for now...
  Vector<LDFB*> phi(numLevels, NULL);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    phi[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);
  }
  copyFirstHOGintoSecond(const_phi, phi, numLevels);

  //cout << " numLevels = " << numLevels
  //   << " sizes: rhs = " << rhs.size()
  //   << " phi = " << phi.size()
  //   << " BB = " << boundingBoxVector.size()
  //   << " Dx = " << DxVector.size()
  //   << " RR = " << refinementRatioVector.size()
  //   << " tags = " << tags.size()
  //   << endl;

  // can this be used for DBLVector instead of passing it in?
  //DBL levleDBL = levdiffabs.getBoxes();


  // create modified_phi and set it equal to RHS
  Vector<LDFB*> modified_phi(numLevels, NULL);

  // modified_phi = rhs
  createHOG(modified_phi, numLevels, DBLVector);
  copyFirstHOGintoSecond(rhs, modified_phi, numLevels);

  // compute modified_phi = phi + k1*RHS
  // can't do this easily -- must go thru each dbl of each level.
  for (int ilev=0; ilev<numLevels; ++ilev)
  {
    const Real k1 = -(DxVector[ilev]*DxVector[ilev])/8.0;
    //cout << " k1 = " << k1 << endl;
    const LDFB& lev_phi = *(phi[ilev]);
    LDFB& lev_modified_phi = *(modified_phi[ilev]);
    for (DataIterator dit = lev_modified_phi.dataIterator(); dit.ok(); ++dit)
    {
      lev_modified_phi[dit()] *= k1;
      lev_modified_phi[dit()] += lev_phi[dit()];
    }
  }

  // make a coarser grid
  Vector<DBL> DBLVector_coarse(numLevels);
  Vector<Box> boundingBoxVector_coarse(numLevels);
  Vector<Real> DxVector_coarse(numLevels);
  Vector<LDFB*> mphi_coarse(numLevels, NULL);

  for (int ilev=0; ilev<numLevels; ilev++)
  {
    coarsen(DBLVector_coarse[ilev], DBLVector[ilev], 2);
    boundingBoxVector_coarse[ilev] = coarsen(boundingBoxVector[ilev], 2);
    DxVector_coarse[ilev] = 2.0*DxVector[ilev];
    mphi_coarse[ilev] = new LDFB(DBLVector_coarse[ilev], 1, unitVector);
  }

  // average soln to coarser levels
  for (int ilev=numLevels-1; ilev>=0; ilev--)
  {
    //cout << DBLVector[ilev] << endl;
    CoarseAverage avgDown(DBLVector[ilev], 1, 2);
    avgDown.averageToCoarse(*mphi_coarse[ilev], *modified_phi[ilev]);
  }

  //outputHDF5(rhs, DBLVector, boundingBoxVector,
  //     refinementRatioVector, numLevels, "rhs.hdf5");
  //outputHDF5(phi, DBLVector, boundingBoxVector,
  //     refinementRatioVector, numLevels, "phi.hdf5");
  //outputHDF5(mphi_coarse, DBLVector_coarse, boundingBoxVector_coarse,
  //     refinementRatioVector, numLevels, "mphi_coarse.hdf5");

  // so mphi_coarse is our 'answer' so far.

  // make new operator and solver on coarse level
  PoissonOp poissonOperator_coarse;
  poissonOperator_coarse.setDomainGhostBC(domghostbc);
  AMRSolver amrSolver_coarse(DBLVector_coarse,
                             boundingBoxVector_coarse,
                             DxVector_coarse,
                             refinementRatioVector,
                             numLevels, baseLevel,
                             &poissonOperator_coarse);

  amrSolver_coarse.setVerbose(verbose);
  amrSolver_coarse.setMaxIter(maxIterations);
  amrSolver_coarse.setNumVCyclesBottom(numberVcyclesBottom);

  // Apply AMR Operator to obtain the Laplacian of phi on coarse grid
  Vector<LDFB*> lap_mphi_coarse(numLevels, NULL);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    lap_mphi_coarse[ilev] = new LDFB(DBLVector_coarse[ilev], 1, unitVector);
    amrSolver_coarse.applyAMROperator(mphi_coarse, *lap_mphi_coarse[ilev], ilev);
  }
  //  outputHDF5(lap_mphi_coarse, DBLVector_coarse, boundingBoxVector_coarse,
  //     refinementRatioVector, numLevels, "coarse_lap.hdf5");

  //  now have the laplacian of mphi_coarse

  // just recreate an AMRSolver since i'm having
  // trouble using the object when i pass it in as a const...
  PoissonOp poissonOperator;
  poissonOperator.setDomainGhostBC(domghostbc);

  AMRSolver amrSolver(DBLVector,  boundingBoxVector,
                      DxVector,   refinementRatioVector,
                      numLevels, baseLevel, &poissonOperator);

  amrSolver.setVerbose(verbose);
  amrSolver.setMaxIter(maxIterations);
  amrSolver.setNumVCyclesBottom(numberVcyclesBottom);


  // Compute the Laplacian of phi
  Vector<LDFB*> phi_lap(numLevels, NULL);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    phi_lap[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);
    amrSolver.applyAMROperator(phi, *phi_lap[ilev], ilev);
  }

  // average Lap of phi to coarser levels
  Vector<LDFB*> lap_phi_coarse_avg(numLevels, NULL);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    lap_phi_coarse_avg[ilev] = new LDFB(DBLVector_coarse[ilev], 1, unitVector);
  }

  for (int ilev=numLevels-1; ilev>=0; ilev--)
  {
    CoarseAverage avgDown(DBLVector[ilev], 1, 2);
    avgDown.averageToCoarse(*lap_phi_coarse_avg[ilev], *phi_lap[ilev]);
  }

  // Now we have two hierachies (HOGs) that can be compared to eachother.
  // Namely, lap_mphi_coarse and lap_phi_coarse_avg

  Vector<LDFB*> diffraw(numLevels, NULL);
  createHOG(diffraw, numLevels, DBLVector_coarse);

  //compareTwoHOGs(lap_mphi_coarse, lap_phi_coarse_avg, numLevels, 0, temp9);

  // Obtain the raw differences (no abs() )
  // Subtract lap_mphi_coarse from lap_phi_coarse_avg and store into diffraw
  subtractTwoHOGs(lap_mphi_coarse, lap_phi_coarse_avg, numLevels, diffraw);


  char filename[80];
  sprintf(filename, "diffraw_level%d.hdf5", numLevels-1);
  outputHDF5(diffraw, DBLVector_coarse, boundingBoxVector_coarse,
             refinementRatioVector, numLevels, filename);

  // input and output is diffraw
  // average down level L+1 onto L for all levels L=0...numLevels-1
  averageDownHOG(diffraw, refinementRatioVector, numLevels);

  char fn4[80];
  sprintf(fn4, "avg_down_diffraw_level%d.hdf5", numLevels-1);
  outputHDF5(diffraw, DBLVector_coarse, boundingBoxVector_coarse,
             refinementRatioVector, numLevels, fn4);

  const int component = 0;

  // input and output is diffraw
  //scaleDownBoundaryCells(diffraw, DBLVector_coarse,
  //                     boundingBoxVector_coarse, DxVector,
  //                     refinementRatioVector,
  //                     baseLevel, numLevels, component);


  // now i need to go thru each cell and tag based on threshold

  int numberCoarseFineBoundaryFaces=0;
  int numberProblemBoundaryFaces=0;
  int numberCells;
  int CFhold=0;
  //const int component = 0;

  Vector<LDFB*> colorLDFB(numLevels, NULL);
  for (int i = 0; i < numLevels; i++)
  {
    colorLDFB[i] = new LDFB(DBLVector_coarse[i], 1, unitVector);
  }

  // loop over levels
  for (int ilev=baseLevel; ilev<numLevels; ilev++)
  {

    printf(" tagging: ilev=%d  baseLevel=%d  numLevels=%d\n",
           ilev, baseLevel, numLevels);

    LDFB& levdiffraw = *diffraw[ilev];
    const DBL currentDBL = levdiffraw.disjointBoxLayout();

    IntVectSet coarseFineSkinSet;
    IntVectSet fineCoarseSkinSet;
    IntVectSet problemBoundarySkinSet;
    numberCoarseFineBoundaryFaces=0;
    numberProblemBoundaryFaces=0;

    // obtain the number of cells in current DBL
    numberCells = cellCount(currentDBL);

    // find the number of faces on boundary and construct a
    // set of all the cells on the problem boundary for this level.
    problemBoundaryCount(currentDBL,
                         boundingBoxVector_coarse[ilev],
                         problemBoundarySkinSet,
                         numberProblemBoundaryFaces);

    // Only compute coarse-fine boundary when there is a finer level
    if (ilev < numLevels-1)
    {
      //cout << " looking for coarse-fine boundaries " << endl;
      LDFB& nextlevdiffraw = *diffraw[ilev+1];
      const DBL nextDBL = nextlevdiffraw.disjointBoxLayout();
      DBL nextDBL_coarse;
      coarsen(nextDBL_coarse, nextDBL, 2);
      coarseFineBoundaryCount(nextDBL_coarse,
                              boundingBoxVector_coarse[ilev],
                              coarseFineSkinSet,
                              numberCoarseFineBoundaryFaces);
    }

    // Only compute fine-coarse boundary when there is a coarser level
    if (ilev > baseLevel)
    {
      //cout << " looking for fine-coarse boundaries " << endl;
      fineCoarseBoundaryCount(currentDBL, boundingBoxVector_coarse[ilev],
                              fineCoarseSkinSet);
    }

    printf("  level%d  numberCoarseFineBoundaryFaces=%d  numberProblemBoundaryFaces=%d  numberCells=%d\n",
           ilev, numberCoarseFineBoundaryFaces, numberProblemBoundaryFaces, numberCells);

    printf("   coarseFine.size=%d  fineCoarse.size=%d  problemBoundary.size=%d\n",
           coarseFineSkinSet.numPts(), fineCoarseSkinSet.numPts(),
           problemBoundarySkinSet.numPts());

    const Real KcoarseFine = (Real)numberCoarseFineBoundaryFaces/(Real)numberCells;
    const Real ND = refinementRatioVector[ilev]*pow(2.0,SpaceDim-2);
    const Real KfineCoarse = ND*CFhold/(Real)numberCells;
    const Real KproblemBoundary = (Real)numberProblemBoundaryFaces/(Real)numberCells
      * DxVector[ilev];
    CFhold = numberCoarseFineBoundaryFaces;

    printf("   KcoarseFine=%12.3e  KfineCoarse=%12.3e  KproblemBoundary=%12.3e\n",
           KcoarseFine, KfineCoarse, KproblemBoundary);


    fillColorFAB(colorLDFB[ilev], coarseFineSkinSet,
                 fineCoarseSkinSet, problemBoundarySkinSet);

    char fn2[80];
    sprintf(fn2, "color%dof%d.hdf5", ilev, numLevels);
    //writeLDFBname(colorLDFB, fn);
    outputHDF5(colorLDFB, DBLVector_coarse, boundingBoxVector_coarse,
               refinementRatioVector, ilev+1, fn2);


    IntVectSet local_tags;

    // loop over boxes in level
    //int boxCount=0;
    for (DataIterator dit=levdiffraw.dataIterator(); dit.ok(); ++dit)
    {

      //printf("   box#%d\n", boxCount);
      //boxCount++;

      const Box& box = levdiffraw.box(dit());
      //FArrayBox grad_fab (b, SpaceDim);
      FArrayBox& fab = levdiffraw[dit()];

      // loop over each cell in box
      for (BoxIterator bit = BoxIterator(box); bit.ok(); ++bit)
      {
        //printf(" %f \n", fab(bit(),component));

        Real Ksmudge = 1.0;

        // if this cell is contained within two skin sets,
        // then use the smallest K?
        if (coarseFineSkinSet.contains(bit()))
        {
          Ksmudge = KcoarseFine;
        }

        if (fineCoarseSkinSet.contains(bit()))
        {
          if (KfineCoarse < Ksmudge) Ksmudge = KfineCoarse;
        }

        if (problemBoundarySkinSet.contains(bit()))
        {
          if (KproblemBoundary < Ksmudge) Ksmudge = KproblemBoundary;
        }

        if ( fabs(Ksmudge*fab(bit(),component)) >= threshold)
        {
          //cout << bit() << endl;
          local_tags |= bit();
        }

        // WACTH OUT!  i'm changing values of diffraw here.  after
        // this loop it will be a scaled down hierarchy and will be
        // used for plotting purposes only
        fab(bit(),component) *= Ksmudge;

        //cout << " num local tags = " << local_tags.numPts() << endl;
      }

    } // end loop on boxes for this level

    //printf("   boxCount=%d\n", boxCount);

    char fn5[80];
    sprintf(fn5, "est_truncation_error_level%d.hdf5", numLevels-1);
    outputHDF5(diffraw, DBLVector_coarse, boundingBoxVector_coarse,
               refinementRatioVector, numLevels, fn5);


    Vector<IntVectSet> all_tags;
    const int dest_proc = uniqueProc (SerialTask::compute);

    //cout << "  gathering.  dest_procs="<< dest_proc << endl;
    gather(all_tags, local_tags, dest_proc);

    if (procID() == dest_proc )
    {
      for (int i = 0; i < all_tags.size(); ++i)
      {
        tags[ilev] |= all_tags[i];
      }
    }

    tags[ilev].refine(2);

    char filename3[80];
    sprintf(filename3, "tags_atLevel%dof%d.hdf5", ilev, numLevels);
    makeHDF5fromIVS(tags[ilev], boundingBoxVector[ilev], filename3);

    cout << "  ilev=" << ilev
         << " numtags = " << tags[ilev].numPts()
         << " num local tags = " << local_tags.numPts()
         << endl;
  } // end loop over levels
}




// this is ridiculous
void makeHDF5fromIVS(const IntVectSet& tags,
                            const Box& boundingBox,
                            const string& name)
{

  if (tags.numPts() <= 0) return;

  const int numLevels = 1;
  Vector<Box> vb(tags.numPts()); // at most?
  vb = tags.boxes();

  Vector<int> pids(vb.size());

  for (int i=0; i<vb.size(); i++)
  {
    pids[i] = 0;
  }

  Vector<DBL> vdbl(numLevels);
  vdbl[0].define(vb,pids);

  Vector<LDFB*> aldf(numLevels, NULL);
  aldf[0] = new LDFB(vdbl[0], 1, unitVector);
  Vector<int> refinementRatioVector(numLevels);
  refinementRatioVector[0] = 2;
  Vector<Box> boundingBoxVector(1);
  boundingBoxVector[0] = boundingBox;
  outputHDF5(aldf, vdbl, boundingBoxVector,
             refinementRatioVector, numLevels, name);
}



void coarseFineBoundaryCount(const DisjointBoxLayout& dbl,
                             const Box& boundingBoxLevel0,
                             IntVectSet& coarseFineSkinSet,
                             int& numberCoarseFineBoundaryFaces)
{

  const bool VV = false;  // verbose debug flag

  const int SX = boundingBoxLevel0.smallEnd()[0];
  const int SY = boundingBoxLevel0.smallEnd()[1];
  const int BX = boundingBoxLevel0.bigEnd()[0];
  const int BY = boundingBoxLevel0.bigEnd()[1];
#if (CH_SPACEDIM == 3)
  const int SZ = boundingBoxLevel0.smallEnd()[2];
  const int BZ = boundingBoxLevel0.bigEnd()[2];
#endif

  numberCoarseFineBoundaryFaces=0;
  int faceCount=0;
  int totalExtra=0;

  // loop over all boxes
  for (LayoutIterator it = dbl.layoutIterator(); it.ok(); ++it)
  {
    //cout << " di()=" << it().intCode() << endl;

    const Box currentBox = dbl[it()];
    //cout << " current box=" << currentBox << endl;

    const int sx = currentBox.smallEnd()[0];
    const int sy = currentBox.smallEnd()[1];
    const int bx = currentBox.bigEnd()[0];
    const int by = currentBox.bigEnd()[1];
#if (CH_SPACEDIM == 3)
    const int sz = currentBox.smallEnd()[2];
    const int bz = currentBox.bigEnd()[2];
#endif

    int a1;
    // Are we within the problem boundary?
    if (sx > SX)
    {
      // yea, so add left side

      // first make a single layer slab on the side (outside the box)
      IntVect siv(D_DECL(sx-1, sy, sz));
      IntVect biv(D_DECL(sx-1, by, bz));
      Box b(siv, biv);

      // count the number of cells -- and since it's
      // only one layer deep, that's also the number of
      // outward facing faces.
      faceCount += b.numPts();

      // Then add this box to the IntVectSet collection.
      // Per IntVectSet rules, identical cells will not be
      // added twice which is what we want here.
      // However, for now, let's keep track of how many
      // "extras" there are.
      a1 = coarseFineSkinSet.numPts() + b.numPts();
      coarseFineSkinSet |= b;
      totalExtra += a1 - coarseFineSkinSet.numPts();
      if (VV) cout << " adding left side = " << b
                  << "  n=" << b.numPts()
                  << "  extra=" << totalExtra
                  << endl;
    }

    if (bx < BX)
    {
      // add right side
      IntVect siv(D_DECL(bx+1, sy, sz));
      IntVect biv(D_DECL(bx+1, by, bz));
      Box b(siv, biv);
      faceCount += b.numPts();

      a1 = coarseFineSkinSet.numPts() + b.numPts();
      coarseFineSkinSet |= b;
      totalExtra += a1 - coarseFineSkinSet.numPts();
      if (VV) cout << " adding right side = " << b
                   << "  n=" << b.numPts()
                   << "  extra=" << totalExtra
                   << endl;
    }

    if (sy > SY)
    {
      // add bottom side
      IntVect siv(D_DECL(sx, sy-1, sz));
      IntVect biv(D_DECL(bx, sy-1, bz));
      Box b(siv, biv);
      faceCount += b.numPts();
      a1 = coarseFineSkinSet.numPts() + b.numPts();
      coarseFineSkinSet |= b;
      totalExtra += a1 - coarseFineSkinSet.numPts();
      if (VV) cout << " adding bottom side = " << b
                  << "  n=" << b.numPts()
                  << "  extra=" << totalExtra
                  << endl;
    }

    if (by < BY)
    {
      // add top side
      IntVect siv(D_DECL(sx, by+1, sz));
      IntVect biv(D_DECL(bx, by+1, bz));
      Box b(siv, biv);
      faceCount += b.numPts();
      a1 = coarseFineSkinSet.numPts() + b.numPts();
      coarseFineSkinSet |= b;
      totalExtra += a1 - coarseFineSkinSet.numPts();
      if (VV) cout << " adding top side = " << b
                  << "  n=" << b.numPts()
                  << "  extra=" << totalExtra
                  << endl;
    }

#if (CH_SPACEDIM == 3)
    if (sz > SZ)
    {
      // add front side
      IntVect siv(D_DECL(sx, sy, sz-1));
      IntVect biv(D_DECL(bx, by, sz-1));
      Box b(siv, biv);
      faceCount += b.numPts();
      a1 = coarseFineSkinSet.numPts() + b.numPts();
      coarseFineSkinSet |= b;
      totalExtra += a1 - coarseFineSkinSet.numPts();
      if (VV) cout << " adding front side = " << b
                  << "  n=" << b.numPts()
                  << "  extra=" << totalExtra
                  << endl;
    }

    if (bz < BZ)
    {
      // add back side
      IntVect siv(D_DECL(sx, sy, bz+1));
      IntVect biv(D_DECL(bx, by, bz+1));
      Box b(siv, biv);
      faceCount += b.numPts();
      a1 = coarseFineSkinSet.numPts() + b.numPts();
      coarseFineSkinSet |= b;
      totalExtra += a1 - coarseFineSkinSet.numPts();
      if (VV) cout << " adding back side = " << b
                  << "  n=" << b.numPts()
                  << "  extra=" << totalExtra
                  << endl;
    }
#endif

  }

  const int totalSkinSetCount = coarseFineSkinSet.numPts();

  // We now have an IntVectSet of all the cells making up
  //  a one-level skin covering all boxes within the problem
  //  (minus any skin going outside of the problem domain)
  //  Some of this skin includes cells of already defined boxes.
  // So now we can go thru each box in the DBL,
  //   and subtract those
  //   cells from our coarseFineSkinSet to obtain only the cells
  //   on a coarse-fine boundary.
  // loop over all boxes
  for (LayoutIterator it = dbl.layoutIterator(); it.ok(); ++it)
  {
    Box b = dbl[it()];
    coarseFineSkinSet -= b;
    //cout << "  subtracting = " << b << endl;
  }

  const int numberSharedFaces =
    totalSkinSetCount - coarseFineSkinSet.numPts();

  if (VV) cout << " totalSkinSetCount " << totalSkinSetCount
              << " numberSharedFaces " << numberSharedFaces
              << " totalExtra= " << totalExtra
              <<endl;

  numberCoarseFineBoundaryFaces = faceCount - numberSharedFaces;

  // And this is a little hack to account for an oversight in
  // my original algorithm.  i'll go back and fix this all later.
  numberCoarseFineBoundaryFaces -= totalExtra;

  // And if you want a set of all the cells comprising the
  //  coarse-fine boundary -- you have it
  //for (IVSIterator it(coarseFineSkinSet); it.ok(); ++it)
  //{
  //cout << " coarseFineSkinSet = " << it()  << endl;
  //}
  //cout << " size = " << coarseFineSkinSet.numPts() << endl;

  //cout << " numberCoarseFineBoundaryFaces " << numberCoarseFineBoundaryFaces << endl;
  //cout << " numberProblemBoundaryFaces    " << numberProblemBoundaryFaces << endl;
}





void problemBoundaryCount(const DisjointBoxLayout& dbl,
                          const Box& boundingBoxLevel0,
                          IntVectSet& problemBoundarySkinSet,
                          int& numberProblemBoundaryFaces)
{

  // Construct an IVS of all cells that touch the problem boundary.
  IntVectSet boundingIVS(boundingBoxLevel0);
  IntVectSet peeledBoundingIVS(boundingBoxLevel0);
  peeledBoundingIVS.grow(-1);
  IntVectSet cellsOnProblemBoundary = boundingIVS - peeledBoundingIVS;

  // loop over all boxes
  for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit)
  {
    const Box currentBox = dbl[lit()];
    problemBoundarySkinSet |= (currentBox & cellsOnProblemBoundary);
  }

  const int SX = boundingBoxLevel0.smallEnd()[0];
  const int SY = boundingBoxLevel0.smallEnd()[1];
  const int BX = boundingBoxLevel0.bigEnd()[0];
  const int BY = boundingBoxLevel0.bigEnd()[1];
#if (CH_SPACEDIM == 3)
  const int SZ = boundingBoxLevel0.smallEnd()[2];
  const int BZ = boundingBoxLevel0.bigEnd()[2];
#endif

  numberProblemBoundaryFaces=0;

  // loop over all boxes
  for (LayoutIterator it = dbl.layoutIterator(); it.ok(); ++it)
  {
    const Box currentBox = dbl[it()];

    const int sx = currentBox.smallEnd()[0];
    const int sy = currentBox.smallEnd()[1];
    const int bx = currentBox.bigEnd()[0];
    const int by = currentBox.bigEnd()[1];
#if (CH_SPACEDIM == 3)
    const int sz = currentBox.smallEnd()[2];
    const int bz = currentBox.bigEnd()[2];
#else
    const int sz = 0;
    const int bz = 0;
#endif

    if (sx <= SX) numberProblemBoundaryFaces += (by-sy+1)*(bz-sz+1);
    if (bx >= BX) numberProblemBoundaryFaces += (by-sy+1)*(bz-sz+1);
    if (sy <= SY) numberProblemBoundaryFaces += (bx-sx+1)*(bz-sz+1);
    if (by >= BY) numberProblemBoundaryFaces += (bx-sx+1)*(bz-sz+1);
#if (CH_SPACEDIM == 3)
    if (sz <= SZ) numberProblemBoundaryFaces += (bx-sx+1)*(by-sy+1);
    if (bz >= BZ) numberProblemBoundaryFaces += (bx-sx+1)*(by-sy+1);
#endif

  }
}


void fineCoarseBoundaryCount(const DisjointBoxLayout& dbl,
                             const Box& boundingBoxLevel0,
                             IntVectSet& fineCoarseSkinSet)
{

  // Construct an IVS of all cells on the fine-coarse boundary...
  //  These are the cells that border a higher level box, not
  //  to be confused with cells that border a lower level box.
  IntVectSet IVS_of_all_boxes_on_level;

  // loop over boxes in level
  for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit)
  {
    const Box currentBox = dbl[lit()];
    IVS_of_all_boxes_on_level |= currentBox;
  }

  const IntVectSet copy_of_IVS_of_all_boxes_on_level = IVS_of_all_boxes_on_level;
  IVS_of_all_boxes_on_level.nestingRegion(1, boundingBoxLevel0);
  fineCoarseSkinSet = copy_of_IVS_of_all_boxes_on_level - IVS_of_all_boxes_on_level;
}


int cellCount(const DisjointBoxLayout& dbl)
{
  int numberCells=0;
  // loop over all boxes
  for (LayoutIterator it = dbl.layoutIterator(); it.ok(); ++it)
  {
    const Box currentBox = dbl[it()];
    numberCells += currentBox.numPts();
  }
  return numberCells;
}


//                                           ....
//                                          W$$$$$u
//                                          $$$$F**+           .oW$$$eu
//                                          ..ueeeWeeo..      e$$$$$$$$$
//                                      .eW$$$$$$$$$$$$$$$b- d$$$$$$$$$$W
//                          ,,,,,,,uee$$$$$$$$$$$$$$$$$$$$$ H$$$$$$$$$$$~
//                       :eoC$$$$$$$$$$$C""?$$$$$$$$$$$$$$$ T$$$$$$$$$$"
//                        $$$*$$$$$$$$$$$$$e "$$$$$$$$$$$$$$i$$$$$$$$F"
//                        ?f"!?$$$$$$$$$$$$$$ud$$$$$$$$$$$$$$$$$$$$*Co
//                        $   o$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                !!!!m.*eeeW$$$$$$$$$$$f?$$$$$$$$$$$$$$$$$$$$$$$$$$$$$U
//                !!!!!! !$$$$$$$$$$$$$$  T$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                 *!!*.o$$$$$$$$$$$$$$$e,d$$$$$$$$$$$$$$$$$$$$$$$$$$$$$:
//                "eee$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
//               b ?$$$$$$$$$$$$$$**$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
//               Tb "$$$$$$$$$$$$$$*uL"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
//                $$o."?$$$$$$$$F" u$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                 $$$$en ```    .e$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
//                  $$$B*  =*"?.e$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$F
//                   $$$W"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
//                    "$$$o#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
//                   R: ?$$$W$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" :!i.
//                    !!n.?$???""``.......,``````"""""""""""``   ...+!!!
//                     !* ,+::!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*`
//                     "!?!!!!!!!!!!!!!!!!!!~ !!!!!!!!!!!!!!!!!!!~`
//                     +!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!?!`
//                   .!!!!!!!!!!!!!!!!!!!!!' !!!!!!!!!!!!!!!, !!!!
//                  :!!!!!!!!!!!!!!!!!!!!!!' !!!!!!!!!!!!!!!!! `!!:
//               .+!!!!!!!!!!!!!!!!!!!!!~~!! !!!!!!!!!!!!!!!!!! !!!.
//              :!!!!!!!!!!!!!!!!!!!!!!!!!.`:!!!!!!!!!!!!!!!!!:: `!!+
//              "~!!!!!!!!!!!!!!!!!!!!!!!!!!.~!!!!!!!!!!!!!!!!!!!!.`!!:
//                  ~~!!!!!!!!!!!!!!!!!!!!!!! ;!!!!~` ..eeeeeeo.`+!.!!!!.
//                :..    `+~!!!!!!!!!!!!!!!!! :!;`.e$$$$$$$$$$$$$u .
//                $$$$$$beeeu..  `````~+~~~~~" ` !$$$$$$$$$$$$$$$$ $b
//                $$$$$$$$$$$$$$$$$$$$$UU$U$$$$$ ~$$$$$$$$$$$$$$$$ $$o
//               !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$. $$$$$$$$$$$$$$$~ $$$u
//               !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! $$$$$$$$$$$$$$$ 8$$$$.
//               !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$X $$$$$$$$$$$$$$`u$$$$$W
//               !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! $$$$$$$$$$$$$".$$$$$$$:
//                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  $$$$$$$$$$$$F.$$$$$$$$$
//                ?$$$$$$$$$$$$$$$$$$$$$$$$$$$$f $$$$$$$$$$$$' $$$$$$$$$$.
//                 $$$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$  $$$$$$$$$$!
//                 "$$$$$$$$$$$$$$$$$$$$$$$$$$$ ?$$$$$$$$$$$$  $$$$$$$$$$!
//                  "$$$$$$$$$$$$$$$$$$$$$$$$Fib ?$$$$$$$$$$$b ?$$$$$$$$$
//                    "$$$$$$$$$$$$$$$$$$$$"o$$$b."$$$$$$$$$$$  $$$$$$$$'
//                   e. ?$$$$$$$$$$$$$$$$$ d$$$$$$o."?$$$$$$$$H $$$$$$$'
//                  $$$W.`?$$$$$$$$$$$$$$$ $$$$$$$$$e. "??$$$f .$$$$$$'
//                 d$$$$$$o "?$$$$$$$$$$$$ $$$$$$$$$$$$$eeeeee$$$$$$$"
//                 $$$$$$$$$bu "?$$$$$$$$$ 3$$$$$$$$$$$$$$$$$$$$*$$"
//                d$$$$$$$$$$$$$e. "?$$$$$:`$$$$$$$$$$$$$$$$$$$$8
//        e$$e.   $$$$$$$$$$$$$$$$$$+  "??f "$$$$$$$$$$$$$$$$$$$$c
//       $$$$$$$o $$$$$$$$$$$$$$$F"          `$$$$$$$$$$$$$$$$$$$$b.
//      M$$$$$$$$U$$$$$$$$$$$$$F"              ?$$$$$$$$$$$$$$$$$$$$$u
//      ?$$$$$$$$$$$$$$$$$$$$F                   "?$$$$$$$$$$$$$$$$$$$$u
//       "$$$$$$$$$$$$$$$$$$"                       ?$$$$$$$$$$$$$$$$$$$$o
//         "?$$$$$$$$$$$$$F                            "?$$$$$$$$$$$$$$$$$$
//            "??$$$$$$$F                                 ""?3$$$$$$$$$$$$F
//                                                      .e$$$$$$$$$$$$$$$$'
//                                                     u$$$$$$$$$$$$$$$$$
//                                                    `$$$$$$$$$$$$$$$$"
//                                                     "$$$$$$$$$$$$F"
//                                                       ""?????""

void fillColorFAB(LDFB* colorLDFB,
                  const IntVectSet& coarseFineSkinSet,
                  const IntVectSet& fineCoarseSkinSet,
                  const IntVectSet& problemBoundarySkinSet)
{

  LDFB& levcolor = *colorLDFB;

  for (DataIterator dit=levcolor.dataIterator(); dit.ok(); ++dit)
  {

    const Box& box = levcolor.box(dit());
    FArrayBox& cfab = levcolor[dit()];

    // loop over each cell in box
    const int component=0;
    for (BoxIterator bit = BoxIterator(box); bit.ok(); ++bit)
    {
      //cout << bit() ;
      cfab(bit(),component) = 0.0;

      if (coarseFineSkinSet.contains(bit()))
      {
        cfab(bit(),component) += 1.0;
      }
      if (fineCoarseSkinSet.contains(bit()))
      {
        cfab(bit(),component) += 3.0;
      }
      if (problemBoundarySkinSet.contains(bit()))
      {
        cfab(bit(),component) += 5.0;
      }
      //printf(" %f \n", cfab(bit(),component));
    }
  }

  // old stuff:
    //this will grab the value from the cell at (3,3)
    //FArrayBox& thisFab = LDFB[dit()];
    //IntVect thisLoc = 3*IntVect::Unit;
    //Real thisVal = thisFab(thisLoc);

//     // now get boxes which make up color sets
//     Vector<Box> CFBoxes = coarseFineSkinSet.boxes();
//     Vector<Box> FCBoxes = fineCoarseSkinSet.boxes();
//     Vector<Box> PBBoxes = problemBoundarySkinSet.boxes();

//     for (adit.reset(); adit.ok(); ++adit)
//     {

//       const Box& box = currentDBL[adit()];

//       // Set all values of boxes labelled as coarse-fine to be 1.0
//       for (int i=0; i<CFBoxes.size(); ++i)
//       {
//         if (CFBoxes[i].intersects(box) )
//         {
//           Box intersectBox(box);
//           intersectBox &= CFBoxes[i];
//           *(colorLDFB[ilev])[adit()].setVal(1.0, intersectBox, 0);
//         }
//       }

//       // Set all values of boxes labelled as fine-coarse to be 2.0
//       for (int i=0; i<FCBoxes.size(); ++i)
//       {
//         if (FCBoxes[i].intersects(box) )
//         {
//           Box intersectBox(box);
//           intersectBox &= FCBoxes[i];
//           colorLDFB[adit()].setVal(2.0, intersectBox, 0);
//         }
//       }

//       // Set all values of boxes labelled as problem-boundary to be 3.0
//       for (int i=0; i<PBBoxes.size(); ++i)
//       {
//         if (PBBoxes[i].intersects(box) )
//         {
//           Box intersectBox(box);
//           intersectBox &= PBBoxes[i];
//           colorLDFB[adit()].setVal(3.0, intersectBox, 0);
//         }
//       }
//     }


}


void scaleDownBoundaryCells(Vector<LDFB*>& hog,
                            const Vector<DBL>& DBLVector,
                            const Vector<Box>& boundingBox,
                            const Vector<Real>& DxVector,
                            const Vector<int>& refinementRatioVector,
                            const int baseLevel,
                            const int numLevels,
                            const int component)
{

  // now i need to go thru each cell and scale appropriately where appropriate
  int numberCoarseFineBoundaryFaces=0;
  int numberProblemBoundaryFaces=0;
  int numberCells;
  int CFhold=0;

  Vector<LDFB*> color(numLevels, NULL);
  createHOG(color, numLevels, DBLVector);

  // loop over levels
  for (int ilev=baseLevel; ilev<numLevels; ilev++)
  {

    printf(" tagging: ilev=%d  baseLevel=%d  numLevels=%d\n",
           ilev, baseLevel, numLevels);

    LDFB& levhog = *hog[ilev];
    const DBL currentDBL = levhog.disjointBoxLayout();

    IntVectSet coarseFineSkinSet;
    IntVectSet fineCoarseSkinSet;
    IntVectSet problemBoundarySkinSet;
    numberCoarseFineBoundaryFaces=0;
    numberProblemBoundaryFaces=0;

    // obtain the number of cells in current DBL
    numberCells = cellCount(currentDBL);

    // find the number of faces on boundary and construct a
    // set of all the cells on the problem boundary for this level.
    problemBoundaryCount(currentDBL,
                         boundingBox[ilev],
                         problemBoundarySkinSet,
                         numberProblemBoundaryFaces);

    // find the number of coarse-fine boundary faces and construct
    // a set of all cells on this level that border such a face.
    // ie all cells that touch a level finer than current level.
    // Only compute coarse-fine boundary when there is a finer level
    if (ilev < numLevels-1)
    {
      //cout << " looking for coarse-fine boundaries " << endl;
      LDFB& nextlevhog = *hog[ilev+1];
      const DBL nextDBL = nextlevhog.disjointBoxLayout();
      DBL nextDBL_coarse;
      coarsen(nextDBL_coarse, nextDBL, 2);  // will this always be 2?
      coarseFineBoundaryCount(nextDBL_coarse,
                              boundingBox[ilev],
                              coarseFineSkinSet,
                              numberCoarseFineBoundaryFaces);
    }

    // construct a set of all cells that border a fine-coarse boundary.
    // ie all cells that touch a level coarser than current level.
    // Only compute fine-coarse boundary when there is a coarser level
    if (ilev > baseLevel)
    {
      //cout << " looking for fine-coarse boundaries " << endl;
      fineCoarseBoundaryCount(currentDBL, boundingBox[ilev],
                              fineCoarseSkinSet);
    }

    printf("  level%d  numberCoarseFineBoundaryFaces=%d  numberProblemBoundaryFaces=%d  numberCells=%d\n",
           ilev, numberCoarseFineBoundaryFaces, numberProblemBoundaryFaces, numberCells);

    printf("   coarseFine.size=%d  fineCoarse.size=%d  problemBoundary.size=%d\n",
           coarseFineSkinSet.numPts(), fineCoarseSkinSet.numPts(),
           problemBoundarySkinSet.numPts());

    const Real KcoarseFine = (Real)numberCoarseFineBoundaryFaces/(Real)numberCells;
    const Real ND = refinementRatioVector[ilev]*pow(2.0,SpaceDim-2);
    const Real KfineCoarse = ND*CFhold/(Real)numberCells;
    const Real KproblemBoundary =
      (Real)numberProblemBoundaryFaces/(Real)numberCells * DxVector[ilev];
    CFhold = numberCoarseFineBoundaryFaces;

    printf("   KcoarseFine=%12.3e  KfineCoarse=%12.3e  KproblemBoundary=%12.3e\n",
           KcoarseFine, KfineCoarse, KproblemBoundary);

    //     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    // for debuggin/plotting purposes only:
    fillColorFAB(color[ilev], coarseFineSkinSet,
                 fineCoarseSkinSet, problemBoundarySkinSet);
    char fn2[80];
    sprintf(fn2, "color%dof%d.hdf5", ilev, numLevels);
    //writeLDFBname(color, fn);
    outputHDF5(color, DBLVector, boundingBox,
               refinementRatioVector, ilev+1, fn2);
    //     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


    IntVectSet local_tags;

    // loop over boxes in level
    //int boxCount=0;
    for (DataIterator dit=levhog.dataIterator(); dit.ok(); ++dit)
    {

      //printf("   box#%d\n", boxCount);
      //boxCount++;

      const Box& box = levhog.box(dit());
      //FArrayBox grad_fab (b, SpaceDim);
      FArrayBox& fab = levhog[dit()];

      // loop over each cell in box
      for (BoxIterator bit = BoxIterator(box); bit.ok(); ++bit)
      {
        //printf(" %f \n", fab(bit(),component));

        Real Ksmudge = 1.0;

        // if this cell is contained within two skin sets,
        // then use the smallest K?
        if (coarseFineSkinSet.contains(bit()))
        {
          Ksmudge = KcoarseFine;
        }

        if (fineCoarseSkinSet.contains(bit()))
        {
          if (KfineCoarse < Ksmudge) Ksmudge = KfineCoarse;
        }

        if (problemBoundarySkinSet.contains(bit()))
        {
          if (KproblemBoundary < Ksmudge) Ksmudge = KproblemBoundary;
        }

        //if ( fabs(Ksmudge*fab(bit(),component)) >= threshold)
        //{
        ////cout << bit() << endl;
        //local_tags |= bit();
        //}

        // WACTH OUT!  i'm changing values of hog here.  after
        // this loop it will be a scaled down hierarchy and will be
        // used for plotting purposes only
        fab(bit(),component) *= Ksmudge;

        //cout << " num local tags = " << local_tags.numPts() << endl;
      }

    } // end loop on boxes for this level

    //printf("   boxCount=%d\n", boxCount);

    char fn5[80];
    sprintf(fn5, "est_truncation_error_level%db.hdf5", numLevels-1);
    outputHDF5(hog, DBLVector, boundingBox,
               refinementRatioVector, numLevels, fn5);

  }

}

