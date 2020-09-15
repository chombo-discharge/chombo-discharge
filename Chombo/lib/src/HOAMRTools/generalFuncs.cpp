#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// generalFuncs.cpp
// petermc, 24 March 2003

#include <cmath>
#include <string>
using std::string;

#include "generalFuncs.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "SPMD.H"
#include "LoadBalance.H"
using std::cout;
using std::cerr;
using std::endl;

// ---------------------------------------------------------
Box
GenFuncs::boxFromVector(const Vector<int>&  a_ints,
                        int                 a_start) // default a_start == 0
{
  IntVect lo(D_DECL6(a_ints[a_start+0], a_ints[a_start+1], a_ints[a_start+2],
                     0, 0, 0));
  IntVect hi(D_DECL6(a_ints[a_start+3], a_ints[a_start+4], a_ints[a_start+5],
                     0, 0, 0));
  Box bx(lo, hi);
  return bx;
}


// ---------------------------------------------------------
int
GenFuncs::readBoxes(Vector< Vector<Box> >&     a_vectVectBoxes,
                    ParmParse&                 a_pp,
                    const ProblemDomain&       a_domain,
                    const int&                 a_maxGridSize,
                    const int&                 a_blockFactor,
                    const int&                 a_numlevels,
                    const Vector<int>&         a_vectRefRatio,
                    const bool&                a_verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int eekflag = 0;

  a_vectVectBoxes.resize(a_numlevels);
  const Box& domain0box = a_domain.domainBox();

  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
    {
      /*
        Level 0
      */
      Vector<Box>& vectBoxes0 = a_vectVectBoxes[0];
      // vectBoxes0 = Vector<Box>(1, domain0box);
      // Break up coarsest level if necessary.
      domainSplit(a_domain, vectBoxes0, a_maxGridSize, a_blockFactor);
      if (a_verbose >= 3)
        {
          pout() << "level 0: ";
          for (int ibox = 0; ibox < vectBoxes0.size(); ibox++)
            {
              pout() << vectBoxes0[ibox] << endl;
            }
        }

      if (a_numlevels > 1)
        {
          /*
            Level 1
          */
          Vector<int> level1domainIndices, level1boxIndices;
          // variables with defaults
          int level1boxcount = 0;

          Vector<Box>& vectBoxes1 = a_vectVectBoxes[1];
          if (a_pp.query("level_1_boxcount", level1boxcount))
            {
              a_pp.getarr("level_1_boxes", level1boxIndices,
                          0, 2*3*level1boxcount);
              a_pp.queryarr("level_1_domain", level1domainIndices,
                            0, 2*3);
            }
          Box domain1box = refine(domain0box, a_vectRefRatio[0]);
          eekflag = getBoxesFromIndices(vectBoxes1, domain1box,
                                        a_maxGridSize,
                                        a_blockFactor,
                                        level1boxcount,
                                        level1boxIndices,
                                        level1domainIndices);
          if (a_numlevels > 2)
            {
              /*
                Level 2
              */
              Vector<int> level2domainIndices, level2boxIndices;
              // variables with defaults
              int level2boxcount = 0;
              Vector<Box>& vectBoxes2 = a_vectVectBoxes[2];
              if (a_pp.query("level_2_boxcount", level2boxcount))
                {
                  a_pp.getarr("level_2_boxes", level2boxIndices,
                              0, 2*3*level2boxcount);
                  a_pp.queryarr("level_2_domain", level2domainIndices,
                                0, 2*3);
                }
              Box domain2box = refine(domain1box, a_vectRefRatio[1]);
              eekflag = getBoxesFromIndices(vectBoxes2, domain2box,
                                            a_maxGridSize,
                                            a_blockFactor,
                                            level2boxcount,
                                            level2boxIndices,
                                            level2domainIndices);
            }
        }
    }
  // Broadcast results to all the processors.
  broadcast(a_vectVectBoxes, uniqueProc(SerialTask::compute));
  return(eekflag);
}

#ifdef CH_USE_PYTHON
//---------------------------------------------------------
int
GenFuncs::readBoxes(Vector< Vector<Box> >&     a_vectVectBoxes,
                    PyParse&                 a_pp,
                    const ProblemDomain&       a_domain,
                    const int&                 a_maxGridSize,
                    const int&                 a_blockFactor,
                    const int&                 a_numlevels,
                    const Vector<int>&         a_vectRefRatio,
                    const bool&                a_verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int eekflag = 0;

  a_vectVectBoxes.resize(a_numlevels);
  const Box& domain0box = a_domain.domainBox();

  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
    {
      /*
        Level 0
      */
      Vector<Box>& vectBoxes0 = a_vectVectBoxes[0];
      // vectBoxes0 = Vector<Box>(1, domain0box);
      // Break up coarsest level if necessary.
      domainSplit(a_domain, vectBoxes0, a_maxGridSize, a_blockFactor);
      if (a_verbose >= 3)
        {
          pout() << "level 0: ";
          for (int ibox = 0; ibox < vectBoxes0.size(); ibox++)
            {
              pout() << vectBoxes0[ibox] << endl;
            }
        }

      if (a_numlevels > 1)
        {
          /*
            Level 1
          */
          Vector<int> level1domainIndices, level1boxIndices;
          // variables with defaults
          int level1boxcount = 0;

          Vector<Box>& vectBoxes1 = a_vectVectBoxes[1];
          if (a_pp.query("level_1_boxcount", level1boxcount))
          {
            a_pp.get("level_1_boxes", level1boxIndices.stdVector());
            a_pp.query("level_1_domain", level1domainIndices.stdVector());
          }
          Box domain1box = refine(domain0box, a_vectRefRatio[0]);
          eekflag = getBoxesFromIndices(vectBoxes1, domain1box,
                                        a_maxGridSize,
                                        a_blockFactor,
                                        level1boxcount,
                                        level1boxIndices,
                                        level1domainIndices);
          if (a_numlevels > 2)
            {
              /*
                Level 2
              */
              Vector<int> level2domainIndices, level2boxIndices;
              // variables with defaults
              int level2boxcount = 0;
              Vector<Box>& vectBoxes2 = a_vectVectBoxes[2];
              if (a_pp.query("level_2_boxcount", level2boxcount))
              {
                a_pp.get("level_2_boxes", level2boxIndices.stdVector());
                a_pp.query("level_2_domain", level2domainIndices.stdVector());
              }
              Box domain2box = refine(domain1box, a_vectRefRatio[1]);
              eekflag = getBoxesFromIndices(vectBoxes2, domain2box,
                                            a_maxGridSize,
                                            a_blockFactor,
                                            level2boxcount,
                                            level2boxIndices,
                                            level2domainIndices);
            }
        }
    }
  // Broadcast results to all the processors.
  broadcast(a_vectVectBoxes, uniqueProc(SerialTask::compute));
  return(eekflag);
}
#endif


// --------------------------------------------------------------
int
GenFuncs::getDiff(Vector<LevelData<FArrayBox>* >& a_vectDiff,
                  const Vector<LevelData<FArrayBox>* >& a_vectPhi1,
                  const Vector<LevelData<FArrayBox>* >& a_vectPhi2)
{
  int numlevels = a_vectPhi1.size();
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      const LevelData<FArrayBox>& phi1 = *a_vectPhi1[ilev];
      const LevelData<FArrayBox>& phi2 = *a_vectPhi2[ilev];
      const int nComp = phi1.nComp();
      CH_assert(nComp == phi2.nComp());
      const DisjointBoxLayout& grids1 = phi1.getBoxes();
      // const DisjointBoxLayout& grids2 = phi2.getBoxes();
      // CH_assert(grids1 == grids2);
      a_vectDiff[ilev] = new
        LevelData<FArrayBox>(grids1, nComp, IntVect::Zero);
      LevelData<FArrayBox>& diffLevel = *a_vectDiff[ilev];
      for (DataIterator dit(grids1.dataIterator()); dit.ok(); ++dit)
        {
          const FArrayBox& phi1Fab = phi1[dit];
          const FArrayBox& phi2Fab = phi2[dit];
          FArrayBox& diffFab = diffLevel[dit];
          diffFab.copy(phi1Fab);
          diffFab -= phi2Fab;
        }
    }
  return(0);
}


// --------------------------------------------------------------
int
GenFuncs::getMagnitude(Vector<LevelData<FArrayBox>* >& a_vectMag,
                       const Vector<LevelData<FArrayBox>* >& a_vectField)
{
  int numlevels = a_vectField.size();
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      const LevelData<FArrayBox>& field = *a_vectField[ilev];
      const int nComp = field.nComp();
      CH_assert(nComp == SpaceDim);
      const DisjointBoxLayout& grids = field.getBoxes();
      // const DisjointBoxLayout& grids2 = phi2.getBoxes();
      // CH_assert(grids1 == grids2);
      a_vectMag[ilev] = new
        LevelData<FArrayBox>(grids, 1, IntVect::Zero);
      LevelData<FArrayBox>& magLevel = *a_vectMag[ilev];
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit)
        {
          const FArrayBox& fieldFab = field[dit];
          FArrayBox& magFab = magLevel[dit];

          Box bx(magFab.box());
          for (BoxIterator bit(bx); bit.ok(); ++bit)
            {
              IntVect iv = bit();

              Real magPoint2 = 0.;
              for (int idim = 0; idim < SpaceDim; idim++)
                magPoint2 += fieldFab(iv, idim) * fieldFab(iv, idim);

              magFab(iv, 0) = sqrt(magPoint2);
            }
        }
    }
  return(0);
}


// --------------------------------------------------------------
int
GenFuncs::project2(Vector<LevelData<FArrayBox>* >& a_vectProj,
                   const Vector<LevelData<FArrayBox>* >& a_vectPhi,
                   const Vector<DisjointBoxLayout>& a_vectGrids)
{
  int numlevels = a_vectPhi.size();
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      const LevelData<FArrayBox>& phiLevel = *a_vectPhi[ilev];
      const int nComp = phiLevel.nComp();
      const DisjointBoxLayout& grids = phiLevel.getBoxes();
      DisjointBoxLayout gridsCoarse;
      coarsen(gridsCoarse, grids, 2);
      // a_vectProj[ilev] = new
      // LevelData<FArrayBox>(a_vectGrids[ilev], nComp, IntVect::Zero);
      LevelData<FArrayBox>* projPtr = new
        LevelData<FArrayBox>(gridsCoarse, nComp, IntVect::Zero);
      LevelData<FArrayBox>& proj = *projPtr;
      for (DataIterator dit(phiLevel.dataIterator()); dit.ok(); ++dit)
        {
          const FArrayBox& phiFab = phiLevel[dit];
          FArrayBox& projFab = proj[dit];
          Box projBox = projFab.box();

          for (BoxIterator bit(projBox); bit.ok(); ++bit)
            {
              IntVect ivc = bit();
              IntVect ivf = 2 * ivc;
              for (int ivar = 0; ivar < nComp; ivar++)
                projFab(ivc, ivar) = phiFab(ivf, ivar);
            }
        }
      // Unfortunately you have to do this copyTo business because
      // of some quirk that doesn't recognize the same indices
      // used in the two LevelDatas when DEBUG==TRUE.
      a_vectProj[ilev] = new
        LevelData<FArrayBox>(a_vectGrids[ilev], nComp, IntVect::Zero);
      proj.copyTo(proj.interval(),
                  *a_vectProj[ilev],
                  a_vectProj[ilev]->interval());
      delete projPtr;
    }
  return(0);
}


// --------------------------------------------------------------
int
GenFuncs::printErrorNorms(const string& a_prefix,
                          Real a_dxFine,
                          int a_deginterp,
                          Real a_normMaxerrFine,
                          Real a_norm1errFine,
                          Real a_norm2errFine,
                          Real a_normMaxerrCoarse,
                          Real a_norm1errCoarse,
                          Real a_norm2errCoarse)
{
  int dxinv = int((1.0 + 1.e-10)/a_dxFine);
  Real logTwo = log(2.);
  const char* cprefix = a_prefix.c_str();
  printf("%s %d & %d & $L_\\infty$ & 1/%-4d & %.5E & %.5E & %.2f \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv,
         a_normMaxerrFine, a_normMaxerrCoarse,
         log(a_normMaxerrCoarse/a_normMaxerrFine)/logTwo);

  printf("%s %d & %d & $L_1     $ & 1/%-4d & %.5E & %.5E & %.2f \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv,
         a_norm1errFine, a_norm1errCoarse,
         log(a_norm1errCoarse/a_norm1errFine)/logTwo);

  printf("%s %d & %d & $L_2     $ & 1/%-4d & %.5E & %.5E & %.2f \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv,
         a_norm2errFine, a_norm2errCoarse,
         log(a_norm2errCoarse/a_norm2errFine)/logTwo);
  return(0);
}


// --------------------------------------------------------------
int
GenFuncs::printDiffNorms(const string& a_prefix,
                         Real a_dxFine,
                         int a_deginterp,
                         Real a_normMaxdiff,
                         Real a_norm1diff,
                         Real a_norm2diff)
{
  int dxinv = int((1.0 + 1.e-10)/a_dxFine);
  int halfdxinv = int((1.0 + 1.e-10)/(0.5 * a_dxFine));
  const char* cprefix = a_prefix.c_str();
  printf("%s %d & %d & $L_\\infty$ & 1/%-4d & %.5E &             &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv, a_normMaxdiff);
  printf("%s %d & %d & $L_\\infty$ & 1/%-4d &             & %.5E &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, halfdxinv, a_normMaxdiff);
  printf("%s %d & %d & $L_1     $ & 1/%-4d & %.5E &             &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv, a_norm1diff);
  printf("%s %d & %d & $L_1     $ & 1/%-4d &             & %.5E &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, halfdxinv, a_norm1diff);
  printf("%s %d & %d & $L_2     $ & 1/%-4d & %.5E &             &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv, a_norm2diff);
  printf("%s %d & %d & $L_2     $ & 1/%-4d &             & %.5E &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, halfdxinv, a_norm2diff);
  return(0);
}


// --------------------------------------------------------------
int
GenFuncs::getBoxesFromIndices(Vector<Box>&        a_vectBoxes,
                              const Box&          a_domainBox,
                              const int&          a_maxGridSize,
                              const int&          a_blockFactor,
                              const int&          a_boxCount,
                              const Vector<int>&  a_boxIndices,
                              const Vector<int>&  a_domainIndices)
{
  if (a_boxCount == 0)
    {
      // Make the level cover the whole domain.
      domainSplit(a_domainBox, a_vectBoxes, a_maxGridSize, a_blockFactor);
    }
  else
    {
      Box subdomain = boxFromVector(a_domainIndices);
      // multiplier from subdomain to domainLevelBox
      int multiplier = 1;
      Box expandedDomain(subdomain);
      Box refinedExpandedDomain = refine(expandedDomain, 2);
      // Keep refining expandedDomain until it's just big enough to fit
      // inside a_domainBox.
      while (a_domainBox.contains(refinedExpandedDomain))
        {
          expandedDomain = refinedExpandedDomain;
          refinedExpandedDomain = refine(expandedDomain, 2);
          multiplier *= 2;
        }
      for (int ibox = 0; ibox < a_boxCount; ibox++)
        {
          Box baseBox = boxFromVector(a_boxIndices, 2*3*ibox);
          if (! subdomain.contains(baseBox))
            {
              cerr << "readGrids: level 1 box " << baseBox
                   << " not in " << subdomain << endl;
              return(1);
            }
          // Now scale thisBox up to size of the domain.
          Box thisBox = refine(baseBox, multiplier);
          a_vectBoxes.push_back(thisBox);
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          breakBoxes(a_vectBoxes, a_maxGridSize, idir);
        }
    }
  return (0);
}
