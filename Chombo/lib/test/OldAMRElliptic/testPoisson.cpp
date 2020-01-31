#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Dan Martin, Fri, Jan 14, 2000

#include "LayoutIterator.H"
#include "REAL.H"
#include "Vector.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Interval.H"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>

#include "parstream.H"
#include "BoxIterator.H"
#include "PoissProbF_F.H"
#include "LoadBalance.H"
#include "AMRIO.H"
#include "PoissonBC.H"
#include "BRMeshRefine.H"
#include "CoarseAverage.H"
#include  <iostream>
#include "LevelOp.H"
#include "CornerCopier.H"

#include "PoissonBC.H"
#include "ParmParse.H"
#include "AMRSolver.H"
#include "PoissonOp.H"

#include "SPMD.H"

#include "PoissProbF_F.H"
#include "ProblemDomain.H"

#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BCFunc.H"
#include "AMRPoissonOp.H"
#include "DebugDump.H"

#include "UsingNamespace.H"
using std::fstream;
using std::string;

// ------------------------------------------------------------
Real
computeSum(const LevelData<FArrayBox>& a_phi,
           const DisjointBoxLayout* a_finerGridsPtr,
           const int a_nRefFine, const Real a_dx,
           const Interval a_comps)
{

  Real sum, sumLevel;
  sum = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids, a_comps.size());

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
    {
      Box copyBox = temp[dit()].box();
      temp[dit()].copy(copyBox, tempComps,copyBox,
                     a_phi[dit()], a_comps);

    if (a_finerGridsPtr != NULL)
    {
      LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

      // now loop over fine boxes and set covered regions to 0
      for (litFine.reset(); litFine.ok(); ++litFine)
      {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= temp[dit()].box();

        if (!coveredBox.isEmpty())
        {
          temp[dit()].setVal(0.0, coveredBox, 0, tempComps.size());
        }
      }  // end loop over fine-grid boxes

    } // end if there is a finer level

    sumLevel = temp[dit()].sum(0, tempComps.size());
    sum += sumLevel;

  } // end loop over this level's grids

  Real scale = a_dx;
  if (SpaceDim == 2)
  {
    scale = a_dx*a_dx;
  }
  else if (SpaceDim == 3)
  {
    scale = a_dx*a_dx*a_dx;
  }

  sum *= scale;

  // do broadcast/gather thing here
#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

  sum = recv;
#endif

  return sum;

}
// --------------------------------------------------------------
Real
computeSum(const Vector<LevelData<FArrayBox>* >& a_phi,
           const Vector<int>& a_nRefFine,
           const Real a_dxCrse, const Interval a_comps,
           const int lBase)
{

  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real sum, sumLevel;

  sum = 0.0;

  // loop over levels
  for (int lev=lBase;  lev< numLevels; lev++)
  {
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
      }
      else
      {
        finerGridsPtr = NULL;
      }

      sumLevel = computeSum(thisPhi, finerGridsPtr, a_nRefFine[lev],
                            dxLevel, a_comps);
      sum += sumLevel;
      dxLevel = dxLevel/a_nRefFine[lev];
    }

  }

  return sum;
}

int outputData(const std::string filename,
               const Vector<LevelData<FArrayBox>* >& vectPhi,
               const Vector<LevelData<FArrayBox>* >& vectRhs,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<ProblemDomain>& vectDomain,
               const Vector<int>& vectRatio,
               Real dxCoarsest,
               int numlevels,
               bool verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int ncomp = vectPhi[0]->nComp();

  string phiNameBase("phi");
  string rhsNameBase("rhs");

  Vector<string> vectName(2*ncomp);
  for (int comp=0; comp<ncomp; comp++)
    {
      char labelString[80];
      sprintf(labelString, "%s%d", phiNameBase.c_str(), comp);
      string phiName(labelString);
      vectName[comp] = phiName;
      sprintf(labelString, "%s%d", rhsNameBase.c_str(), comp);
      string rhsName(labelString);
      vectName[ncomp+comp] = rhsName;
    }
  vectName.push_back("procID");
  Box domain = vectDomain[0].domainBox();
  //XXX -- not used
  //Real dx = 1.;
#ifdef CH_USE_HDF5
  Real dt = 1.;
  Real time = 1.;
#endif

  // problem is that invalid regions of phi contain whatever garbage
  // is left over from the solve; this can make viewing output difficult
  // since the garbage may affect the dynamic range in the plotfile.
  // fix this by doing an average down of finer-level solution(if it
  // exists) onto covered regions
  for (int ilev = numlevels-2; ilev >= 0; ilev--)
    {
      const DisjointBoxLayout& fineGrids = vectPhi[ilev+1]->getBoxes();
      CoarseAverage avgObject(fineGrids,
                              ncomp,
                              vectRatio[ilev]);

      avgObject.averageToCoarse(*vectPhi[ilev],
                                *vectPhi[ilev+1]);
    }

  //IntVect ghostVect = IntVect::Unit;
  IntVect ghostVect = IntVect::Zero;
  Vector<LevelData<FArrayBox>* > vectPhiAndRHS(numlevels, NULL);
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      vectPhiAndRHS[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], 2*ncomp+1, ghostVect);

      Interval phiInterval(0,ncomp-1);
      Interval rhsInterval(ncomp, 2*ncomp-1);

      vectPhi[ilev]->copyTo(vectPhi[ilev]->interval(),
                            *vectPhiAndRHS[ilev],
                            phiInterval);

      if (ghostVect != IntVect::Zero)
        {
          // do fab-by-fab copy to also get ghost cells
          LevelData<FArrayBox>& levelPhi = *vectPhi[ilev];
          LevelData<FArrayBox>& levelOut = *vectPhiAndRHS[ilev];
          DataIterator dit=levelPhi.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              levelOut[dit].copy(levelPhi[dit],0,0,ncomp);
            }
        }

      vectRhs[ilev]->copyTo(vectRhs[ilev]->interval(),
                            *vectPhiAndRHS[ilev],
                            rhsInterval);

      if (ghostVect != IntVect::Zero)
        {
          // use cornerCopier to make things a bit better
          CornerCopier cornerCop(vectGrids[ilev], vectGrids[ilev],
                                 vectDomain[ilev], IntVect::Unit, true);

          vectPhiAndRHS[ilev]->exchange(cornerCop);
        }

      DataIterator dit = vectPhiAndRHS[ilev]->dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        vectPhiAndRHS[ilev]->operator[](dit()).setVal(procID(), 2*ncomp);
      }
    }
#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhiAndRHS,
                        vectName,
                        domain,
                        dxCoarsest, dt, time,
                        vectRatio,
                        numlevels);
#endif

  for (int ilev = 0; ilev < numlevels; ilev++)
    delete vectPhiAndRHS[ilev];

  return 0;
}

/*
  Set RHS on hierarchy from input file
 */
int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
           Vector<Real>&                   vectDx,
           Vector<ProblemDomain>&          vectDomain,
           int numlevels, bool verbose)
{

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Real rhono, rno;
  int iprob;

  ParmParse ppLocal("main");
  bool useConstRHS = false;
  Real constRHS = -1.0;
  if (ppLocal.contains("constRhs"))
    {
      useConstRHS = true;
      ppLocal.get("constRhs", constRHS);
    }

  // problem 1 is better for dirichlet problem, while
  // problem 2 is for the periodic case
  if (vectDomain[0].isPeriodic())
    {
      iprob = 2;
    }
  else
    {
      iprob = 1;
    }

  rhono = 0.75;
  rno = 0.5;
  if (verbose)
    pout()
      << " rhono  = " << rhono
      << " rno  = " << rno
      << " iprob  = " << iprob
      << endl;
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      if (verbose)
        pout()
          << " ilev = " << ilev
          << " dom = " << vectDomain[ilev]
          << " dx  = " << vectDx[ilev]
        << endl;
      Real dxlev = vectDx[ilev];
      Box domlev =vectDomain[ilev].domainBox();
      LevelData<FArrayBox>& rhsLD = *vectRhs[ilev];
      DataIterator dit =  rhsLD.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& rhsFab = rhsLD[dit()];
          //kluge to get things the same in 0.2
          rhsFab.setVal(constRHS);
          /**/
          if (!useConstRHS)
            {
              FORT_GETRHSPOIS(CHF_FRA(rhsFab),
                              CHF_BOX(rhsFab.box()),
                              CHF_BOX(domlev),
                              CHF_CONST_REAL(dxlev),
                              CHF_CONST_REAL(rhono),
                              CHF_CONST_REAL(rno),
                              CHF_CONST_INT(iprob));
            }
        }
#ifdef CH_MPI
      MPI_Barrier(Chombo_MPI::comm);
#endif
    }
  return 0;
}
/*
  tag cells for refinement based on magnitude(RHS)
*/
void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels)
{
  for (int lev=baseLevel; lev!= numLevels; lev++)
    {
      IntVectSet local_tags;
      LevelData<FArrayBox> & levelRhs = *vectRHS[lev];
      DisjointBoxLayout level_domain = levelRhs.getBoxes();
      DataIterator dit = levelRhs.dataIterator();

      Real maxRHS = 0;

      maxRHS = norm(levelRhs, levelRhs.interval(), 0);

      Real tagVal = maxRHS * refine_thresh;

      // now loop through grids and tag cells where RHS > tagVal
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box thisBox = level_domain.get(dit());
          const FArrayBox& thisRhs = levelRhs[dit()];
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (abs(thisRhs(iv)) >= tagVal)
                local_tags |= iv;
            }
        } // end loop over grids on this level

      local_tags.grow(tags_grow);
      const Box& domainBox = vectDomain[lev].domainBox();
      local_tags &= domainBox;
      //  changed by bvs Fri Jan 16 14:40:45 PST 2004
      //
      //  somehow I had managed to not update the relationship between MeshRefine::regrid
      //  change and the AMRPoisson example code.
      tagVect[lev] = local_tags;
 //      int numTags = local_tags.numPts();
//       Vector<IntVectSet> all_tags;
//       const int dest_proc = uniqueProc(SerialTask::compute);
//       gather (all_tags, local_tags, dest_proc);
//       if (procID() == uniqueProc (SerialTask::compute)) {
//         for (int i=0; i< all_tags.size(); ++i)
//           {
//             tagVect[lev] |= all_tags[i];
//           }
//       }
//       numTags = tagVect[lev].numPts();
    } // end loop over levels
}

/*
  Set grid hierarchy from input file
 */
int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             Vector<ProblemDomain>&     vectDomain,
             Vector<Real>&              vectDx,
             Vector<int>&               vectRefRatio,
             int& numlevels, bool verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int max_level;
  Real fillRat = 0.77;
  std::vector<int> ancells(SpaceDim), refs;
  std::vector<Real>  prob_loa(SpaceDim);
  std::vector<Real>  prob_hia(SpaceDim);
  // set default to non-periodic
  std::vector<int> is_periodica(SpaceDim,0);
  bool is_periodic[SpaceDim];
  int temp;
  bool predefinedGrids = false;
  bool gridsFromFile = false;
  string gridFile;
  int maxboxsize = 32;
  Real refine_thresh = 0.5;

  ParmParse pp("main");
  pp.get("max_level", max_level);
  numlevels = max_level + 1;
  //kluge to get things the same in 0.2
  //max_level = 1;
  pp.query("fill_ratio", fillRat);
  pp.getarr("n_cell", ancells, 0, SpaceDim);
  pp.getarr("ref_ratio", refs, 0, numlevels);
  vectRefRatio = refs;
  pp.getarr("prob_lo",prob_loa,0,SpaceDim);
  pp.getarr("prob_hi",prob_hia,0,SpaceDim);
  pp.queryarr("is_periodic", is_periodica,0, SpaceDim);
  pp.query("maxboxsize",maxboxsize);

  numlevels = max_level + 1;

  temp = predefinedGrids;
  pp.query("presetGrids", temp);
  predefinedGrids = (temp == 1);

  if (pp.contains("grids_file"))
    {
      gridsFromFile = true;
      pp.get("grids_file", gridFile);
    }

  bool useConstantRHS = false;
  Real constRHS = -1.0;
  if (pp.contains("constRHS"))
    {
      useConstantRHS = true;
      pp.get("constRHS", constRHS);
    }

  pp.query("refinement_threshold", refine_thresh);

  // grid generation parameters
  int maxrat = 2;
  for (int irat = 0; irat < vectRefRatio.size(); irat++)
    maxrat = Max(vectRefRatio[irat], maxrat);
  int blockFactor = maxrat*4;
  if (maxboxsize/blockFactor < 1) blockFactor = maxboxsize;
  int nesting_radius = 1;

  IntVect ivlo = IntVect::Zero;
  IntVect ivhi;
  for (int idir = 0; idir < SpaceDim; idir++)
    ivhi[idir] = ancells[idir] - 1;

  Box basedom = Box(ivlo, ivhi);

  vectGrids.resize(numlevels);
  vectDomain.resize(numlevels);
  vectDx.resize(numlevels);

  vectDx.resize(numlevels,0.0);
  vectDx[0] = (prob_hia[0]-prob_loa[0])/ancells[0];
  for (int ilev = 1; ilev < numlevels; ilev++)
    {
      CH_assert(vectRefRatio[ilev-1] > 0);
      vectDx[ilev] = vectDx[ilev-1]/vectRefRatio[ilev-1];
    }

  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      is_periodic[dim] = (is_periodica[dim] == 1);
      if (is_periodic[dim] && verbose && procID() == 0)
        pout() << "Using Periodic BCs in direction: " << dim << endl;
    }

  vectDomain[0] = ProblemDomain(basedom,is_periodic);
  for (int ilev = 1;ilev < numlevels; ilev++)
    {
      vectDomain[ilev] = refine(vectDomain[ilev-1],vectRefRatio[ilev-1]);
    }

  int maxLevel = numlevels-1;
  Vector<Vector<Box> > newBoxes(numlevels);
  Vector<Vector<Box> > oldBoxes(numlevels);

  // use predefined grid configuration designed to test quadCFInterp configs
  if (predefinedGrids)
    {
      int nc = ancells[0];
      int nmi =  nc    / 2; // 16
      int nqu =  nc    / 4; //  8
      int ntf = (nc*3) / 4; // 24
#if CH_SPACEDIM > 1
      int nte = (nc*3) / 8; // 12
#endif
      int nfe = (nc*5) / 8; // 20

      Box boxf1, boxf2, boxf3, boxf4;

      boxf1.define(IntVect(D_DECL(0, nqu,nqu)),
                   IntVect(D_DECL(nmi-1,ntf-1,ntf-1)));

      boxf2.define(IntVect(D_DECL(nmi,nte,nte)),
                   IntVect(D_DECL(ntf-1,nfe-1,nfe-1)));

      boxf3.define(IntVect(D_DECL(nqu,0,0)  ),
                   IntVect(D_DECL(nfe-1,nqu-1,nqu-1)));

      boxf4.define(IntVect(D_DECL(nfe,nqu,nqu)),
                   IntVect(D_DECL(nc-1,nte-1,nte-1)));


      //comment out for kluge
      /**/
      IntVectSet tags;

      tags |= boxf2;
      tags |= boxf3;
      tags |= boxf4;
      tags |= boxf1;

      for (int ilev = 0; ilev <numlevels; ilev++)
        {
          oldBoxes[ilev].push_back(vectDomain[ilev].domainBox());
        }
      int baseLevel = 0;
      int topLevel  = numlevels - 2;
      int eekflag = 0;
      if (topLevel >= 0)
        {
          int newFinestLev;
          BRMeshRefine meshrefine(vectDomain[0],vectRefRatio, fillRat,
                                  blockFactor, nesting_radius, maxboxsize);
          newFinestLev = meshrefine.regrid(newBoxes, tags, baseLevel,
                                           topLevel, oldBoxes);
        }
      else
        {
          newBoxes = oldBoxes;
        }

      Vector< Vector<int> > procAssign;
      Real effRatio = 0.75;
      Vector< Vector<long> > loads(numlevels);
      for (int ilev = 0; ilev <numlevels; ilev++)
        {
          loads[ilev].resize(newBoxes[ilev].size());
          for (int ibox = 0; ibox < newBoxes[ilev].size() ; ibox++)
            {
              loads[ilev][ibox] = newBoxes[ilev][ibox].numPts();
            }
        }
      LoadBalance(procAssign, effRatio, newBoxes, loads, vectRefRatio);

      if (eekflag != 0)
        {
          pout() << "setGrids: loadBalance returned error code "
               << eekflag << endl;
          return(1);
        }
      for (int ilev = 0; ilev <numlevels; ilev++)
        {
          vectGrids[ilev].define(newBoxes[ilev], procAssign[ilev],
                                 vectDomain[ilev]);
          vectGrids[ilev].close();
        }
    }
  else if (gridsFromFile)
    {
      Vector<Vector<Box> > amrGrids(numlevels);
      Vector<Vector<int> > amrProcAssign(numlevels);

#ifdef CH_MPI
      MPI_Barrier(Chombo_MPI::comm);
      if (procID() == 0)
        {
          // additional check info for MPI case
          int numProcs = numProc();
          int maxProc = 0;
#endif
      // read grids (and proc assignments, if MPI) from file
      ifstream is(gridFile.c_str(), ios::in);
      if (is.fail())
        {
          pout() << "Cannot open grids file " << gridFile << endl;
        }
#ifdef CH_MPI
      if (verbose)
        {
          pout() << "procID: " << procID() << "  opening gridfile \n"
                 << endl;
        }
#endif
      // format of file -- number of levels, then for each level starting with
      // level 1, number of grids on level, list of boxes (along with proc
      // assignments if mpi)
      int in_numLevels;
      while (is.get() != '=');
      is >> in_numLevels;
      if (verbose)
        {
          pout() << "NumLevels in grids file = " << in_numLevels << endl;
        }

      numlevels = min(in_numLevels, numlevels);
      while (is.get() != '\n');
      amrGrids.resize(numlevels);
      int ngrid;
      for (int lev=0; lev<numlevels; lev++)
        {
          while (is.get() != ':');
          is >> ngrid;
          if (verbose)
            {
              pout() << "level " << lev << " numGrids = " << ngrid << endl;
              pout() << "grids = " ;
            }
          while (is.get() != '\n');
          amrGrids[lev].resize(ngrid);
          amrProcAssign[lev].resize(ngrid);
          Box bx;
          int this_proc=0;

          for (int i=0; i<ngrid; i++)
            {
              is >> bx;
#ifdef CH_MPI
              while (is.get() != '[');
              is >> this_proc;
              if (this_proc >= numProcs)
                {
                  pout() << "invalid processor assignment: assigned "
                         << this_proc << " for a "
                         << numProcs << " processor run" << endl;
                  MayDay::Error("Aborting");
                }
              if (this_proc > maxProc) maxProc = this_proc;
#endif
              // advance to next box
              while (char ch = is.get())
                {
                  if (ch == '#') break;
                  if (ch == '\n') break;
                }

              // assume that if we're defining the grids, we know what
              // we're doing (so comment out this test)
              /*
              // quick check on box size
              Box bxRef(bx);
              if (bxRef.longside() > maxboxsize)
              {
                pout() << "Grid " << bx << " too large" << endl;
                MayDay::Error();
              }
              */

              if (verbose)
              {
                pout() << bx << endl;
              }
              amrGrids[lev][i] = bx;
              amrProcAssign[lev][i] = this_proc;
            } // end loop over boxes on this level
        } // end loop over levels
#ifdef CH_MPI

        // warn if we're not using all processors
        if (maxProc < numProcs-1)
          {
            MayDay::Warning("Not using all processors in read-in grids file");
            pout() << "max procssor in grids file = " << maxProc
                   << "num processors = " << numProcs << endl;
          }
        } // end if procID = 0

      broadcast(amrGrids,0);
      broadcast(amrProcAssign,0);
#endif
      // now define disjointBoxLayouts
      for (int lev=0; lev<amrGrids.size(); lev++)
        {
          const DisjointBoxLayout newDBL(amrGrids[lev],
                                         amrProcAssign[lev],
                                         vectDomain[lev]);
          vectGrids[lev] = newDBL;
        }
    }
  else
    {
      // determine grids dynamically, based on grad(RHS)
      // will need temp storage for RHS
      Vector<LevelData<FArrayBox>* > vectRHS(maxLevel+1,NULL);
      int ncomps = 1;

      // define base level first
      Vector< Vector<int> > procAssign(maxLevel+1);
      domainSplit(vectDomain[0], oldBoxes[0], maxboxsize, blockFactor);
      procAssign[0].resize(oldBoxes[0].size());
      LoadBalance(procAssign[0],oldBoxes[0]);

      vectGrids[0].define(oldBoxes[0],procAssign[0],vectDomain[0]);

      vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], ncomps,
                                            IntVect::Zero);

      int topLevel = 0;

      bool moreLevels = (maxLevel > 0);

      // create grid generation object
      BRMeshRefine meshrefine(vectDomain[0], vectRefRatio, fillRat,
                              blockFactor, nesting_radius, maxboxsize);

      while (moreLevels)
        {
          // default is moreLevels = false
          // (only repeat loop in the case where a new level
          // is generated which is still less than maxLevel)
          moreLevels = false;

          int baseLevel = 0;
          int oldTopLevel = topLevel;

          // now initialize RHS for this existing hierarchy
          setRHS(vectRHS, vectDx, vectDomain, topLevel+1, verbose);

          Vector<IntVectSet> tagVect(topLevel+1);
          int tags_grow = 1;
          tagCells(vectRHS, tagVect, vectDx, vectDomain, refine_thresh,
                   tags_grow, baseLevel, topLevel+1);

//           if (procID() == uniqueProc(SerialTask::compute))
//             {
//               if (topLevel >= 0 && ! tagVect[topLevel].isEmpty())
//                 {

                  int new_finest = meshrefine.regrid(newBoxes, tagVect,
                                                     baseLevel,
                                                     topLevel, oldBoxes);

                  if (new_finest > topLevel)
                    topLevel++;

       //          } // end if there was a new level generated
//             } // end if proc is serial node
//           broadcast(topLevel, uniqueProc(SerialTask::compute) );
//           broadcast(newBoxes, uniqueProc(SerialTask::compute) );

          oldBoxes = newBoxes;

          //  no need to do this for the base level (already done)
          for (int lev=1; lev<= topLevel; lev++)
            {
              // do load balancing
              procAssign[lev].resize(newBoxes[lev].size());
              LoadBalance(procAssign[lev], newBoxes[lev]);
              const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                             vectDomain[lev]);
              vectGrids[lev] = newDBL;
              delete vectRHS[lev];
              vectRHS[lev] = new LevelData<FArrayBox>(vectGrids[lev], ncomps,
                                                      IntVect::Zero);
            } // end loop over levels for initialization

          // figure out whether we need another pass through grid generation
          if ((topLevel<maxLevel) && (topLevel > oldTopLevel))
            moreLevels = true;

        } // end while moreLevels loop
      // clean up temp storage
      for (int ilev=0; ilev <vectRHS.size(); ilev++)
        {
          if (vectRHS[ilev] != NULL)
            {
              delete vectRHS[ilev];
              vectRHS[ilev] = NULL;
            }
        }

    } // end if grids generated dynamically

  return 0;
}

static IntVect bclo = IntVect::Zero;
static IntVect bchi = IntVect::Zero;

extern "C"
{
  void Parabola_neum(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    switch (*dir)
    {
    case 0:
      a_values[0]=2*(*side)*(pos[0]);
      return;
    case 1:
      a_values[0]=2*(*side)*pos[1];
    return;
    case 2:
      a_values[0]=2*(*side)*pos[2];
      return;
    default:
      MayDay::Error("no such dimension");
    };
  }

  void Parabola_diri(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    a_values[0] = D_TERM((pos[0])*(pos[0]),+pos[1]*pos[1],+pos[2]*pos[2]);
  }

  void zeroFunc(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
  {
    a_values[0] = 0;
  }

  void FuncBC(FArrayBox& a_state,
              const Box& valid,
              const ProblemDomain& a_domain,
              Real a_dx,
              bool a_homogeneous)
  {
    const Box& d = a_domain.domainBox();
    for (int i=0; i<CH_SPACEDIM; ++i)
      {
        if ( valid.loVect()[i] == d.loVect()[i])
        {
          if (bclo[i] == 0)
            {
              DiriBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     zeroFunc,
                     i,
                     Side::Lo);
            }
          else
          {
            NeumBC(a_state,
                   valid,
                   a_dx,
                   a_homogeneous,
                   zeroFunc,
                   i,
                   Side::Lo);
          }
        }
        if (valid.hiVect()[i] == d.hiVect()[i])
        {
          if (bchi[i] == 0)
            {
              DiriBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     zeroFunc,
                     i,
                     Side::Hi);
            }
          else
            {
              NeumBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     zeroFunc,
                     i,
                     Side::Hi);
            }
        }
      }
  }

}

/*
  set domain boundary conditions from input file
  note that periodic BC's will override whichever physical
  boundary conditions are set here.
 */
int
setDomainBC(DomainGhostBC& domghostbc, bool verbose, int ncomp)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  std::vector<int> ibclo;
  std::vector<int> ibchi;
  ParmParse pp("main");
  pp.getarr("bc_lo", ibclo, 0, SpaceDim);
  pp.getarr("bc_hi", ibchi, 0, SpaceDim);

  Interval comps(0,ncomp-1);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //lo bc
      {
        NeumannBC neumbc(idir, Side::Lo,comps);
        DirichletBC dircbc(idir, Side::Lo,comps);
        if (ibclo[idir] == 0)
          {
            bclo[idir] = 0;
            if (verbose&& procID()==0)
              pout() << "Dirichlet bcs in low direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(dircbc);
          }
        else if (ibclo[idir] == 1)
          {
            bclo[idir]=1;
            if (verbose&& procID()==0)
              pout() << "neumann bcs in low direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(neumbc);
          }
        else
          {
            if (verbose)
              pout() << "setDomainBC:: bogus input bc_lo flag" << endl;
            return(1);
          }
      }
      //hi bc
      {
        NeumannBC neumbc(idir, Side::Hi,comps);
        DirichletBC dircbc(idir, Side::Hi,comps);
        if (ibchi[idir] == 0)
          {
            bchi[idir] = 0;
            if (verbose&& procID()==0)
              pout() << "Dirichlet bcs in high direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(dircbc);
          }
        else if (ibchi[idir] == 1)
          {
            bchi[idir]=1;
            if (verbose&& procID()==0)
              pout() << "neumann bcs in high direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(neumbc);
          }
        else
          {
            if (verbose)
              pout() << "setDomainBC:: bogus input bc_hi flag" << endl;
            return(2);
          }
      }
    }
  return(0);
}

/// Global variables for handling output:
static const char* pgmname = "testPoisson" ;
//static const char* indent = "   ";  //not used
static const char* indent2 = "      " ;
static bool verbose = false ;

///
// Parse the standard test options (-v -q) out of the command line.
///
void
parseTestOptions( int* argc ,char* argv[] )
{
  int i = 1;
  while ( i < *argc )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // move the rest of the arguments down by one
        for ( int j = i ; j+1 < *argc ; ++j ) argv[j] = argv[j+1] ;
        --(*argc);
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // move the rest of the arguments down by one
        for ( int j = i ; j+1 < *argc ; ++j ) argv[j] = argv[j+1] ;
        --(*argc);
      }
      else
      {
        ++i;
      }
    }
    else
    {
      ++i;
    }
  }
  return ;
}


int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //scoping trick
  {
  barrier();

  parseTestOptions( &argc ,argv ) ;

  char* in_file;
  char in_default[100];

  if (argc < 2)
  {
    sprintf(in_default,"inputs.testPoisson.%dd",CH_SPACEDIM);
    in_file = in_default;
  }
  else
  {
    in_file = argv[1];
  }

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << endl ;

  ParmParse  pp(argc-2,argv+2,NULL,in_file);

  int lbase;

  int iverbose;
  ParmParse pp2("main");
  pp2.get("lbase", lbase);
  pp2.get("verbose", iverbose);
  verbose |= (iverbose <= 1);
  int ncomp = 1;

  DomainGhostBC domghostbc;
  int eekflag = setDomainBC(domghostbc, verbose,ncomp);
  if (eekflag != 0)
    {
      MayDay::Error("error: setDomainBC returned error code",eekflag);
    }
  PoissonOp levelop;
  levelop.setDomainGhostBC(domghostbc);

  int numlevels;
  Vector<DisjointBoxLayout> vectGrids;
  Vector<ProblemDomain> vectDomain;
  Vector<int> vectRefRatio;
  Vector<Real> vectDx;
  eekflag = setGrids(vectGrids, vectDomain, vectDx,
                     vectRefRatio, numlevels, verbose);
  if (eekflag != 0)
    {
      MayDay::Error("error: setGrids returned error code",eekflag);
    }

  Vector<LevelData<FArrayBox>* > phi(numlevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numlevels, NULL);
  Vector<LevelData<FArrayBox>* > phi_new(numlevels, NULL);

  bool constPhi = false;
  Real initialPhi = 0.0;
  if (pp2.contains("initialGuess"))
    {
      constPhi = true;
      pp2.get("initialGuess", initialPhi);
    }

  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      phi[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Unit);
      rhs[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Zero);
      phi_new[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Unit);

      LevelData<FArrayBox>& phifabs= *phi[ilev];
      LevelData<FArrayBox>& p      = *phi_new[ilev];
      Real rval = Real(3*ilev);
      for (DataIterator dit = phifabs.dataIterator(); dit.ok(); ++dit)
      {
        rval += 1.0;
        if (constPhi)
          {
            // option to set phi to constant-value initial guess
            phifabs[dit()].setVal(initialPhi);
            p[dit].setVal(initialPhi);

          }
        else
          {
            phifabs[dit()].setVal(rval);
            p[dit].setVal(rval);
          }
      }
    }
  eekflag = setRHS(rhs, vectDx, vectDomain, numlevels, verbose);
  if (eekflag != 0)
    {
      MayDay::Error("error: setRHS returned error code",eekflag);
    }

  // now compute sum(RHS) as a diagnostic
  // (should be 0 for all-Neumann or periodic BC's for solvability).
  for (int comp=0; comp<ncomp; comp++)
    {
      Interval comps(comp,comp);
      Real sumRHS = computeSum(rhs,vectRefRatio,
                               vectDx[0], comps, 0);

      if (verbose && procID()==0)
        pout() << "Component " << comp << ",  Sum(RHS) = " << sumRHS << endl;

    } // end loop over components

  RealVect coarseSpacing(IntVect::Unit);
  coarseSpacing*=vectDx[0];

  AMRPoissonOpFactory factory;
  factory.define(vectDomain[0], vectGrids,
                 vectRefRatio, coarseSpacing[0],
                 FuncBC);
  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact =
    (  AMRLevelOpFactory<LevelData<FArrayBox> >&) factory;

  BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = 0;
  AMRMultiGrid<LevelData<FArrayBox> > newSolver;

  newSolver.define(vectDomain[0],
                   castFact,
                   &bottomSolver,
                   numlevels);

  AMRSolver amrSolver(vectGrids,  vectDomain,
                      vectDx,     vectRefRatio,
                      numlevels, lbase, &levelop,
                      ncomp);
  int maxiter;
  pp2.get("max_iterations", maxiter);


  amrSolver.setNumSmoothUp(4);

  amrSolver.setNumSmoothDown(3);

  amrSolver.setVerbose(verbose);
  pp2.get("max_iterations", maxiter);
  amrSolver.setMaxIter(maxiter);
  newSolver.m_iterMax = maxiter;

  int numvbot;
  pp2.get("num_v_cycles_bottom", numvbot);
  amrSolver.setNumVCyclesBottom(numvbot);
  Real solverTol = 1.0e-9;
  if (pp2.contains("solverTol"))
  {
    pp2.get("solverTol", solverTol);
    amrSolver.setTolerance(solverTol);
    newSolver.m_eps = solverTol;
  }


  int nsmooth = 4;
  if (pp2.contains("numSmoothUp"))
    {
      int numSmoothUp;
      pp2.get("numSmoothUp", numSmoothUp);
      amrSolver.setNumSmoothUp(numSmoothUp);
      newSolver.m_post = numSmoothUp;
      pout() << "numSmoothUp = " << numSmoothUp << endl;
      nsmooth = numSmoothUp;
    }

  if (pp2.contains("numSmoothDown"))
    {
      int numSmoothDown;
      pp2.get("numSmoothDown", numSmoothDown);
      amrSolver.setNumSmoothDown(numSmoothDown);
      newSolver.m_pre = numSmoothDown;
      pout() << "numSmoothDown = " << numSmoothDown << endl;
    }

  int nummg = 1;
  newSolver.setSolverParameters(nsmooth, nsmooth, nsmooth,
                                nummg, maxiter,
                                solverTol, solverTol, solverTol);

  if (pp2.contains("normType"))
    {
      int normType = 0;
      pp2.query("normType", normType);
      amrSolver.setNormType(normType);
    }

  // Add the ability to loop over the solve (the same solve) for benchmarking
  //  purposes only.  (ndk)  note that it will default to one time thru...
  int solve_iterations = 1;
  pp2.get("iterations", solve_iterations);
  for (int i=0; i<solve_iterations; i++)
    {
      newSolver.solve(phi_new, rhs, numlevels-1, 0);
      amrSolver.solveAMR(phi, rhs, false);
  }

//XXX should check that final residual is reasonable <dbs>

  eekflag = outputData("poisson_old.hdf5",phi, rhs, vectGrids, vectDomain,
                       vectRefRatio, vectDx[0], numlevels, verbose);

  if (eekflag != 0)
    {
      MayDay::Error("error: 1st outputData returned error code",eekflag);
    }
  eekflag = outputData("poisson_new.hdf5",phi_new, rhs, vectGrids, vectDomain,
                       vectRefRatio, vectDx[0], numlevels, verbose);
  if (eekflag != 0)
    {
      MayDay::Error("error: 2nd outputData returned error code",eekflag);
    }

  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      delete phi[ilev];
      delete rhs[ilev];
      delete phi_new[ilev];
    }
  }//end scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}
