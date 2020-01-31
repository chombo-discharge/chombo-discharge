#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelHelmholtzSolver.H"
#include "NamespaceHeader.H"

// default constructor
LevelHelmholtzSolver::LevelHelmholtzSolver()
{
  setDefaultValues();
}

// destructor
LevelHelmholtzSolver::~LevelHelmholtzSolver()
{
  clearMemory();
}

LevelHelmholtzSolver::LevelHelmholtzSolver(const DisjointBoxLayout& a_grids,
                                           const DisjointBoxLayout* a_baseGrids,
                                           const ProblemDomain& a_domain,
                                           Real a_dxLevel,
                                           int a_nRefCrse,
                                           const BaseHelmholtzOp* const a_opin,
                                           int a_ncomp)
{

  setDefaultValues();
  LevelSolver::define(a_grids, a_baseGrids, a_domain, a_dxLevel,
                      a_nRefCrse, a_opin, a_ncomp);
}


LevelHelmholtzSolver::LevelHelmholtzSolver(const DisjointBoxLayout& a_grids,
                                           const DisjointBoxLayout* a_baseGrids,
                                           const Box& a_domain,
                                           Real a_dxLevel,
                                           int a_nRefCrse,
                                           const BaseHelmholtzOp* const a_opin,
                                           int a_ncomp)
{

  setDefaultValues();
  ProblemDomain physdomain(a_domain);
  LevelSolver::define(a_grids, a_baseGrids, physdomain, a_dxLevel,
                      a_nRefCrse, a_opin, a_ncomp);
}


void
LevelHelmholtzSolver::define(const DisjointBoxLayout& a_grids,
                             const DisjointBoxLayout* a_baseGrids,
                             const Box& a_domain,
                             Real a_dxLevel,
                             int a_nRefCrse,
                             const BaseHelmholtzOp* const a_opin,
                             int a_ncomp)
{
  CH_TIME("LevelHelmholtzSolver::define");
  ProblemDomain physdomain(a_domain);
  LevelSolver::define(a_grids, a_baseGrids, physdomain, a_dxLevel,
                      a_nRefCrse, a_opin, a_ncomp);
}


void
LevelHelmholtzSolver::define(const DisjointBoxLayout& a_grids,
                             const DisjointBoxLayout* a_baseGrids,
                             const ProblemDomain& a_domain,
                             Real a_dxLevel,
                             int a_nRefCrse,
                             const BaseHelmholtzOp* const a_opin,
                             int a_ncomp)
{
  LevelSolver::define(a_grids, a_baseGrids, a_domain, a_dxLevel,
                      a_nRefCrse, a_opin, a_ncomp);
}


// this is just here to ensure that we don't try to get things going
// with a non-helmholtzOp
void
LevelHelmholtzSolver::define(const DisjointBoxLayout& a_grids,
                             const DisjointBoxLayout* a_baseGrids,
                             const Box& a_domain,
                             Real a_dxLevel,
                             int a_nRefCrse,
                             const LevelOp* const a_opin,
                             int a_ncomp)
{
  MayDay::Error("Tried to define a LevelHelmholtzSolver with a non-BaseHelmholtzOp");

}

// this is just here to ensure that we don't try to get things going
// with a non-helmholtzOp
void
LevelHelmholtzSolver::define(const DisjointBoxLayout& a_grids,
                             const DisjointBoxLayout* a_baseGrids,
                             const ProblemDomain& a_domain,
                             Real a_dxLevel,
                             int a_nRefCrse,
                             const LevelOp* const a_opin,
                             int a_ncomp)
{
  MayDay::Error("Tried to define a LevelHelmholtzSolver with a non-BaseHelmholtzOp");

}


// allows coefficient to be changed on the fly -- multiply beta by a_scale
void
LevelHelmholtzSolver::scaleBeta(Real a_scale)
{

  CH_assert (m_levelOpPtr != NULL);

  BaseHelmholtzOp* helmOpPtr = (BaseHelmholtzOp*) m_levelOpPtr;

  if (helmOpPtr != NULL)
    {
      helmOpPtr->scaleBeta(a_scale);
    }
  else
    {
      MayDay::Error("Bad levelOp->helmholtzOp conversion!");
    }

  LevelMG* levelmgPtr = &m_levelMG;

  while (levelmgPtr != NULL)
    {
      helmOpPtr = (BaseHelmholtzOp*) levelmgPtr->levelOpPtr();
      if (helmOpPtr != NULL)
        {
          helmOpPtr->scaleBeta(a_scale);
        }
      else
        {
          MayDay::Error("Bad levelOp->helmholtzOp conversion!");
        }

      levelmgPtr = levelmgPtr->lCoarsePtr();
    }
}


/// allows BC to be changed on the fly
void
LevelHelmholtzSolver::setBC(const DomainGhostBC& a_BC)
{
  CH_assert (m_levelOpPtr != NULL);

  BaseHelmholtzOp* helmOpPtr = (BaseHelmholtzOp*) m_levelOpPtr;

  if (helmOpPtr != NULL)
    {
      helmOpPtr->setDomainGhostBC(a_BC);
    }
  else
    {
      MayDay::Error("Bad levelOp->helmholtzOp conversion!");
    }

  LevelMG* levelmgPtr = &m_levelMG;

  while (levelmgPtr != NULL)
    {
      helmOpPtr = (BaseHelmholtzOp*) levelmgPtr->levelOpPtr();
      if (helmOpPtr != NULL)
        {
          helmOpPtr->setDomainGhostBC(a_BC);
        }
      else
        {
          MayDay::Error("Bad levelOp->helmholtzOp conversion!");
        }

      levelmgPtr = levelmgPtr->lCoarsePtr();
    }
}

#include "NamespaceFooter.H"

