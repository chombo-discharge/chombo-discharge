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
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "MFBackwardEuler.H"
#include "MFLevelDataOps.H"
#include "NamespaceHeader.H"

void MFBackwardEuler::
resetAlphaAndBeta(const Real& a_alpha,
                  const Real& a_beta)
{
  Vector<MGLevelOp< LevelData<MFCellFAB> >* > ops = m_solver->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
    {
      TGAHelmOp< LevelData<MFCellFAB> >* helmop = (TGAHelmOp< LevelData<MFCellFAB> >*) ops[iop];
      helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
}

void MFBackwardEuler::setTime(Real a_time)
{
  Vector<MGLevelOp< LevelData<MFCellFAB> >* > ops = m_solver->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
    {
      TGAHelmOp< LevelData<MFCellFAB> >* helmop = (TGAHelmOp< LevelData<MFCellFAB> >*) ops[iop];
      helmop->setTime(a_time);
    }
}

/*****/
TGAHelmOp<LevelData<MFCellFAB> >*
MFBackwardEuler::
newOp(const ProblemDomain&                             a_indexSpace,
      const AMRLevelOpFactory<LevelData<MFCellFAB> >&  a_opFact)
{
  AMRLevelOpFactory<LevelData<MFCellFAB> >& opFact = (AMRLevelOpFactory<LevelData<MFCellFAB> >&) a_opFact;
  TGAHelmOp<LevelData<MFCellFAB> >* retval = (TGAHelmOp<LevelData<MFCellFAB> >*) opFact.AMRnewOp(a_indexSpace);
  return retval;
}
/*****/
MFBackwardEuler::
MFBackwardEuler(const RefCountedPtr<AMRMultiGrid< LevelData<MFCellFAB> > >& a_solver,
                const AMRLevelOpFactory<LevelData<MFCellFAB> > &            a_opFact,
                const ProblemDomain&                                        a_level0Domain,
                const Vector<int>&                                          a_refRat,
                int                                                         a_numLevels,
                int                                                         a_verbosity)
{
  m_verbosity = a_verbosity;
  m_level0Domain = a_level0Domain;
  m_refRat = a_refRat;
  m_solver  = a_solver;
  m_numLevels = a_numLevels;
  if (m_numLevels < 0)
    {
      m_numLevels = a_refRat.size();
    }

  m_ops.resize(m_numLevels);

  AMRLevelOpFactory<LevelData<MFCellFAB> >& opFact  =
    (AMRLevelOpFactory<LevelData<MFCellFAB> >&) a_opFact;
  ProblemDomain curDom = m_level0Domain;
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      m_ops[ilev] = RefCountedPtr<TGAHelmOp<LevelData<MFCellFAB> > >(newOp(curDom, opFact));
      curDom.refine(a_refRat[ilev]);
    }

  m_dataCreated = false;
}

/*****/
MFBackwardEuler::~MFBackwardEuler()
{
  for (int ilev = 0; ilev < m_rhst.size(); ilev++)
    {
      if (m_rhst[ilev] != NULL)
        {
          delete m_rhst[ilev];
          m_rhst[ilev] = NULL;
        }
    }
}
/***/
void
MFBackwardEuler::
createData(Vector<LevelData<MFCellFAB>* >&       a_source,
           int                                   a_lbase,
           int                                   a_lmax)
{
  m_dataCreated = true;
  m_rhst.resize(a_source.size(), NULL);
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_rhst[ilev] = new LevelData<MFCellFAB>();
      m_ops[ilev]->create(*m_rhst[ilev], *a_source[ilev]);
    }
}
/***/
void MFBackwardEuler::
oneStep(Vector<LevelData<MFCellFAB>* >&       a_phiNew,
        Vector<LevelData<MFCellFAB>* >&       a_phiOld,
        Vector<LevelData<MFCellFAB>* >&       a_source,
        const Real&                           a_dt,
        int                                   a_lbase,
        int                                   a_lmax,
        Real                                  a_told,
        bool                                  a_zeroPhi)
{
  if (!m_dataCreated)
    {
      createData(a_source, a_lbase, a_lmax);
    }

  if (m_verbosity >= 3)
    {
      pout() << "  MFBackwardEuler:: making rhs" << std::endl;
    }
  createEulerRHS(m_rhst, a_source, a_phiOld, a_lbase, a_lmax, a_dt);

  if (m_verbosity >= 3)
    {
      pout() << "  MFBackwardEuler:: solving" << std::endl;
    }
  Real tnew = a_told + a_dt;
  setTime(tnew);
  solveHelm(a_phiNew, m_rhst, a_lbase, a_lmax,a_dt, a_zeroPhi);
}
/***/
void MFBackwardEuler::
residual(Vector<LevelData<MFCellFAB>* >&       a_error,
         Vector<LevelData<MFCellFAB>* >&       a_phiNew,
         Vector<LevelData<MFCellFAB>* >&       a_phiOld,
         Vector<LevelData<MFCellFAB>* >&       a_source,
         const Real&                           a_dt,
         int                                   a_lbase,
         int                                   a_lmax,
         Real                                  a_told)
{
  if (!m_dataCreated)
    {
      createData(a_source, a_lbase, a_lmax);
    }

  createEulerRHS(m_rhst, a_source, a_phiOld, a_lbase, a_lmax, a_dt);
  Real tnew = a_told + a_dt;
  setTime(tnew);
  applyHelm(a_error, a_phiNew, a_lbase, a_lmax, -1.0, a_dt, false);
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      MFLevelDataOps::incr(*a_error[ilev], *m_rhst[ilev], -1.0);
      MFLevelDataOps::scale(*a_error[ilev], -1.0);
    }
}
/*******/
//fills a_ans = dt*kappa*a_source + kappa*acoef*phiOld
void MFBackwardEuler::
createEulerRHS(Vector<LevelData<MFCellFAB>* >&   a_ans,
               Vector<LevelData<MFCellFAB>* >&   a_rho,
               Vector<LevelData<MFCellFAB>* >&   a_phiOld,
               int                               a_lbase,
               int                               a_lmax,
               Real                              a_dt)

{
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      DisjointBoxLayout grids = a_phiOld[ilev]->disjointBoxLayout();
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          //this makes rhs  = 0
          (*a_ans[ilev])[dit()].setVal(0.);
          //this makes rhs = phiOld
          (*a_ans[ilev])[dit()] +=  ((*a_phiOld[ilev])[dit()]);
        }
      //this makes rhs = kappa*acoef*(phi^n)
      m_ops[ilev]->diagonalScale(*a_ans[ilev]);

      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          const MFCellFAB& rho  = (*a_rho[ilev])[dit()];
          int nphase = rho.numPhases();
          Vector<EBISBox> ebisboxv(nphase);
          Vector<int> nvar(nphase);
          for (int iphase=0; iphase<nphase; iphase++)
            {
              nvar[iphase] = rho.nComp(iphase);
              ebisboxv[iphase] = rho.getPhase(iphase).getEBISBox();
            }
          MFCellFAB scaleRho(ebisboxv, rho.box(), nvar);
          scaleRho.setVal(0.);

          scaleRho += rho;
          //this makes scalerho = kappa*a_rho
          MFLevelDataOps::kappaWeight(scaleRho);

          //this makes scalerho = dt*kappa*a_rho
          scaleRho *= a_dt;

          //this makes rhs = kappa*acoef*(phi^n) + dt*kappa*a_rho
          (*a_ans[ilev])[dit()] += scaleRho;
        }
    }
}
/*******/
void MFBackwardEuler::
solveHelm(Vector<LevelData<MFCellFAB>*>&       a_ans,
          Vector<LevelData<MFCellFAB>*>&       a_rhs,
          int                                  a_lbase,
          int                                  a_lmax,
          Real                                 a_dt,
          bool                                 a_zeroPhi)
{
  Real factor  = -a_dt;

  resetAlphaAndBeta(1.0, factor);

  m_solver->solveNoInit(a_ans, a_rhs, a_lmax, a_lbase, a_zeroPhi);

  if (m_solver->m_exitStatus==1 && m_verbosity>3)
    {
      pout() << "MFBackwardEuler:: WARNING: solver exitStatus == 1" << std::endl;
    }
}
/*******/
void MFBackwardEuler::
applyHelm(Vector<LevelData<MFCellFAB>*>&       a_ans,
          Vector<LevelData<MFCellFAB>*>&       a_phi,
          int                                  a_lbase,
          int                                  a_lmax,
          Real                                 a_mu,
          Real                                 a_dt,
          bool                                 a_homogeneousBC)
{
  Real factor  = a_mu*a_dt;

  resetAlphaAndBeta(1.0, factor);

  m_solver->computeAMROperator(a_ans, a_phi, a_lmax, a_lbase, a_homogeneousBC);

//   for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
//     {
//       m_ops[ilev]->scale(*a_ans[ilev], -1);
//     }
}

#include "NamespaceFooter.H"
