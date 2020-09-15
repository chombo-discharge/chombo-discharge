#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FORT_PROTO.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CH_OpenMP.H"
#include "AMRMultiGrid.H"
#include "Misc.H"
#include "computeSum.H"

#include "AMRNonLinearPoissonOp.H"
#include "AMRNonLinearPoissonOpF_F.H"
#include "CCProjectorF_F.H"
#include "MACProjectorF_F.H"

#include "NamespaceHeader.H"

int AMRNonLinearPoissonOp::s_exchangeMode = 1; // 1: no overlap (default); 0: ...
//int AMRNonLinearPoissonOp::s_relaxMode = 0;
int AMRNonLinearPoissonOp::s_relaxMode = 1; // 1: GSRB (not implemented); 4: Jacobi(not implemented), 6: GS (not red black)
int AMRNonLinearPoissonOp::s_maxCoarse = 2;

// ---------------------------------------------------------
static void
amrpgetMultiColors(Vector<IntVect>& a_colors)
{

#if CH_SPACEDIM==2
  a_colors.resize(4);
  a_colors[0] = IntVect::Zero;             // (0,0)
  a_colors[1] = IntVect::Unit;             // (1,1)
  a_colors[2] = IntVect::Zero + BASISV(1); // (0,1)
  a_colors[3] = IntVect::Zero + BASISV(0); // (1,0)
#elif CH_SPACEDIM==3
  a_colors.resize(8);
  a_colors[0] = IntVect::Zero;                         // (0,0,0)
  a_colors[1] = IntVect::Zero + BASISV(0) + BASISV(1); // (1,1,0)
  a_colors[2] = IntVect::Zero + BASISV(1) + BASISV(2); // (0,1,1)
  a_colors[3] = IntVect::Zero + BASISV(0) + BASISV(2); // (1,0,1)
  a_colors[4] = IntVect::Zero + BASISV(1);             // (0,1,0)
  a_colors[5] = IntVect::Zero + BASISV(0);             // (1,0,0)
  a_colors[6] = IntVect::Zero + BASISV(2);             // (0,0,1)
  a_colors[7] = IntVect::Unit;                         // (1,1,1)
#endif
}

// ---------------------------------------------------------
/** full define function for AMRLevelOp with both coarser and finer levels */
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                          const DisjointBoxLayout& a_gridsFiner,
                          const DisjointBoxLayout& a_gridsCoarser,
                          const Real&              a_dxLevel,
                          int                      a_refRatio,
                          int                      a_refRatioFiner,
                          const ProblemDomain&     a_domain,
                          BCHolder                 a_bc,
                          const Copier&            a_exchange,
                          const CFRegion&          a_cfregion,
                          const int                a_nComp)

{
  CH_TIME("AMRNonLinearPoissonOp::define1");

  amrpgetMultiColors(m_colors);

  this->define(a_grids, a_gridsCoarser, a_dxLevel, a_refRatio, a_domain, a_bc,
               a_exchange, a_cfregion, a_nComp);
  m_refToFiner = a_refRatioFiner;

  if (a_gridsFiner.isClosed())
    {
      ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
      m_levfluxreg.define(a_gridsFiner,
                          a_grids,
                          fineDomain,
                          m_refToFiner,
                          a_nComp /* ncomp*/);
    }
}

// ---------------------------------------------------------
/** full define function for AMRLevelOp<LevelData<FArrayBox> > with finer levels, but no coarser */
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                          const DisjointBoxLayout& a_gridsFiner,
                          const Real&              a_dxLevel,
                          int                      a_refRatio, // dummy arg
                          int                      a_refRatioFiner,
                          const ProblemDomain&     a_domain,
                          BCHolder                 a_bc,
                          const Copier&            a_exchange,
                          const CFRegion&          a_cfregion,
                          const int                a_nComp)
{
  CH_TIME("AMRNonLinearPoissonOp::define2");

  CH_assert(a_refRatio == 1);

  amrpgetMultiColors(m_colors);

  // calls the MG version of define
  this->define(a_grids, a_dxLevel, a_domain, a_bc, a_exchange, a_cfregion);
  m_refToFiner = a_refRatioFiner;

  ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
  m_levfluxreg.define(a_gridsFiner,
                      a_grids,
                      fineDomain,
                      m_refToFiner,
                      a_nComp /* ncomp*/);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                          const DisjointBoxLayout& a_coarse,
                          const Real&              a_dxLevel,
                          int                      a_refRatio,
                          const ProblemDomain&     a_domain,
                          BCHolder                 a_bc,
                          const Copier&            a_exchange,
                          const CFRegion&          a_cfregion,
                          int a_numComp)
{
  CH_TIME("AMRNonLinearPoissonOp::define3");

  amrpgetMultiColors(m_colors);

  this->define(a_grids, a_dxLevel, a_domain, a_bc, a_exchange, a_cfregion);
  m_refToCoarser = a_refRatio;

  m_dxCrse = a_refRatio*a_dxLevel;
  m_refToFiner = 1;

  m_interpWithCoarser.define(a_grids, &a_coarse, a_dxLevel,
                             m_refToCoarser, a_numComp, m_domain);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                          const Real&              a_dx,
                          const ProblemDomain&     a_domain,
                          BCHolder                 a_bc,
                          const Copier&            a_exchange,
                          const CFRegion&          a_cfregion)
{
  CH_TIME("AMRNonLinearPoissonOp::define4");

  amrpgetMultiColors(m_colors);

  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
  m_dxCrse = 2*a_dx;

  // redefined in AMRLevelOp<LevelData<FArrayBox> >::define virtual function.
  m_refToCoarser = 2;
  m_refToFiner   = 2;

  // these get set again after define is called
  m_alpha = 0.0;
  m_beta  = 1.0;

  m_gamma = 0.0;

  m_exchangeCopier = a_exchange;
  // m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  // m_exchangeCopier.trimEdges(a_grids, IntVect::Unit);

  m_cfregion = a_cfregion;

  m_homoCFinterp = true;
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                          const Real&              a_dx,
                          const ProblemDomain&     a_domain,
                          BCHolder                 a_bc)
{
  CH_TIME("AMRNonLinearPoissonOp::define5");

  amrpgetMultiColors(m_colors);

  // 
  Copier copier;
  copier.define(a_grids, a_grids, IntVect::Unit, true);

  CFRegion cfregion(a_grids, a_domain);

  this->define(a_grids, a_dx, a_domain, a_bc, copier, cfregion);
}

void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                          const DisjointBoxLayout* a_baseBAPtr,
                          Real                     a_dx,
                          int                      a_refRatio,
                          const ProblemDomain&     a_domain,
                          BCHolder                 a_bc)
{
  Copier exchangeCopier;
  exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  exchangeCopier.trimEdges(a_grids, IntVect::Unit);

  CFRegion cfregion;
  cfregion.define(a_grids, a_domain);

  if (a_baseBAPtr != NULL)
    {
      define(a_grids, *a_baseBAPtr, a_dx, a_refRatio, a_domain, a_bc, exchangeCopier, cfregion);
    }
  else
    {
      define(a_grids, a_dx, a_domain, a_bc, exchangeCopier, cfregion);
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::residual(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_phi,
                            const LevelData<FArrayBox>& a_rhs,
                            bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::residual");

  if (a_homogeneous)
  {
          homogeneousCFInterp((LevelData<FArrayBox>&)a_phi);
  }
  residualI(a_lhs,a_phi,a_rhs,a_homogeneous);
}

void AMRNonLinearPoissonOp::residualNF(LevelData<FArrayBox>& a_lhs,
                  LevelData<FArrayBox>& a_phi,
                  const LevelData<FArrayBox>* a_phiCoarse,
                  const LevelData<FArrayBox>& a_rhs,
                  bool a_homogeneous)
{
         CH_TIME("AMRNonLinearPoissonOp::residualNF");

          if (a_homogeneous)
          {
                  homogeneousCFInterp((LevelData<FArrayBox>&)a_phi);
          }
          else if (a_phiCoarse != NULL)
          {
                  m_interpWithCoarser.coarseFineInterp(a_phi, *a_phiCoarse);
          }

          residualI(a_lhs,a_phi,a_rhs,a_homogeneous);
}


// ---------------------------------------------------------
void AMRNonLinearPoissonOp::residualI(LevelData<FArrayBox>&       a_lhs,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_rhs,
                             bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::residualI");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  if (s_exchangeMode == 0)
    phi.exchange(phi.interval(), m_exchangeCopier);
  else if (s_exchangeMode == 1)
    phi.exchangeNoOverlap(m_exchangeCopier);
  else
    MayDay::Abort("exchangeMode");

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  {
    CH_TIME("AMRNonLinearPoissonOp::BCs");

    for (dit.begin(); dit.ok(); ++dit)
      {
        m_bc(phi[dit], dbl[dit], m_domain, m_dx, a_homogeneous);
      }
  }

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit];
      FORT_OPERATORLAPRESNL(CHF_FRA(a_lhs[dit]),
                          CHF_CONST_FRA(phi[dit]),
                          CHF_CONST_FRA(a_rhs[dit]),
                          CHF_BOX(region),
                          CHF_CONST_REAL(m_dx),
                          CHF_CONST_REAL(m_alpha),
                          CHF_CONST_REAL(m_beta),
                          CHF_CONST_REAL(m_gamma));
    }
}

// ---------------------------------------------------------
/**************************/
// Preconditionior is not used for FAS solve, but retained so
// regular multigrid can still be used (which should be identical to FAS
// multigrid for linear problems)
void AMRNonLinearPoissonOp::preCond(LevelData<FArrayBox>&       a_phi,
                           const LevelData<FArrayBox>& a_rhs)
{

  CH_TIME("AMRNonLinearPoissonOp::preCond");
//  MayDay::Error("AMRNonLinearPoissonOp::preCond not implemented");

  // diagonal term of this operator is (alpha - 4 * beta/h/h) in 2D,
  // (alpha - 6 * beta/h/h) in 3D,
  // so inverse of this is our initial multiplier

  CH_assert(a_phi.nComp() == a_rhs.nComp());

  Real mult = 1.0 / (m_alpha - 2.0*SpaceDim * m_beta / (m_dx*m_dx));

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phi.dataIterator();
  int nbox = dit.size();

#pragma omp parallel for
    for(int ibox=0; ibox<nbox; ibox++)
      {
      a_phi[dit[ibox]].copy(a_rhs[dit[ibox]]);
      a_phi[dit[ibox]] *= mult;
      }
 //end pragma
  relax(a_phi, a_rhs, 2);
}

void AMRNonLinearPoissonOp::applyOpMg(LevelData<FArrayBox>& a_lhs,
                        LevelData<FArrayBox>& a_phi,
                        LevelData<FArrayBox>* a_phiCoarse,
                        bool a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOpMg");

  // Do CF stuff if we have a coarser level that's not just a single grid cell
   if (a_phiCoarse != NULL)
   {
     const ProblemDomain& probDomain = a_phiCoarse->disjointBoxLayout().physDomain();
     const Box& domBox = probDomain.domainBox();
     //    IntVect hi = domBox.b
     if (domBox.bigEnd() != domBox.smallEnd())
     {
       m_interpWithCoarser.coarseFineInterp(a_phi, *a_phiCoarse);
     }
   }

   applyOpI(a_lhs, a_phi, a_homogeneous);

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::applyOp(LevelData<FArrayBox>&       a_lhs,
                           const LevelData<FArrayBox>& a_phi,
                           bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOp");

  homogeneousCFInterp((LevelData<FArrayBox>&)a_phi);
  applyOpI(a_lhs,a_phi,a_homogeneous);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::applyOpI(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_phi,
                            bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOpI");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  if (s_exchangeMode == 0)
    phi.exchange(phi.interval(), m_exchangeCopier);
  else if (s_exchangeMode == 1)
    phi.exchangeNoOverlap(m_exchangeCopier);
  else
    MayDay::Abort("exchangeMode");

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  int nbox=dit.size();
#pragma omp parallel   default (shared)
  {
    CH_TIME("AMRNonLinearPoissonOp::applyOpIBC");
#pragma omp for 
    for (int ibox=0;ibox<nbox; ibox++)
      {
        m_bc(phi[dit[ibox]], dbl[dit[ibox]], m_domain, m_dx, a_homogeneous);
      }
  }// end pragma

#pragma omp parallel 
  {
#pragma omp for 
    for (int ibox=0;ibox<nbox; ibox++)
      {
      const Box& region = dbl[dit[ibox]];

      FORT_OPERATORLAPNL(CHF_FRA(a_lhs[dit[ibox]]),
                       CHF_CONST_FRA(phi[dit[ibox]]),
                       CHF_BOX(region),
                       CHF_CONST_REAL(m_dx),
                       CHF_CONST_REAL(m_alpha),
                       CHF_CONST_REAL(m_beta),
                       CHF_CONST_REAL(m_gamma));
    }
  }//end pragma
}

void AMRNonLinearPoissonOp::applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                     const LevelData<FArrayBox>& a_phi)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOpNoBoundary");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  int nbox=dit.size();
  phi.exchange(phi.interval(), m_exchangeCopier);

#pragma omp parallel 
  {
#pragma omp for 
    for (int ibox = 0; ibox < nbox; ibox++)
    {
      const Box& region = dbl[dit[ibox]];
      FORT_OPERATORLAPNL(CHF_FRA(a_lhs[dit[ibox]]),
                       CHF_CONST_FRA(phi[dit[ibox]]),
                       CHF_BOX(region),
                       CHF_CONST_REAL(m_dx),
                       CHF_CONST_REAL(m_alpha),
                       CHF_CONST_REAL(m_beta),
                       CHF_CONST_REAL(m_gamma));
    }
  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::create(LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::create");

  m_levelOps.create(a_lhs, a_rhs);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::createCoarsened(LevelData<FArrayBox>&       a_lhs,
                                   const LevelData<FArrayBox>& a_rhs,
                                   const int &                 a_refRat)
{
  CH_TIME("AMRNonLinearPoissonOp::createCoarsened");

  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  const DisjointBoxLayout& dbl = a_rhs.disjointBoxLayout();
  CH_assert(dbl.coarsenable(a_refRat));

  DisjointBoxLayout dblCoarsenedFine;

  if (a_refRat == 2)
    {
      if (m_coarsenedMGrids.size() == 0)
        coarsen(m_coarsenedMGrids, dbl, 2);
      dblCoarsenedFine = m_coarsenedMGrids;
    }
  else
    {
      coarsen(dblCoarsenedFine, dbl, a_refRat);
    }

  a_lhs.define(dblCoarsenedFine, ncomp, ghostVect );
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::assign(LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::assign");

  m_levelOps.assign(a_lhs, a_rhs);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::assignLocal(LevelData<FArrayBox>&        a_lhs,
                               const  LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::assignLocal");

  for (DataIterator dit= a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      a_lhs[dit].copy(a_rhs[dit]);
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::buildCopier(Copier&                      a_copier,
                               const  LevelData<FArrayBox>& a_lhs,
                               const  LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::buildCopier");

  const DisjointBoxLayout& dbl=a_lhs.disjointBoxLayout();

  a_copier.define(a_rhs.disjointBoxLayout(), dbl, IntVect::Zero);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::assignCopier( LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_rhs,
                                 const Copier&               a_copier)
{
  CH_TIME("AMRNonLinearPoissonOp::assignCopier");
  a_rhs.copyTo(a_rhs.interval(), a_lhs, a_lhs.interval(), a_copier);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::zeroCovered( LevelData<FArrayBox>& a_lhs,
                                LevelData<FArrayBox>& a_rhs,
                                const Copier&         a_copier)
{
  CH_TIME("AMRNonLinearPoissonOp::zeroCovered");

  m_levelOps.copyToZero(a_lhs, a_copier);
}

// ---------------------------------------------------------
Real AMRNonLinearPoissonOp::dotProduct(const LevelData<FArrayBox>& a_1,
                              const LevelData<FArrayBox>& a_2)
{
  CH_TIME("AMRNonLinearPoissonOp::dotProduct");

  return m_levelOps.dotProduct(a_1, a_2);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::mDotProduct(const LevelData<FArrayBox>& a_1,
                               const int a_sz,
                               const LevelData<FArrayBox> a_2[],
                               Real a_mdots[])
{
  CH_TIME("AMRNonLinearPoissonOp::mDotProduct");

  m_levelOps.mDotProduct(a_1, a_sz, a_2, a_mdots);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::incr( LevelData<FArrayBox>&       a_lhs,
                         const LevelData<FArrayBox>& a_x,
                         Real                        a_scale)
{
  CH_TIME("AMRNonLinearPoissonOp::incr");

  m_levelOps.incr(a_lhs, a_x, a_scale);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::axby( LevelData<FArrayBox>&       a_lhs,
                         const LevelData<FArrayBox>& a_x,
                         const LevelData<FArrayBox>& a_y,
                         Real                        a_a,
                         Real                        a_b)
{
  CH_TIME("AMRNonLinearPoissonOp::axby");

  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::scale(LevelData<FArrayBox>& a_lhs,
                         const Real&           a_scale)
{
  CH_TIME("AMRNonLinearPoissonOp::scale");

  m_levelOps.scale(a_lhs, a_scale);
}

// ---------------------------------------------------------
Real AMRNonLinearPoissonOp::norm(const LevelData<FArrayBox>& a_x,
                        int                         a_ord)
{
  CH_TIME("AMRNonLinearPoissonOp::norm");

  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}

// ---------------------------------------------------------
Real AMRNonLinearPoissonOp::localMaxNorm(const LevelData<FArrayBox>& a_x)
{
  CH_TIME("AMRNonLinearPoissonOp::localMaxNorm");

  Real localMax = 0;
  int nComp=a_x.nComp();
  for (DataIterator dit=a_x.dataIterator(); dit.ok(); ++dit)
    {
      localMax = Max(localMax, a_x[dit].norm(a_x.box(dit()), 0, 0, nComp));
    }
  return localMax;
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::setToZero(LevelData<FArrayBox>& a_lhs)
{
  CH_TIME("AMRNonLinearPoissonOp::setToZero");

  m_levelOps.setToZero(a_lhs);
}

void AMRNonLinearPoissonOp::relaxNF(LevelData<FArrayBox>&       a_e,
                         const LevelData<FArrayBox>* a_eCoarse,
                         const LevelData<FArrayBox>& a_residual,
                         int                         a_iterations)
{
  if (a_eCoarse != NULL)
  {
    m_interpWithCoarser.coarseFineInterp(a_e, *a_eCoarse);
  }
  else
  {
    homogeneousCFInterp(a_e);
  }

  relax(a_e, a_residual, a_iterations);

}
// ---------------------------------------------------------
// Only option 6 - Gauss Seidel - is currently implemented
void AMRNonLinearPoissonOp::relax(LevelData<FArrayBox>&       a_e,
                         const LevelData<FArrayBox>& a_residual,
                         int                         a_iterations)
{
  CH_TIME("AMRNonLinearPoissonOp::relax");

  for (int i = 0; i < a_iterations; i++)
    {
      switch (s_relaxMode)
        {
        case 0:
          looseGSRB(a_e, a_residual);
          break;
        case 1:
          levelGSRB(a_e, a_residual);
          break;
        case 2:
          overlapGSRB(a_e, a_residual);
          break;
        case 3:
          levelGSRBLazy(a_e, a_residual);
          break;
        case 4:
          levelJacobi(a_e, a_residual);
          break;
        case 5:
          levelMultiColor(a_e, a_residual);
          break;
        case 6:
          levelGS(a_e, a_residual);
          break;
        default:
          MayDay::Abort("unrecognized relaxation mode");
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::createCoarser(LevelData<FArrayBox>&       a_coarse,
                                 const LevelData<FArrayBox>& a_fine,
                                 bool                        a_ghosted)
{
  CH_TIME("AMRNonLinearPoissonOp::createCoarser");

  // CH_assert(!a_ghosted);
  IntVect ghost = a_fine.ghostVect();

  CH_assert(a_fine.disjointBoxLayout().coarsenable(2));
  if (m_coarsenedMGrids.size() == 0)
    coarsen(m_coarsenedMGrids, a_fine.disjointBoxLayout(), 2); //multigrid, so coarsen by 2
  a_coarse.define(m_coarsenedMGrids, a_fine.nComp(), ghost);
}


void AMRNonLinearPoissonOp::restrictR(LevelData<FArrayBox>& a_phiCoarse,
                                      const LevelData<FArrayBox>& a_phiFine)
{
  //    a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox&       phiFine = a_phiFine[dit];
    FArrayBox&       phiCoarse = a_phiCoarse[dit];

    Box region = dblFine.get(dit());
    const IntVect& iv = region.smallEnd();
    IntVect civ = coarsen(iv, 2);

    phiCoarse.setVal(0.0);

    FORT_RESTRICT(CHF_FRA_SHIFT(phiCoarse, civ),
                  CHF_CONST_FRA_SHIFT(phiFine, iv),
                  CHF_BOX_SHIFT(region, iv),
                  CHF_CONST_REAL(m_dx));
  }
}


// ---------------------------------------------------------
void AMRNonLinearPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                    LevelData<FArrayBox>&       a_phiFine,
                                    const LevelData<FArrayBox>& a_rhsFine)
{
  CH_TIME("AMRNonLinearPoissonOp::restrictResidual");

  homogeneousCFInterp(a_phiFine);

  if (s_exchangeMode == 0)
    a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);
  else if (s_exchangeMode == 1)
    a_phiFine.exchangeNoOverlap(m_exchangeCopier);
  else
    MayDay::Abort("exchangeMode");

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  DataIterator dit = a_phiFine.dataIterator();
  int nbox=dit.size();
#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi = a_phiFine[dit[ibox]];
        m_bc(phi, dblFine[dit[ibox]], m_domain, m_dx, true);
      }
  }//end pragma

#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox&       phi = a_phiFine[dit[ibox]];
        const FArrayBox& rhs = a_rhsFine[dit[ibox]];
        FArrayBox&       res = a_resCoarse[dit[ibox]];
        
        Box region = dblFine[dit[ibox]];
        const IntVect& iv = region.smallEnd();
        IntVect civ = coarsen(iv, 2);
        
        res.setVal(0.0);
        
        FORT_RESTRICTRESNL(CHF_FRA_SHIFT(res, civ),
                           CHF_CONST_FRA_SHIFT(phi, iv),
                           CHF_CONST_FRA_SHIFT(rhs, iv),
                           CHF_CONST_REAL(m_alpha),
                           CHF_CONST_REAL(m_beta),
                           CHF_CONST_REAL(m_gamma),
                           CHF_BOX_SHIFT(region, iv),
                           CHF_CONST_REAL(m_dx));
      }
  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                    const LevelData<FArrayBox>& a_correctCoarse)
{
  CH_TIME("AMRNonLinearPoissonOp::prolongIncrement");
  
  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  int mgref = 2; //this is a multigrid func
  DataIterator dit = a_phiThisLevel.dataIterator();
  int nbox=dit.size();
  
#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_phiThisLevel[dit[ibox]];
        const FArrayBox& coarse = a_correctCoarse[dit[ibox]];
        Box region = dbl[dit[ibox]];
        const IntVect& iv = region.smallEnd();
        IntVect civ=coarsen(iv, 2);
        
        FORT_PROLONGNL(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(mgref));
      }
  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRResidual(LevelData<FArrayBox>&              a_residual,
                               const LevelData<FArrayBox>&        a_phiFine,
                               const LevelData<FArrayBox>&        a_phi,
                               const LevelData<FArrayBox>&        a_phiCoarse,
                               const LevelData<FArrayBox>&        a_rhs,
                               bool                               a_homogeneousPhysBC,
                               AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRResidual");

  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse,
              a_homogeneousPhysBC, a_finerOp);

  axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRResidualNC(LevelData<FArrayBox>&              a_residual,
                                 const LevelData<FArrayBox>&        a_phiFine,
                                 const LevelData<FArrayBox>&        a_phi,
                                 const LevelData<FArrayBox>&        a_rhs,
                                 bool                               a_homogeneousPhysBC,
                                 AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRResidualNC");

   AMROperatorNC(a_residual, a_phiFine, a_phi,
                 a_homogeneousPhysBC, a_finerOp);

   axby(a_residual, a_residual, a_rhs, -1.0, 1.0);

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRResidualNF(LevelData<FArrayBox>&       a_residual,
                                 const LevelData<FArrayBox>& a_phi,
                                 const LevelData<FArrayBox>& a_phiCoarse,
                                 const LevelData<FArrayBox>& a_rhs,
                                 bool                        a_homogeneousPhysBC)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRResidualNF");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined())
  {
    m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  }

  //apply boundary conditions
  this->residualI(a_residual, a_phi, a_rhs, a_homogeneousPhysBC);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMROperator(LevelData<FArrayBox>&              a_LofPhi,
                               const LevelData<FArrayBox>&        a_phiFine,
                               const LevelData<FArrayBox>&        a_phi,
                               const LevelData<FArrayBox>&        a_phiCoarse,
                               bool                               a_homogeneousPhysBC,
                               AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMROperator");

  CH_assert(a_phi.isDefined());

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined())
  {
    m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  }

  // apply physical boundary conditions in applyOpI
  applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);

  if (a_phiFine.isDefined())
  {
    CH_assert(a_finerOp != NULL);
    reflux(a_phiFine, a_phi, a_LofPhi, a_finerOp);
  }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMROperatorNC(LevelData<FArrayBox>&              a_LofPhi,
                                 const LevelData<FArrayBox>&        a_phiFine,
                                 const LevelData<FArrayBox>&        a_phi,
                                 bool                               a_homogeneousPhysBC,
                                 AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMROperatorNC");

  CH_assert(a_phi.isDefined());

  // apply physical boundary conditions in applyOpI
  applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);

//  Real sum = ::computeSum(a_LofPhi, NULL, -1, m_dx);
//  pout() << "LofPhi sum before reflux = " << sum << endl;

  if (a_phiFine.isDefined())
  {
    CH_assert(a_finerOp != NULL);
    reflux(a_phiFine, a_phi, a_LofPhi, a_finerOp);

//    sum = ::computeSum(a_LofPhi, NULL, -1, m_dx);
//    pout() << "LofPhi sum after reflux = " << sum << endl;
  }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                                 const LevelData<FArrayBox>& a_phi,
                                 const LevelData<FArrayBox>& a_phiCoarse,
                                 bool                        a_homogeneousPhysBC)
{
  CH_TIME("AMRNonLinearPoissonOp::AMROperatorNF");

  CH_assert(a_phi.isDefined());

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined())
  {
    m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  }

  // apply physical boundary conditions in applyOpI
  this->applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
                                        const LevelData<FArrayBox>& a_residual,
                                        const LevelData<FArrayBox>& a_correction,
                                        const LevelData<FArrayBox>& a_coarseCorrection,
                                        bool a_skip_res )
{
  CH_TIME("AMRNonLinearPoissonOp::AMRRestrict");

  LevelData<FArrayBox> r;
  create(r, a_residual);

  AMRRestrictS(a_resCoarse, a_residual, a_correction, a_coarseCorrection, r, a_skip_res );
}

// hacks to do BCs right
// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRRestrictS(LevelData<FArrayBox>&       a_resCoarse, // output
                                         const LevelData<FArrayBox>& a_residual,  // input
                                         const LevelData<FArrayBox>& a_correction, // for residual
                                         const LevelData<FArrayBox>& a_coarseCorrection, // C-F interp.
                                         LevelData<FArrayBox>&       a_scratch,   // temp buffer
                                         bool a_skip_res // flag skipping computing residual, used for lhs restrictions (FAS)
                                         ) 
{
  CH_TIME("AMRNonLinearPoissonOp::AMRRestrictS");

  if ( !a_skip_res )
    {
      AMRResidualNF( a_scratch, a_correction, a_coarseCorrection, a_residual, true );
    }
  else 
    {
      // just copy data (phi in this case, even if its called residual)
      assignLocal( a_scratch, a_residual );
    }

  DisjointBoxLayout dblCoar = a_resCoarse.disjointBoxLayout();

#pragma omp parallel 
  {
    DataIterator dit = a_residual.dataIterator();
    int nbox=dit.size();
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& coarse = a_resCoarse[dit[ibox]];
        const FArrayBox& fine = a_scratch[dit[ibox]];
        const Box& b = dblCoar[dit[ibox]];
        Box refbox(IntVect::Zero,
                   (m_refToCoarser-1)*IntVect::Unit);
        FORT_AVERAGE( CHF_FRA(coarse),
                      CHF_CONST_FRA(fine),
                      CHF_BOX(b),
                      CHF_CONST_INT(m_refToCoarser),
                      CHF_BOX(refbox)
                      );
      }
  }//end pragma
}

// ---------------------------------------------------------
/** a_correction += I[2h->h](a_coarseCorrection) */
void AMRNonLinearPoissonOp::AMRProlong(LevelData<FArrayBox>&       a_correction,
                              const LevelData<FArrayBox>& a_coarseCorrection)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRProlong");

  DisjointBoxLayout c;
  coarsen(c, a_correction.disjointBoxLayout(), m_refToCoarser);

  LevelData<FArrayBox> eCoar(c, a_correction.nComp(), a_coarseCorrection.ghostVect());
  a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  DataIterator dit = a_correction.dataIterator();
  int nbox=dit.size();

#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_correction[dit[ibox]];
        const FArrayBox& coarse = eCoar[dit[ibox]];
        
        Box region = dbl[dit[ibox]];
        const IntVect& iv = region.smallEnd();
        IntVect civ = coarsen(iv, m_refToCoarser);
        
        FORT_PROLONGNL(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(m_refToCoarser));
      }
  }// end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRProlongS(LevelData<FArrayBox>&       a_correction,
                               const LevelData<FArrayBox>& a_coarseCorrection,
                               LevelData<FArrayBox>&       a_temp,
                               const Copier&               a_copier)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRProlongS");

  a_coarseCorrection.copyTo(a_temp.interval(), a_temp, a_temp.interval(), a_copier);

  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  DataIterator dit = a_correction.dataIterator();
  int nbox=dit.size();
  
#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_correction[dit[ibox]];
        const FArrayBox& coarse = a_temp[dit[ibox]];
        
        Box region = dbl[dit[ibox]];
        const IntVect& iv =  region.smallEnd();
        IntVect civ= coarsen(iv, m_refToCoarser);
        
        FORT_PROLONGNL(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(m_refToCoarser));
      }
  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRProlongS_2(LevelData<FArrayBox>&       a_correction,
                                 const LevelData<FArrayBox>& a_coarseCorrection,
                                 LevelData<FArrayBox>&       a_temp,
                                 const Copier&               a_copier,
                                 const Copier&         a_cornerCopier,
                                 const AMRLevelOp<LevelData<FArrayBox> >* a_crsOp )
{
  CH_TIME("AMRNonLinearPoissonOp::AMRProlongS_2");
  
  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  DisjointBoxLayout cdbl = a_temp.disjointBoxLayout();
  AMRNonLinearPoissonOp* coarserAMRPOp = (AMRNonLinearPoissonOp*) a_crsOp;
  
  a_coarseCorrection.copyTo( a_temp.interval(), a_temp, a_temp.interval(), a_copier );

  // I think we should be using a coarse data iterator for applying the coarse BC?
  for (DataIterator cdit = a_temp.dataIterator(); cdit.ok(); ++cdit)
  {
    FArrayBox& coarse = a_temp[cdit];
    coarserAMRPOp->m_bc( coarse, cdbl[cdit], coarserAMRPOp->m_domain, coarserAMRPOp->m_dx, true );
  }

  // The corner copier passed in as an argument doesn't always work, whilst this one seems better
  CornerCopier cornerCopy(cdbl, cdbl, coarserAMRPOp->m_domain, a_temp.ghostVect(), true);

  a_temp.exchange(cornerCopy);
  //a_temp.exchange( a_temp.interval(), a_cornerCopier ); -- needed for AMR
#pragma omp parallel
  {
    DataIterator dit = a_correction.dataIterator();
    int nbox = dit.size();
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_correction[dit[ibox]];
        FArrayBox& coarse = a_temp[dit[ibox]];
        
        Box region = dbl[dit[ibox]];
        const IntVect& iv = region.smallEnd();
        IntVect civ = coarsen(iv, m_refToCoarser);
        
#if 0
        FORT_PROLONGNL( CHF_FRA_SHIFT(phi, iv),
                      CHF_CONST_FRA_SHIFT(coarse, civ),
                      CHF_BOX_SHIFT(region, iv),
                      CHF_CONST_INT(m_refToCoarser));
#else
        FORT_PROLONG_2_NL( CHF_FRA_SHIFT(phi, iv),
                        CHF_CONST_FRA_SHIFT(coarse, civ),
                        CHF_BOX_SHIFT(region, iv),
                        CHF_CONST_INT(m_refToCoarser) 
                        );
#endif
      }
  }//end pragma

  // debug
  //write(&a_temp,"z_src.hdf5"); 
  //write(&a_correction,"z_prol2.hdf5"); exit(12);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
                                     const LevelData<FArrayBox>& a_correction,
                                     const LevelData<FArrayBox>& a_coarseCorrection)
{
//   LevelData<FArrayBox> r;
//   this->create(r, a_residual);
//   this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
//   this->assign(a_residual, r);
  this->AMRResidualNF(a_residual, a_correction, a_coarseCorrection, a_residual, true);
}

// ---------------------------------------------------------
///compute norm over all cells on coarse not covered by finer
Real AMRNonLinearPoissonOp::AMRNorm(const LevelData<FArrayBox>& a_coarResid,
                           const LevelData<FArrayBox>& a_fineResid,
                           const int&                  a_refRat,
                           const int&                  a_ord)

{
  CH_TIME("AMRNonLinearPoissonOp::AMRNorm");

  //create temp and zero out under finer grids
  LevelData<FArrayBox> coarTemp;
  m_levelOps.create(coarTemp, a_coarResid);
  m_levelOps.assign(coarTemp, a_coarResid);

  if (a_fineResid.isDefined())
  {
    const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();
    const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();

    int ncomp = coarTemp.nComp();

    for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& coarTempFAB = coarTemp[dit];
        LayoutIterator litFine = fineGrids.layoutIterator();

        for (litFine.reset(); litFine.ok(); ++litFine)
          {
            Box overlayBox = coarTempFAB.box();
            Box coarsenedGrid = coarsen(fineGrids[litFine], a_refRat);

            overlayBox &= coarsenedGrid;

            if (!overlayBox.isEmpty())
              {
                coarTempFAB.setVal(0.0, overlayBox, 0, ncomp);
              }
          }
      }
  }

  // return norm of temp
  return norm(coarTemp, a_ord);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::setAlphaAndBeta(const Real& a_alpha,
                                   const Real& a_beta)
{
  m_alpha = a_alpha * m_aCoef;
  m_beta  = a_beta  * m_bCoef;
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::setBC(const BCHolder& a_bc)
{
  m_bc = a_bc;
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                          const LevelData<FArrayBox>&        a_phi,
                          LevelData<FArrayBox>&              a_residual,
                          AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIMERS("AMRNonLinearPoissonOp::reflux");

  m_levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  CH_TIMER("AMRNonLinearPoissonOp::reflux::incrementCoarse", t2);
  CH_START(t2);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab = a_phi[dit];

      if (m_levfluxreg.hasCF(dit()))
        {
          //pout() << " AMRNonLinearPoissonOp::reflux: " << coarfab.box() << endl;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox coarflux;
              getFlux(coarflux, coarfab, idir);

              Real scale = 1.0;
              m_levfluxreg.incrementCoarse(coarflux, scale, dit(),
                                           interv, interv, idir);
            }
        }
    }

  CH_STOP(t2);

  // const cast:  OK because we're changing ghost cells only
  LevelData<FArrayBox>& phiFineRef = ( LevelData<FArrayBox>&)a_phiFine;

  AMRNonLinearPoissonOp* finerAMRPOp = (AMRNonLinearPoissonOp*) a_finerOp;
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(phiFineRef, a_phi);
  // I'm pretty sure this is not necessary. bvs -- flux calculations use
  // outer ghost cells, but not inner ones
  // phiFineRef.exchange(a_phiFine.interval());
  IntVect phiGhost = phiFineRef.ghostVect();
  int ncomps = a_phiFine.nComp();

  CH_TIMER("AMRNonLinearPoissonOp::reflux::incrementFine", t3);
  CH_START(t3);

  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = a_phiFine[ditf];
      const Box& gridbox = dblFine[ditf];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //int normalGhost = phiGhost[idir];
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              if (m_levfluxreg.hasCF(ditf(), sit()))
                {
                  Side::LoHiSide hiorlo = sit();
                  Box fluxBox = bdryBox(gridbox,idir,hiorlo,1);

                  FArrayBox fineflux(fluxBox,ncomps);
                  getFlux(fineflux, phifFab, fluxBox, idir, m_refToFiner);

                  Real scale = 1.0;
                  m_levfluxreg.incrementFine(fineflux, scale, ditf(),
                                             interv, interv, idir, hiorlo);
                }
            }
        }
    }

  CH_STOP(t3);

  Real scale = 1.0/m_dx;
  m_levfluxreg.reflux(a_residual, scale);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::write(const LevelData<FArrayBox>* a_data,
                         const char*                 a_filename)
{
#ifdef CH_USE_HDF5
  writeLevelname(a_data, a_filename);
#else
  MayDay::Warning("AMRNonLinearPoissonOp::write unimplemented since CH_USE_HDF5 undefined");
#endif
}

/***/
// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelGSRB( LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs )
{
  CH_TIME("AMRNonLinearPoissonOp::levelGSRB");
//  MayDay::Error("AMRNonLinearPoissonOp::levelGSRB not implemented");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  const DisjointBoxLayout& dbl = a_rhs.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();
  int nbox=dit.size();
  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {
      CH_TIME("AMRNonLinearPoissonOp::levelGSRB::Compute");

      // fill in intersection of ghostcells and a_phi's boxes
      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::homogeneousCFInterp");
        // For FAS, we don't want to do this?
//        homogeneousCFInterp(a_phi);
      }

      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::exchange");
        if (s_exchangeMode == 0)
          a_phi.exchange( a_phi.interval(), m_exchangeCopier );
        else if (s_exchangeMode == 1)
          a_phi.exchangeNoOverlap(m_exchangeCopier);
        else
          MayDay::Abort("exchangeMode");
      }
#pragma omp parallel
      {
#pragma omp for
        for (int ibox=0; ibox < nbox; ibox++)
          {
            const Box& region = dbl[dit[ibox]];
            FArrayBox& phiFab = a_phi[dit[ibox]];

            m_bc( phiFab, region, m_domain, m_dx, true );

            if (m_alpha == 0.0 && m_beta == 1.0 )
              {
                FORT_GSRBLAPLACIANNL(CHF_FRA(phiFab),
                                   CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                   CHF_BOX(region),
                                   CHF_CONST_REAL(m_dx),
                                   CHF_CONST_INT(whichPass),
                                   CHF_CONST_REAL(m_gamma));
              }
            else
              {
              MayDay::Error("AMRNonLinearPoissonOp: FORTRAN routine not tested");
                FORT_GSRBHELMHOLTZNL(CHF_FRA(phiFab),
                                   CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                   CHF_BOX(region),
                                   CHF_CONST_REAL(m_dx),
                                   CHF_CONST_REAL(m_alpha),
                                   CHF_CONST_REAL(m_beta),
                                   CHF_CONST_REAL(m_gamma),
                                   CHF_CONST_INT(whichPass));
              }
          } // end loop through grids
      }//end pragma
    } // end loop through red-black
}

void AMRNonLinearPoissonOp::levelGS( LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs )
{

  CH_TIME("AMRNonLinearPoissonOp::levelGS");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  const DisjointBoxLayout& dbl = a_rhs.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();
  int nbox=dit.size();
  // do first red, then black passes
  int whichPass = 0;

      //CH_TIME("AMRNonLinearPoissonOp::levelGSRB::Compute");

      // fill in intersection of ghostcells and a_phi's boxes
      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::homogeneousCFInterp");
//        homogeneousCFInterp(a_phi);
      }

      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::exchange");
        if (s_exchangeMode == 0)
          a_phi.exchange( a_phi.interval(), m_exchangeCopier );
        else if (s_exchangeMode == 1)
          a_phi.exchangeNoOverlap(m_exchangeCopier);
        else
          MayDay::Abort("exchangeMode");
      }
#pragma omp parallel
      {
#pragma omp for
        for (int ibox=0; ibox < nbox; ibox++)
          {
            const Box& region = dbl[dit[ibox]];
            FArrayBox& phiFab = a_phi[dit[ibox]];

            m_bc( phiFab, region, m_domain, m_dx, true );

            if (m_alpha == 0.0 && m_beta == 1.0 )
              {
                FORT_GSLAPLACIANNL(CHF_FRA(phiFab),
                                   CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                   CHF_BOX(region),
                                   CHF_CONST_REAL(m_dx),
                                   CHF_CONST_INT(whichPass),
                                   CHF_CONST_REAL(m_gamma));
              }
            else
              {
              MayDay::Error("AMRNonLinearPoissonOp: FORT_GSHELMHOLTZNL not implemented");
               /* FORT_GSHELMHOLTZNL(CHF_FRA(phiFab),
                                   CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                   CHF_BOX(region),
                                   CHF_CONST_REAL(m_dx),
                                   CHF_CONST_REAL(m_alpha),
                                   CHF_CONST_REAL(m_beta),
                                   CHF_CONST_REAL(m_gamma),
                                   CHF_CONST_INT(whichPass));*/
              }
          } // end loop through grids
      }//end pragma

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelMultiColor(LevelData<FArrayBox>&       a_phi,
                                   const LevelData<FArrayBox>& a_rhs)
{
  MayDay::Error("AMRNonLinearPoissonOp::levelMultiColor not implemented");

//  CH_assert(a_phi.isDefined());
//  CH_assert(a_rhs.isDefined());
//  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
//  CH_assert(a_phi.nComp() == a_rhs.nComp());
//  // do first red, then black passes
//  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
//  for (int icolor = 0; icolor < m_colors.size(); icolor++)
//    {
//      const IntVect& color= m_colors[icolor];
//      homogeneousCFInterp(a_phi);
//
//      a_phi.exchange(a_phi.interval(), m_exchangeCopier);
//      //after this lphi = L(phi)
//      //this call contains bcs and exchange
//
//      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
//        {
//          Box dblBox  = dbl[dit];
//
//          IntVect loIV = dblBox.smallEnd();
//          IntVect hiIV = dblBox.bigEnd();
//
//          for (int idir = 0; idir < SpaceDim; idir++)
//            {
//              if (loIV[idir] % 2 != color[idir])
//                {
//                  loIV[idir]++;
//                }
//            }
//
//          m_bc(a_phi[dit], dbl[dit], m_domain, m_dx, true);
//
//          if (loIV <= hiIV)
//            {
//              Box coloredBox(loIV, hiIV);
//              MayDay::Error("AMRNonLinearPoissonOp: FORT_GSMCAMRPOPNL not tested");
//              FORT_GSMCAMRPOPNL(CHF_FRA1(a_phi[dit], 0),
//                              CHF_CONST_FRA1(a_rhs[dit], 0),
//                              CHF_BOX(coloredBox),
//                              CHF_REAL(m_dx),
//                              CHF_REAL(m_alpha),
//                              CHF_REAL(m_beta),
//                              CHF_CONST_REAL(m_gamma) );
//            }
//        }
//    } // end loop through red-black
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::looseGSRB(LevelData<FArrayBox>&       a_phi,
                             const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::looseGSRB");

  MayDay::Error("AMRNonLinearPoissonOp::looseGSRB not implemented");

//  CH_assert(a_phi.isDefined());
//  CH_assert(a_rhs.isDefined());
//  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
//  CH_assert(a_phi.nComp() == a_rhs.nComp());
//
//  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
//
//  DataIterator dit = a_phi.dataIterator();
//  int nbox=dit.size();
//  //fill in intersection of ghostcells and a_phi's boxes
//  {
//    CH_TIME("AMRNonLinearPoissonOp::looseGSRB::homogeneousCFInterp");
//    homogeneousCFInterp(a_phi);
//  }
//
//  {
//    CH_TIME("AMRNonLinearPoissonOp::looseGSRB::exchange");
//    if (s_exchangeMode == 0)
//      a_phi.exchange(a_phi.interval(), m_exchangeCopier);
//    else if (s_exchangeMode == 1)
//      a_phi.exchangeNoOverlap(m_exchangeCopier);
//    else
//      MayDay::Abort("exchangeMode");
//  }
//
//  // now step through grids...
//#pragma omp parallel
//  {
//#pragma omp for
//    for(int ibox = 0; ibox < nbox; ibox++)
//      {
//      // invoke physical BC's where necessary
//      {
//
//        m_bc(a_phi[dit[ibox]], dbl[dit[ibox]], m_domain, m_dx, true);
//      }
//
//      const Box& region = dbl[dit[ibox]];
//
//      if (m_alpha == 0.0 && m_beta == 1.0)
//        {
//        MayDay::Error("AMRNonLinearPoissonOp: FORTRAN routine not tested");
//          int whichPass = 0;
//          FORT_GSRBLAPLACIANNL(CHF_FRA(a_phi[dit[ibox]]),
//                             CHF_CONST_FRA(a_rhs[dit[ibox]]),
//                             CHF_BOX(region),
//                             CHF_CONST_REAL(m_dx),
//                             CHF_CONST_INT(whichPass),
//                                CHF_CONST_REAL(m_gamma));
//
//          whichPass = 1;
//          FORT_GSRBLAPLACIANNL(CHF_FRA(a_phi[dit[ibox]]),
//                             CHF_CONST_FRA(a_rhs[dit[ibox]]),
//                             CHF_BOX(region),
//                             CHF_CONST_REAL(m_dx),
//                             CHF_CONST_INT(whichPass),
//                                CHF_CONST_REAL(m_gamma));
//        }
//      else
//        {
//          int whichPass = 0;
//          MayDay::Error("AMRNonLinearPoissonOp: FORTRAN routine not tested");
//          FORT_GSRBHELMHOLTZNL(CHF_FRA(a_phi[dit[ibox]]),
//                             CHF_CONST_FRA(a_rhs[dit[ibox]]),
//                             CHF_BOX(region),
//                             CHF_CONST_REAL(m_dx),
//                             CHF_CONST_REAL(m_alpha),
//                             CHF_CONST_REAL(m_beta),
//                                CHF_CONST_REAL(m_gamma),
//                             CHF_CONST_INT(whichPass));
//
//          whichPass = 1;
//          FORT_GSRBHELMHOLTZNL(CHF_FRA(a_phi[dit[ibox]]),
//                             CHF_CONST_FRA(a_rhs[dit[ibox]]),
//                             CHF_BOX(region),
//                             CHF_CONST_REAL(m_dx),
//                             CHF_CONST_REAL(m_alpha),
//                             CHF_CONST_REAL(m_beta),
//                                CHF_CONST_REAL(m_gamma),
//                             CHF_CONST_INT(whichPass));
//        }
//      } // end loop through grids
//  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::overlapGSRB(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::overlapGSRB");
  MayDay::Error("AMRNonLinearPoissonOp::overlapGSRB not implemented");
//
//  CH_assert(a_phi.isDefined());
//  CH_assert(a_rhs.isDefined());
//  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
//  CH_assert(a_phi.nComp() == a_rhs.nComp());
//  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
//
//  if (dbl.size() == 1 || m_alpha != 1 || m_beta != 0)
//    looseGSRB(a_phi, a_rhs);
//
//  a_phi.exchangeBegin(m_exchangeCopier);
//
//  homogeneousCFInterp(a_phi);
//
//  // now step through grids...
//  DataIterator dit = a_phi.dataIterator();
//  int nbox=dit.size();
//#pragma omp parallel
//  {
//#pragma omp for
//    for(int ibox = 0; ibox < nbox; ibox++)
//      {
//      MayDay::Error("AMRNonLinearPoissonOp: FORTRAN routine not tested");
//
//      // invoke physical BC's where necessary
//      m_bc(a_phi[dit[ibox]], dbl[dit[ibox]], m_domain, m_dx, true);
//      Box region = dbl[dit[ibox]];
//      region.grow(-1); // just do the interior on the first run through
//      int whichPass = 0;
//      FORT_GSRBLAPLACIANNL(CHF_FRA(a_phi[dit[ibox]]),
//                         CHF_CONST_FRA(a_rhs[dit[ibox]]),
//                         CHF_BOX(region),
//                         CHF_CONST_REAL(m_dx),
//                         CHF_CONST_INT(whichPass),
//                                CHF_CONST_REAL(m_gamma));
//
//      whichPass = 1;
//      FORT_GSRBLAPLACIANNL(CHF_FRA(a_phi[dit[ibox]]),
//                         CHF_CONST_FRA(a_rhs[dit[ibox]]),
//                         CHF_BOX(region),
//                         CHF_CONST_REAL(m_dx),
//                         CHF_CONST_INT(whichPass),
//                                CHF_CONST_REAL(m_gamma));
//      }
//  }//end pragma
//
//  a_phi.exchangeEnd();
//
//  // now step through grid edges
//
//  for (dit.begin(); dit.ok(); ++dit)
//    {
//
//      const Box& b = dbl[dit];
//      for (int dir=0; dir<CH_SPACEDIM; ++dir)
//      {
//        for (SideIterator side; side.ok(); ++side)
//        {
//          MayDay::Error("AMRNonLinearPoissonOp: FORTRAN routine not tested");
//
//          int whichPass = 0;
//          Box region = adjCellBox(b, dir, side(), 1);
//          region.shift(dir, -sign(side()));
//          for (int i=0; i<dir; i++) region.grow(i, -1);
//          FORT_GSRBLAPLACIANNL(CHF_FRA(a_phi[dit]),
//                             CHF_CONST_FRA(a_rhs[dit]),
//                             CHF_BOX(region),
//                             CHF_CONST_REAL(m_dx),
//                             CHF_CONST_INT(whichPass),
//                             CHF_CONST_REAL(m_gamma));
//
//          whichPass = 1;
//          FORT_GSRBLAPLACIANNL(CHF_FRA(a_phi[dit]),
//                             CHF_CONST_FRA(a_rhs[dit]),
//                             CHF_BOX(region),
//                             CHF_CONST_REAL(m_dx),
//                             CHF_CONST_INT(whichPass),
//                             CHF_CONST_REAL(m_gamma));
//        }
//      }
//    }
}

/***/
// ---------------------------------------------------------
//static bool nextColorLoc(IntVect&       a_color,
//                         const IntVect& a_limit)
//{
//  a_color[0]++;
//
//  for (int i=0; i<CH_SPACEDIM-1; ++i)
//    {
//      if (a_color[i] > a_limit[i])
//        {
//          a_color[i] = 0;
//          a_color[i+1]++;
//        }
//    }
//  if (a_color[CH_SPACEDIM-1] > a_limit[CH_SPACEDIM-1])
//    {
//      return false;
//    }
//
//  return true;
//}

/***/
// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelGSRBLazy( LevelData<FArrayBox>&       a_phi,
                                  const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::levelGSRBLazy");
  MayDay::Error("AMRNonLinearPoissonOp::levelGSRBLazy not implemented");

//  CH_assert(a_phi.isDefined());
//  CH_assert(a_rhs.isDefined());
//  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
//  CH_assert(a_phi.nComp() == a_rhs.nComp());
//  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
//
//  IntVect color = IntVect::Zero;
//  IntVect limit = IntVect::Unit;
//  color[0]=-1;
//
//  // Loop over all possibilities (in all dimensions)
//  while (nextColorLoc(color, limit))
//    {
//      LevelData<FArrayBox> lphi;
//      m_levelOps.create(lphi, a_rhs);
//      //after this lphi = L(phi)
//      //this call contains bcs and exchange
//      applyOp(lphi, a_phi, true);
//
//      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
//        {
//          Box dblBox  = dbl[dit];
//
//          IntVect loIV = dblBox.smallEnd();
//          IntVect hiIV = dblBox.bigEnd();
//
//          for (int idir = 0; idir < SpaceDim; idir++)
//            {
//              if (loIV[idir] % 2 != color[idir])
//                {
//                  loIV[idir]++;
//                }
//            }
//
//          if (loIV <= hiIV)
//            {
//            MayDay::Error("AMRNonLinearPoissonOp: FORTRAN routine not tested");
//
//              Box coloredBox(loIV, hiIV);
//              FORT_GSRBLAZYNL(CHF_FRA1(a_phi[dit], 0),
//                            CHF_CONST_FRA1( lphi[dit], 0),
//                            CHF_CONST_FRA1(a_rhs[dit], 0),
//                            CHF_BOX(coloredBox),
//                            CHF_REAL(m_alpha),
//                            CHF_REAL(m_beta),
//                            CHF_CONST_REAL(m_gamma),
//                            CHF_REAL(m_dx));
//            }
//        }
//    } // end loop through red-black
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelJacobi(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::levelJacobi");
  MayDay::Error("AMRNonLinearPoissonOp::levelJacobi not implemented");

//  LevelData<FArrayBox> resid;
//  create(resid, a_rhs);
//
//  // Get the residual
//  residual(resid,a_phi,a_rhs,true);
//
//  // Create the weight
//  Real weight = m_alpha - 2.0*SpaceDim * m_beta / (m_dx*m_dx);
//
//  // Do the Jacobi relaxation
//  incr(a_phi, resid, -0.666/weight);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{
  CH_TIME("AMRNonLinearPoissonOp::CFInterp");

  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

//  bool FASmode = true;
  if (!m_homoCFinterp)
  {
    return;
  }

  DataIterator dit = a_phif.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterp(a_phif,datInd,idir,sit());
            }
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                                       const DataIndex&      a_datInd,
                                       int                   a_idir,
                                       Side::LoHiSide        a_hiorlo)
{
  // CH_TIME("AMRNonLinearPoissonOp::homogeneousCFInterp");

  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  //  CH_assert (m_ncomp == a_phif.nComp());

  const CFIVS* cfivs_ptr = NULL;
  if (a_hiorlo == Side::Lo)
    cfivs_ptr = &(m_cfregion.loCFIVS(a_datInd, a_idir));
  else
    cfivs_ptr = &(m_cfregion.hiCFIVS(a_datInd, a_idir));

  if (cfivs_ptr->isPacked())
    {
      int ihiorlo = sign(a_hiorlo);
      FArrayBox& phiFab = a_phif[a_datInd];
      if (phiFab.box().size(a_idir) == 3)
        {

          FORTNT_INTERPHOMOLINEAR(CHF_FRA(phiFab),
                                  CHF_BOX(cfivs_ptr->packedBox()),
                                  CHF_CONST_REAL(m_dx),
                                  CHF_CONST_REAL(m_dxCrse),
                                  CHF_CONST_INT(a_idir),
                                  CHF_CONST_INT(ihiorlo));
        }
      else
        {
          FORTNT_INTERPHOMO(CHF_FRA(phiFab),
                            CHF_BOX(cfivs_ptr->packedBox()),
                            CHF_CONST_REAL(m_dx),
                            CHF_CONST_REAL(m_dxCrse),
                            CHF_CONST_INT(a_idir),
                            CHF_CONST_INT(ihiorlo));
        }
    }
  else
    {
      const IntVectSet& interp_ivs = cfivs_ptr->getFineIVS();
      if (!interp_ivs.isEmpty())
        {
          // Assuming homogenous, interpolate on fine ivs
          interpOnIVSHomo(a_phif, a_datInd, a_idir,
                          a_hiorlo, interp_ivs);
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::singleBoxCFInterp(FArrayBox& a_phi)
{
  CH_TIME("AMRNonLinearPoissonOp::singleBoxCFInterp");

  Box region = a_phi.box();
  region.grow(-1);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
        {
          Box edge = adjCellBox(region, idir, sit(), 1);
          int ihiorlo = sign(sit());
          FORT_INTERPHOMO(CHF_FRA(a_phi),
                      CHF_BOX(edge),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_REAL(m_dxCrse),
                      CHF_CONST_INT(idir),
                      CHF_CONST_INT(ihiorlo));
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                                   const DataIndex&      a_datInd,
                                   const int             a_idir,
                                   const Side::LoHiSide  a_hiorlo,
                                   const IntVectSet&     a_interpIVS)
{
  // CH_TIME("AMRNonLinearPoissonOp::interpOnIVSHomo");

  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  IVSIterator fine_ivsit(a_interpIVS);
  FArrayBox& a_phi = a_phif[a_datInd];
  int ihilo = sign(a_hiorlo);

  if (a_phi.box().size(a_idir)==3) // we are in a 1-wide  box
    {
      IntVect iv;
      Real pa;
      Real factor = 1-2*m_dx/(m_dx+m_dxCrse);
      for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
        {
          iv = fine_ivsit();
          iv[a_idir]-=ihilo;
          // linear interpolation
          for (int ivar = 0; ivar < a_phif.nComp(); ivar++)
            {
              pa = a_phi(iv, ivar);
              a_phi(fine_ivsit(), ivar) = factor*pa;
            }
        }
    }
  else if (false) // the old expensive way of computing this  bvs
    {

      // much of these scalar values can be precomputed and stored if
      // we ever need to speed-up this function (ndk)
      Real x1 = m_dx;
      Real x2 = 0.5*(3. * x1 + m_dxCrse);
      Real denom = 1.0-((x1+x2)/x1);
      Real idenom = 1/(denom); // divide is more expensive usually
      Real x = 2.*x1;
      Real xsquared = x*x;

      Real m1 = 1/(x1*x1);
      Real m2 = 1/(x1*(x1-x2));

      Real q1 = 1/(x1-x2);
      Real q2 = x1+x2;

      Real pa,pb,a,b;
      IntVect ivf;
      for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
        {
          ivf = fine_ivsit();
          // quadratic interpolation
          for (int ivar = 0; ivar < a_phif.nComp(); ivar++)
            {
              ivf[a_idir]-=2*ihilo;
              pa = a_phi(ivf, ivar);
              ivf[a_idir]+=ihilo;
              pb = a_phi(ivf, ivar);

              a = ((pb-pa)*m1 - (pb)*m2)*idenom;
              b = (pb)*q1 - a*q2;

              ivf[a_idir]+=ihilo;
              a_phi(fine_ivsit(), ivar) = a*xsquared + b*x + pa;

            } //end loop over components
        } //end loop over fine intvects
    }
  else // symbolic reduced version of CF quadratic stencil
    {
      Real pa, pb;
      Real c1 = 2*(m_dxCrse-m_dx)/(m_dxCrse+  m_dx); //first inside point
      Real c2 =  -(m_dxCrse-m_dx)/(m_dxCrse+3*m_dx); // next point inward
      IntVect ivf;
      for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
        {
          ivf = fine_ivsit();
          // quadratic interpolation
          for (int ivar = 0; ivar < a_phif.nComp(); ivar++)
            {
              ivf[a_idir]-=2*ihilo;
              pa = a_phi(ivf, ivar);
              ivf[a_idir]+=ihilo;
              pb = a_phi(ivf, ivar);

              ivf[a_idir]+=ihilo;
              a_phi(fine_ivsit(), ivar) = c1*pb + c2*pa;

            } //end loop over components
        } //end loop over fine intvects
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::getFlux(FArrayBox&       a_flux,
                           const FArrayBox& a_data,
                           const Box&       a_edgebox,
                           int              a_dir,
                           int              a_ref) const
{
  // In this version of getFlux, the edgebox is passed in, and the flux array
  // is already defined.

  CH_TIME("AMRNonLinearPoissonOp::getFlux1");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  // if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!a_edgebox.isEmpty());

  Real scale = m_beta * a_ref / m_dx;

  FORT_NEWGETFLUXNL(CHF_FRA(a_flux),
                  CHF_CONST_FRA(a_data),
                  CHF_BOX(a_edgebox),
                  CHF_CONST_REAL(scale),
                  CHF_CONST_INT(a_dir));
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::getFlux(FArrayBox&       a_flux,
                           const FArrayBox& a_data,
                           int              a_dir,
                           int              a_ref) const
{
  CH_TIME("AMRNonLinearPoissonOp::getFlux2");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());

  Box edgebox = surroundingNodes(a_data.box(),a_dir);
  edgebox.grow(a_dir, -1);

  // if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!edgebox.isEmpty());

  a_flux.define(edgebox, a_data.nComp());

  //FArrayBox fflux(edgebox, a_data.nComp());
  Real scale = m_beta * a_ref / m_dx;

  FORT_NEWGETFLUXNL(CHF_FRA(a_flux),
                  CHF_CONST_FRA(a_data),
                  CHF_BOX(edgebox),
                  CHF_CONST_REAL(scale),
                  CHF_CONST_INT(a_dir));
}

// Factory

// ---------------------------------------------------------
//  AMR Factory define function
void AMRNonLinearPoissonOpFactory::define(const ProblemDomain&             a_coarseDomain,
                                 const Vector<DisjointBoxLayout>& a_grids,
                                 const Vector<int>&               a_refRatios,
                                 const Real&                      a_coarsedx,
                                 BCHolder                         a_bc,
                                 Real                             a_alpha,
                                 Real                             a_beta,
                                 Real                             a_gamma)
{
  CH_TIME("AMRNonLinearPoissonOpFactory::define");

  m_boxes = a_grids;

  m_refRatios = a_refRatios;

  m_bc = a_bc;

  m_dx.resize(a_grids.size());
  m_dx[0] = a_coarsedx;

  m_domains.resize(a_grids.size());
  m_domains[0] = a_coarseDomain;

  m_exchangeCopiers.resize(a_grids.size());
  m_exchangeCopiers[0].exchangeDefine(a_grids[0], IntVect::Unit);
  m_exchangeCopiers[0].trimEdges(a_grids[0], IntVect::Unit);

  m_cfregion.resize(a_grids.size());
  m_cfregion[0].define(a_grids[0], m_domains[0]);

  for (int i = 1; i < a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1] / m_refRatios[i-1];

      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);

      if (a_grids[i].isClosed())
        {
          m_exchangeCopiers[i].exchangeDefine(a_grids[i], IntVect::Unit);
          m_exchangeCopiers[i].trimEdges(a_grids[i], IntVect::Unit);
          m_cfregion[i].define(a_grids[i], m_domains[i]);
        }
    }

  m_alpha = a_alpha;
  m_beta = a_beta;
  m_gamma = a_gamma;
}

// ---------------------------------------------------------
// MultiGrid define function
void AMRNonLinearPoissonOpFactory::define(const ProblemDomain&     a_domain,
                                 const DisjointBoxLayout& a_grid,
                                 const Real&              a_dx,
                                 BCHolder                 a_bc,
                                 int                      a_maxDepth,
                                 Real                     a_alpha,
                                 Real                     a_beta,
                                 Real                     a_gamma)
{

  Vector<DisjointBoxLayout> grids(1);
  grids[0] = a_grid;
  Vector<int> refRatio(1, 2);
  define(a_domain, grids, refRatio, a_dx, a_bc, a_alpha, a_beta, a_gamma);
}

// ---------------------------------------------------------
MGLevelOp<LevelData<FArrayBox> >* AMRNonLinearPoissonOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                                                int                  a_depth,
                                                                bool                 a_homoOnly)
{
  CH_TIME("AMRNonLinearPoissonOpFactory::MGnewOp");

  Real dxCrse = -1.0;

  int ref;

  for (ref = 0; ref< m_domains.size(); ref++)
    {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox())
      {
        break;
      }
    }

  CH_assert(ref != m_domains.size()); // didn't find domain

  if (ref > 0)
    {
      dxCrse = m_dx[ref-1];
    }

  ProblemDomain domain(m_domains[ref]);
  Real dx = m_dx[ref];
  int coarsening = 1;

  for (int i = 0; i < a_depth; i++)
    {
      coarsening *= 2;
      domain.coarsen(2);
    }

  if (coarsening > 1 && !m_boxes[ref].coarsenable(coarsening*AMRNonLinearPoissonOp::s_maxCoarse))
  {
    return NULL;
  }

  dx *= coarsening;

  DisjointBoxLayout layout;
  coarsen_dbl(layout, m_boxes[ref], coarsening);

  Copier ex = m_exchangeCopiers[ref];
  CFRegion cfregion = m_cfregion[ref];

  if (coarsening > 1)
    {
      ex.coarsen(coarsening);
      cfregion.coarsen(coarsening);
    }

  AMRNonLinearPoissonOp* newOp = new AMRNonLinearPoissonOp;
  newOp->define(layout, dx, domain, m_bc, ex, cfregion);

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;
  newOp->m_gamma = m_gamma;

  newOp->m_aCoef = m_alpha;
  newOp->m_bCoef = m_beta;

  newOp->m_dxCrse = dxCrse;

  return (MGLevelOp<LevelData<FArrayBox> >*)newOp;
}

// ---------------------------------------------------------
AMRLevelOp<LevelData<FArrayBox> >* AMRNonLinearPoissonOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  CH_TIME("AMRNonLinearPoissonOpFactory::AMRnewOp");

  AMRNonLinearPoissonOp* newOp = new AMRNonLinearPoissonOp;
  Real dxCrse = -1.0;

  int ref;

  for (ref = 0; ref< m_domains.size(); ref++)
    {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox())
      {
        break;
      }
    }

  if (ref == 0)
    {
      // coarsest AMR level
      if ((m_domains.size() == 1) || (!m_boxes[1].isClosed()) )
        {
          // no finer level
          newOp->define(m_boxes[0], m_dx[0],
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0]);
        }
      else
        {
          // finer level exists but no coarser
          int dummyRat = 1;  // argument so compiler can find right function
          int refToFiner = m_refRatios[0]; // actual refinement ratio
          newOp->define(m_boxes[0], m_boxes[1], m_dx[0],
                        dummyRat, refToFiner,
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0]);
        }
    }
  else if ((ref ==  m_domains.size()-1) || (!m_boxes[ref+1].isClosed()))
    {
      dxCrse = m_dx[ref-1];

      // finest AMR level
      newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1],
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[ref], m_cfregion[ref]);
    }
  else if ( ref == m_domains.size())
    {
      MayDay::Abort("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");

    }
  else
    {
      dxCrse = m_dx[ref-1];

      // intermediate AMR level, full define
      newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1], m_refRatios[ref],
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[ref], m_cfregion[ref]);
    }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;
  newOp->m_gamma = m_gamma;

  newOp->m_aCoef = m_alpha;
  newOp->m_bCoef = m_beta;

  newOp->m_dxCrse = dxCrse;

  return (AMRLevelOp<LevelData<FArrayBox> >*)newOp;
}

// ---------------------------------------------------------
int AMRNonLinearPoissonOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if (m_domains[ilev].domainBox() == a_domain.domainBox())
        {
          retval = m_refRatios[ilev];
          found = true;
        }
    }

  if (!found)
    {
      MayDay::Abort("Domain not found in AMR hierarchy");
    }

  return retval;
}

#include "NamespaceFooter.H"
