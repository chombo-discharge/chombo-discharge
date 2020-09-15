#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRPoissonOpF_F.H"
#include "PoissonOp4.H"
#include "NewPoissonOp4F_F.H"
#include "NamespaceHeader.H"

void PoissonOp4::define(const RealVect& a_dx,
                        const ProblemDomain& a_domain,
                        BCFunc a_bc)
{
  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
}

void PoissonOp4::residual(LevelData<FArrayBox>& a_lhs, 
                          const LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>& a_rhs, 
                          bool a_homogeneous)
{
  applyOp(a_lhs, a_phi, a_homogeneous);
  m_levelOps.scale(a_lhs, -1.0);
  incr(a_lhs, a_rhs, 1);
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void PoissonOp4::preCond(LevelData<FArrayBox>& a_phi, 
                         const LevelData<FArrayBox>& a_rhs)
{
  
// diagonal term of this operator is 4/h/h in 2D, 6/h/h in 3D,
  // so inverse of this is our initial multiplier
  Real mult = -m_dx[0]*m_dx[0]/(2.0*SpaceDim);

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

void PoissonOp4::applyOp(LevelData<FArrayBox>& a_lhs, 
                         const LevelData<FArrayBox>& a_phi,
                         bool a_homogeneous )
{

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  phi.exchange();
  Real dx = m_dx[0];

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  int nbox=dit.size();
#pragma omp parallel   default (shared)
  {
    CH_TIME("AMRPoissonOp::applyOpIBC");
#pragma omp for
    for (int ibox=0;ibox<nbox; ibox++)
      {
        m_bc(phi[dit[ibox]], dbl[dit[ibox]], m_domain, m_dx[0], a_homogeneous);
      }
  }
// end pragma


  Real alpha = 0.0;
  Real beta = 1.0;
#pragma omp parallel
  {
#pragma omp for
    for (int ibox=0;ibox<nbox; ibox++)
      {
        const Box& region = dbl[dit[ibox]];

        FORT_OPERATORLAP4(CHF_FRA(a_lhs[dit[ibox]]),
                          CHF_CONST_FRA(phi[dit[ibox]]),
                          CHF_BOX(region),
                          CHF_CONST_REAL(dx),
                          CHF_CONST_REAL(alpha),
                          CHF_CONST_REAL(beta));
      }
  }//end pragma
}

void PoissonOp4::createCoarsened(LevelData<FArrayBox>& a_lhs,
                                 const LevelData<FArrayBox>& a_rhs,
                                 const int& a_refRat)
{
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

  a_lhs.define(dblCoarsenedFine, ncomp, ghostVect);
}

void PoissonOp4::create(LevelData<FArrayBox>& a_lhs, 
                        const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.create(a_lhs, a_rhs);
}

void PoissonOp4::assign(LevelData<FArrayBox>& a_lhs, 
                        const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.assign(a_lhs, a_rhs);
}

Real PoissonOp4::dotProduct(const LevelData<FArrayBox>& a_1, 
                            const LevelData<FArrayBox>& a_2)
{
  return m_levelOps.dotProduct(a_1, a_2);
}

void PoissonOp4::incr(LevelData<FArrayBox>& a_lhs, 
                      const LevelData<FArrayBox>& a_x, 
                      Real a_scale)
{
  m_levelOps.incr(a_lhs, a_x, a_scale);
}

void PoissonOp4::axby(LevelData<FArrayBox>& a_lhs, 
                      const LevelData<FArrayBox>& a_x,
                      const LevelData<FArrayBox>& a_y, Real a_a, Real a_b)
{
  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}

void PoissonOp4::scale(LevelData<FArrayBox>& a_lhs, 
                       const Real& a_scale)
{
  m_levelOps.scale(a_lhs, a_scale);
}

void PoissonOp4::setToZero(LevelData<FArrayBox>& a_lhs)
{
  m_levelOps.setToZero(a_lhs);
}

Real PoissonOp4::norm(const LevelData<FArrayBox>& a_x, int a_ord)
{
  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}

void PoissonOp4::relax(LevelData<FArrayBox>& a_e,
                       const LevelData<FArrayBox>& a_residual,
                       int a_iterations)
{
  for (int i=0; i<a_iterations; i++)
  {
    levelGSRB(a_e, a_residual);
  }
}

void PoissonOp4::createCoarser(LevelData<FArrayBox>& a_coarse,
                               const LevelData<FArrayBox>& a_fine,
                               bool ghosted)
{
  IntVect ghost = a_fine.ghostVect();
  CH_assert(a_fine.disjointBoxLayout().coarsenable(2));
  if (m_coarsenedMGrids.size() == 0)
    coarsen(m_coarsenedMGrids, a_fine.disjointBoxLayout(), 2); // multigrid, so coarsen by 2
  a_coarse.define(m_coarsenedMGrids, a_fine.nComp(), ghost);
}

void PoissonOp4::restrictResidual(LevelData<FArrayBox>& a_resCoarse,
                                  LevelData<FArrayBox>& a_phiFine,
                                  const LevelData<FArrayBox>& a_rhsFine)
{

  a_phiFine.exchange();

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  DataIterator dit = a_phiFine.dataIterator();
  int nbox=dit.size();
#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi = a_phiFine[dit[ibox]];
        m_bc(phi, dblFine[dit[ibox]], m_domain, m_dx[0], true);
      }
  }//end pragma

  Real alpha = 0.0, beta = 1.0;
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

        FORT_RESTRICTRES4(CHF_FRA_SHIFT(res, civ),
                          CHF_CONST_FRA_SHIFT(phi, iv),
                          CHF_CONST_FRA_SHIFT(rhs, iv),
                          CHF_CONST_REAL(alpha),
                          CHF_CONST_REAL(beta),
                          CHF_BOX_SHIFT(region, iv),
                          CHF_CONST_REAL(m_dx[0]));
      }
  }//end pragma
}

void PoissonOp4::prolongIncrement(LevelData<FArrayBox>& a_phiThisLevel,
                                  const LevelData<FArrayBox>& a_correctCoarse)
{

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
        
        FORT_PROLONG(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(mgref));
      }
  }//end pragma
}

/***/
void PoissonOp4::levelGSRB(LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  const DisjointBoxLayout& dbl = a_rhs.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();
  int nbox=dit.size();

  LevelData<FArrayBox> tmp;
  create(tmp, a_phi);
  assign(tmp, a_phi);

  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {
      a_phi.exchange();
      tmp.exchange();
#pragma omp parallel
      {
#pragma omp for 
        for (int ibox=0; ibox < nbox; ibox++)
          {
            const Box& region = dbl[dit[ibox]];
            FArrayBox& phiFab = a_phi[dit[ibox]];
            
            if (whichPass==0)
              m_bc( phiFab, region, m_domain, m_dx[0], true );
            else
              m_bc(tmp[dit[ibox]],  region, m_domain,  m_dx[0], true);
            FORT_GSRBLAPLACIAN4(CHF_FRA(phiFab),
                                CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                CHF_BOX(region),
                                CHF_CONST_REAL(m_dx[0]),
                                CHF_FRA(tmp[dit[ibox]]),
                                CHF_CONST_INT(whichPass));
          } // end loop through grids
      }//end pragma
    } // end loop through red-black
}

// Factory

//----------------------------------------------------------
PoissonOp4Factory::PoissonOp4Factory()
{
  m_dx  = RealVect::Unit;
  m_bc  = NULL;
}

PoissonOp4Factory::PoissonOp4Factory(RealVect& a_dx, BCFunc a_bc)
  :m_dx(a_dx), m_bc(a_bc)
{
}

void PoissonOp4Factory::define(const RealVect& a_dx,  BCFunc a_bc)
{
  m_dx=a_dx;
  m_bc=a_bc;
}

PoissonOp4* PoissonOp4Factory::MGnewOp(const ProblemDomain& a_FineindexSpace,
                                       int   a_depth,
                                       bool  a_homoOnly)
{
  CH_assert(a_depth >= 0 );
  CH_assert(m_bc != NULL);
  PoissonOp4* newOp = new PoissonOp4();
  RealVect dx = m_dx;
  ProblemDomain domain = a_FineindexSpace;
  for (int i=0; i<a_depth; i++)
    {
      Box d = domain.domainBox();
      d.coarsen(8);
      d.refine(8);
      if (domain.domainBox() == d)
        {
          dx*=2;
          domain.coarsen(2);
        }
      else
        {
          return NULL;
        }
    }
  newOp->define(dx, domain, m_bc);
  return newOp;
}

void PoissonOp4Factory::MGreclaim(PoissonOp4* a_reclaim)
{
  delete a_reclaim;
}

#include "NamespaceFooter.H"
