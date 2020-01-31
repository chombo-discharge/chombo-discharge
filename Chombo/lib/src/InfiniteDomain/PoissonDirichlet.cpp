#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// PoissonDirichlet.cpp
// petermc, 19 Aug 2004

#include "PoissonDirichlet.H"
#include "InfiniteDomainF.H"
#include "MayDay.H"
#include "EvalOperatorF_F.H"

using std::cout;
using std::endl;
using std::cerr;

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
PoissonDirichlet::PoissonDirichlet()
{
  m_isDefined = false;
}


// ---------------------------------------------------------
PoissonDirichlet::PoissonDirichlet(const Box&                         a_bx,
                                   const RealVect&                    a_dx,
                                   OperatorType                       a_op)
{
  define(a_bx, a_dx, a_op);
}


// ---------------------------------------------------------
PoissonDirichlet::~PoissonDirichlet()
{
}


// ---------------------------------------------------------
void
PoissonDirichlet::define(const Box&                         a_bx,
                         const RealVect&                    a_dx,
                         OperatorType                       a_op)
{
  m_bx = a_bx;
  m_dx = a_dx;
  m_op = a_op;
  m_verbose = 0;
  m_isDefined = true;
}

// ---------------------------------------------------------
bool PoissonDirichlet::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void PoissonDirichlet::setVerbose(int a_verbose)
{
  m_verbose = a_verbose;
}


// ---------------------------------------------------------
void PoissonDirichlet::solveHomogeneous(FArrayBox&         a_phi,
                                        const FArrayBox&   a_rhs)
{
  CH_assert(isDefined());
  CH_assert(a_phi.box() == m_bx);
  CH_assert(a_rhs.box() == m_bx);

  int n0 = m_bx.size(0) - 1;
  int n1 = m_bx.size(1) - 1;
  int n2 = m_bx.size(2) - 1;

  int intOp = (int) m_op;
  a_phi.copy(a_rhs);
  FORT_POISSON3D(a_phi.dataPtr(), &n0, &n1, &n2, m_dx.dataPtr(), &intOp);
}


// ---------------------------------------------------------
void PoissonDirichlet::solveHomogeneousInPlace(FArrayBox&   a_phirhs)
{
  CH_assert(isDefined());
  CH_assert(a_phirhs.box() == m_bx);

  int n0 = m_bx.size(0) - 1;
  int n1 = m_bx.size(1) - 1;
  int n2 = m_bx.size(2) - 1;

  int intOp = (int) m_op;
  FORT_POISSON3D(a_phirhs.dataPtr(), &n0, &n1, &n2, m_dx.dataPtr(), &intOp);
}


// ---------------------------------------------------------
void PoissonDirichlet::solveInhomogeneous(FArrayBox&         a_phi,
                                          const FArrayBox&   a_rhs)
{
  CH_assert(isDefined());
  CH_assert(a_phi.box() == m_bx);
  CH_assert(a_rhs.box() == m_bx);

  // Copy a_rhs to a_phi in interior.
  Box bxInner = grow(m_bx, -1);
  a_phi.copy(a_rhs, bxInner);

  int intOp = (int) m_op;
  int errcode = 0;
  const int *bxLo = m_bx.loVect();
  const int *bxHi = m_bx.hiVect();
  FORT_INHOMODIRICHLET3D(a_phi.dataPtr(),
                         &bxLo[0], &bxLo[1], &bxLo[2],
                         &bxHi[0], &bxHi[1], &bxHi[2],
                         &m_dx[0], &m_dx[1], &m_dx[2],
                         &intOp, &errcode);
  if (errcode != 0)
    {
      cerr << "inhomodirichlet3d returned error code " << errcode << endl;
      MayDay::Error("returning");
      return;
    }
}


// ---------------------------------------------------------
void PoissonDirichlet::solveInhomogeneousInPlace(FArrayBox&   a_phirhs)
{
  CH_assert(isDefined());
  CH_assert(a_phirhs.box() == m_bx);

  int intOp = (int) m_op;
  int errcode = 0;
  const int *bxLo = m_bx.loVect();
  const int *bxHi = m_bx.hiVect();
  FORT_INHOMODIRICHLET3D(a_phirhs.dataPtr(),
                         &bxLo[0], &bxLo[1], &bxLo[2],
                         &bxHi[0], &bxHi[1], &bxHi[2],
                         &m_dx[0], &m_dx[1], &m_dx[2],
                         &intOp, &errcode);
  if (errcode != 0)
    {
      cerr << "inhomodirichlet3d returned error code " << errcode << endl;
      MayDay::Error("returning");
      return;
    }
}


// ---------------------------------------------------------
void PoissonDirichlet::evalOperator(FArrayBox&         a_eval,
                                    const FArrayBox&   a_phi)
{
  CH_assert(isDefined());
  CH_assert(a_eval.box() == m_bx);
  CH_assert(a_phi.box() == m_bx);

  int intOp = (int) m_op;

  a_eval.setVal(0.);
  Box bxInner = grow(m_bx, -1);
  FORT_EVALOP(CHF_FRA1(a_eval, 0),
              CHF_CONST_INT(intOp),
              CHF_CONST_REALVECT(m_dx),
              CHF_BOX(bxInner),
              CHF_CONST_FRA1(a_phi, 0));
}

#include "NamespaceFooter.H"
