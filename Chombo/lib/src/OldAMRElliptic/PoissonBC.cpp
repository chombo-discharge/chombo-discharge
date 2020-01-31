#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Box.H"
#include "FArrayBox.H"
#include "REAL.H"
#include "SPACE.H"
#include "Tuple.H"
#include "Vector.H"
#include "LoHiSide.H"

#include "PoissonBC.H"
#include "MayDay.H"
#include "NamespaceHeader.H"

NeumannBC::NeumannBC() : BoxGhostBC(), m_inhomVal(0.0)
{
}

NeumannBC::~NeumannBC()
{
}
void NeumannBC::setInhomVal(Real a_inhomVal)
{
  m_inhomVal = a_inhomVal;
}
void DirichletBC::setInhomVal(Real a_inhomVal)
{
  m_inhomVal = a_inhomVal;
}
void NeumannBC::fillBCValues(FArrayBox& a_neumfac,
                             FArrayBox& a_dircfac,
                             FArrayBox& a_inhmval,
                             Real       a_dx,
                             const Box& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());

  a_neumfac.setVal(1.0);
  a_dircfac.setVal(0.0);
  a_inhmval.setVal(m_inhomVal);
}

void NeumannBC::fillBCValues(FArrayBox&           a_neumfac,
                             FArrayBox&           a_dircfac,
                             FArrayBox&           a_inhmval,
                             Real                 a_dx,
                             const ProblemDomain& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());

  a_neumfac.setVal(1.0);
  a_dircfac.setVal(0.0);
  a_inhmval.setVal(m_inhomVal);
}

void DirichletBC::fillBCValues(FArrayBox& a_neumfac,
                               FArrayBox& a_dircfac,
                               FArrayBox& a_inhmval,
                               Real       a_dx,
                               const Box& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());

  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(m_inhomVal);
}

void DirichletBC::fillBCValues(FArrayBox&           a_neumfac,
                               FArrayBox&           a_dircfac,
                               FArrayBox&           a_inhmval,
                               Real                 a_dx,
                               const ProblemDomain& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());

  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(m_inhomVal);
}

BoxGhostBC* NeumannBC::new_boxghostbc() const
{
  NeumannBC* newop = new NeumannBC();

  if (newop == NULL)
    {
      MayDay::Error("Out of Memory in NeumannBC::new_boxghostbc");
    }
  newop->setInhomVal(m_inhomVal);
  return static_cast<BoxGhostBC*>(newop);
}

BoxGhostBC* DirichletBC::new_boxghostbc() const
{
  DirichletBC* newop = new DirichletBC();

  if (newop == NULL)
    {
      MayDay::Error("Out of Memory in DirichletBC::new_boxghostbc");
    }

  newop->setInhomVal(m_inhomVal);
  return static_cast<BoxGhostBC*>(newop);
}

NeumannBC::NeumannBC(int             a_dir,
                     Side::LoHiSide  a_side,
                     const Interval& a_comps)
  :BoxGhostBC(a_dir,a_side,a_comps),m_inhomVal(0.0)
{
}

NeumannBC::NeumannBC(int a_dir,
                     Side::LoHiSide a_side)
  :BoxGhostBC(a_dir,a_side), m_inhomVal(0.0)
{
}

DirichletBC::DirichletBC()
  :BoxGhostBC(), m_inhomVal(0.0)
{
}

DirichletBC::~DirichletBC()
{
}

DirichletBC::DirichletBC(int dir, Side::LoHiSide side)
  :BoxGhostBC(dir,side), m_inhomVal(0.0)
{
}

DirichletBC::DirichletBC(int             a_dir,
                         Side::LoHiSide  a_side,
                         const Interval& a_comps)
  :BoxGhostBC(a_dir,a_side,a_comps), m_inhomVal(0.0)
{
}
#include "NamespaceFooter.H"
