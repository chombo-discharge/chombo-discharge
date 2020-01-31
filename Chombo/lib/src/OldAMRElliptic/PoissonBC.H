#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Mon, July 19, 1999

#ifndef _POISSONBC_H_
#define _POISSONBC_H_

#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "GhostBC.H"
#include "NamespaceHeader.H"

/// Class to enforce homogeneous Neumann boundary conditions
/**
    BoxGhostBC class implentation to enforce homogeneous neumann
    boundary conditions
    Recall that  BoxGhostBC is a class to encapsulate the
    operations of ghost-cell boundary conditions at a face.
    If the solution is phi and the face normal direction is x,
    the boundary conditions usually used can be expressed as
    A*phi + B*dphi/dx = C.
    For homogeneous Neumann BCs, A=C=0, B= 1
    These functions are only meant to be called by DomainGhostBC
 */
class NeumannBC: public BoxGhostBC
{
public:
  friend class DomainGhostBC;

  /// Null constructor
  NeumannBC();

  /// Basic constructor -- sets Interval to (0,0)
  NeumannBC(int dir,
            Side::LoHiSide sd);

  /// Full constructor
  NeumannBC(int             dir,
            Side::LoHiSide  sd,
            const Interval& a_comps);

  /// Destructor
  ~NeumannBC();

  void setInhomVal(Real a_inhomVal);
protected:
  ///
  /**
     Fills coefficient arrays --
     fills neumann factor with 1.0
     and other two fields with 0.0.
   */
  virtual void fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real       a_dx,
                            const Box& a_domain) const;

  ///
  /**
     Fills coefficient arrays --
     fills neumann factor with 1.0
     and other two fields with 0.0.
   */
  virtual void fillBCValues(FArrayBox&           a_neumfac,
                            FArrayBox&           a_dircfac,
                            FArrayBox&           a_inhmval,
                            Real                 a_dx,
                            const ProblemDomain& a_domain) const;

  ///
  /**
     Virtual constructor workaround.
  */
  BoxGhostBC* new_boxghostbc() const;

private:
  Real m_inhomVal;
  NeumannBC(const NeumannBC& a_neuin)
  {
  }

  virtual void operator=(const NeumannBC& a_neuin)
  {
  }
};

/// Class to enforce homogeneous Dirichlet boundary conditions
/**
    BoxGhostBC class implentation to enforce homogeneous dirchelet
    boundary conditions
    Recall that
    BoxGhostBC is a class to encapsulate the
    operations  of ghost-cell boundary conditions at a face.
    If the solution is phi and the face normal direction is x,
    the boundary conditions usually used can be expressed as
    A*phi + B*dphi/dx = C.
    For homogeneous Dirichlet BCs, B=C=0, A= 1
    These functions are only meant to be called by DomainGhostBC
 */
class DirichletBC: public BoxGhostBC
{
public:
  friend class DomainGhostBC;

  /// Null constructor
  DirichletBC();

  /// Simple constructor (sets Interval to (0,0)
  DirichletBC(int dir, Side::LoHiSide sd);

  /// Full constructor
  DirichletBC(int dir, Side::LoHiSide sd,
              const Interval& a_comps);

  /// Destructor
  ~DirichletBC();

  void setInhomVal(Real a_inhomVal);
private:
  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with 0.0
   */
  virtual void fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real       a_dx,
                            const Box& a_domain) const;

  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with 0.0
   */
  virtual void fillBCValues(FArrayBox&           a_neumfac,
                            FArrayBox&           a_dircfac,
                            FArrayBox&           a_inhmval,
                            Real                 a_dx,
                            const ProblemDomain& a_domain) const;

  ///
  BoxGhostBC* new_boxghostbc() const;

private:
  Real m_inhomVal;

  DirichletBC(const DirichletBC& a_dircin)
  {
  }

  virtual void operator=(const DirichletBC& a_dircin)
  {
  }
};

#include "NamespaceFooter.H"
#endif
