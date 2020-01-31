#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Gradient class: utility for computing grad(phi) for difference stencils
// written by F. Miniati

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "ProblemDomain.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "Box.H"

enum StencilType
{
  TwoPts, FourPts, TenPts
};

class Gradient
{
public:
  /// Constructor
  Gradient()
  {
  }

  /// Destructor
  virtual ~Gradient()
  {
  }

  // computed cell centered gradient
  virtual void gradient(FArrayBox&           a_grad,
                        const FArrayBox&     a_phi,
                        const ProblemDomain& a_domain,
                        const Real&          a_dx,
                        const Box&           a_box) = 0;
};

class TwoPtsGradient : public Gradient
{
public:
  /// Constructor
  TwoPtsGradient()
    :
    Gradient()
  {
  }

  /// Destructor
  virtual ~TwoPtsGradient()
  {
  }

  // computed cell centered gradient
  virtual void gradient(FArrayBox&           a_grad,
                        const FArrayBox&     a_phi,
                        const ProblemDomain& a_domain,
                        const Real&          a_dx,
                        const Box&           a_box);
};

class FourPtsGradient : public Gradient
{
public:
  /// Constructor
  FourPtsGradient()
    :
    Gradient()
  {
  }

  /// Destructor
  virtual ~FourPtsGradient()
  {
  }

  // computed cell centered gradient
  virtual void gradient(FArrayBox&           a_grad,
                        const FArrayBox&     a_phi,
                        const ProblemDomain& a_domain,
                        const Real&          a_dx,
                        const Box&           a_box);
};

class TenPtsGradient : public Gradient
{
public:
  /// Constructor
  TenPtsGradient()
    :
    Gradient()
  {
  }

  /// Destructor
  virtual ~TenPtsGradient()
  {
  }

  // computed cell centered gradient
  virtual void gradient(FArrayBox&           a_grad,
                        const FArrayBox&     a_phi,
                        const ProblemDomain& a_domain,
                        const Real&          a_dx,
                        const Box&           a_box);
};

#endif
