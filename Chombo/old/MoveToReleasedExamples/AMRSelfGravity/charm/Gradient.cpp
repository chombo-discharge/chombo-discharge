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

#include "parstream.H"
#include "Gradient.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"
#include "GradientF_F.H"

// compute cell centered gradient
void TwoPtsGradient::gradient(FArrayBox&           a_grad,
                              const FArrayBox&     a_phi,
                              const ProblemDomain& a_domain,
                              const Real&          a_dx,
                              const Box&           a_box)
{
#ifndef NDEBUG
  Box box = grow(a_box,1);
  CH_assert(a_phi.box().contains(box));
#endif

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  a_grad.setVal(0.0);
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    Box gbox = grow(a_box, BASISV(dir));
    //     Box gbox = a_box;
    //     gbox.grow(dir,1);
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
               gbox,a_domain,dir);

    FORT_TWOPTSGRAD(CHF_FRA1(a_grad,dir),
                     CHF_CONST_FRA1(a_phi,0),
                     CHF_CONST_INT(dir),
                     CHF_CONST_REAL(a_dx),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));
  }
}

// compute cell centered gradient
void FourPtsGradient::gradient(FArrayBox&           a_grad,
                               const FArrayBox&     a_phi,
                               const ProblemDomain& a_domain,
                               const Real&          a_dx,
                               const Box&           a_box)
{
#ifndef NDEBUG
  Box box = grow(a_box,2);
  CH_assert(a_phi.box().contains(box));
#endif

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  a_grad.setVal(0.0);
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    // loHiCenter assumes three point stencil; fix this "by hand"
    //     Box gbox = a_box;
    //     gbox.grow(dir,1);
    Box gbox = grow(a_box, BASISV(dir));
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
               gbox,a_domain,dir);


     FORT_FOURPTSGRAD(CHF_FRA1(a_grad,dir),
                     CHF_CONST_FRA1(a_phi,0),
                     CHF_CONST_INT(dir),
                     CHF_CONST_REAL(a_dx),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));
  }
}

// compute cell centered gradient
void TenPtsGradient::gradient(FArrayBox&           a_grad,
                              const FArrayBox&     a_phi,
                              const ProblemDomain& a_domain,
                              const Real&          a_dx,
                              const Box&           a_box)
{
#ifndef NDEBUG
  Box box = grow(a_box,1);
  CH_assert(a_phi.box().contains(box));
#endif

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  a_grad.setVal(0.0);
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    Box gbox = grow(a_box, BASISV(dir));
    //     Box gbox = a_box;
    //     gbox.grow(dir,1);
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
               gbox,a_domain,dir);

    FORT_TENPTSGRAD(CHF_FRA1(a_grad,dir),
                     CHF_CONST_FRA1(a_phi,0),
                     CHF_CONST_INT(dir),
                     CHF_CONST_REAL(a_dx),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));
  }
}


