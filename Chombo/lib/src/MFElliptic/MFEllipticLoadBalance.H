#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFELLIPTICLOADBALANCE_H_
#define _MFELLIPTICLOADBALANCE_H_

#include "REAL.H"
#include "Vector.H"
#include "Box.H"
#include "ProblemDomain.H"
#include "RefCountedPtr.H"
#include "MFIndexSpace.H"
#include "NamespaceHeader.H"

class DisjointBoxLayout;
///
/**
    Load balance a Vector of Boxes.
    We  make the loads bases on a few ad hoc ideas.
    we shall make fully covered boxes = constant load = covered load
    we shall say that irregular points get irregular factor more load than
    regular points.
 */
extern int
MFEllipticLoadBalance(Vector<int>&                 a_proc,
                      const Vector<Box>&           a_boxes,
                      RefCountedPtr<MFIndexSpace>& a_mfis,
                      const ProblemDomain&         a_domain,
                      bool                         a_verbose = false);

///
/**
 */
extern void
getPoissonLoadsAndBoxes(Vector<unsigned long long>&   a_loads,
                        Vector<Box> &                 a_boxes,
                        RefCountedPtr<MFIndexSpace>&  a_mfis,
                        const DisjointBoxLayout&      a_dblOrig,
                        const ProblemDomain&          a_domain);

#include "NamespaceFooter.H"
#endif
