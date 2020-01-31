#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NWOEBVTOLoadBalance_H_
#define _NWOEBVTOLoadBalance_H_

#include "Box.H"
#include "ProblemDomain.H"
#include "Vector.H"
#include "REAL.H"
#include "EBIndexSpace.H"
#include "NWOEBViscousTensorOpFactory.H"
#include "DisjointBoxLayout.H"
#include "NamespaceHeader.H"

///
/**
    Load balance a Vector of Boxes.
    The order of boxes can change so that is why it is output as well as input.
    This is introspective load balancing using NWOEBViscousTensorOp::relax times
    for loads.
 */
extern int
NWOEBVTORelaxLoadBalance(Vector<int>                              &  a_procs,
                         Vector<Box>                              &  a_boxes,
                         const RefCountedPtr<BaseDomainBCFactory> &  a_domBC,
                         const RefCountedPtr<BaseEBBCFactory>     &  a_ebBC,
                         const ProblemDomain                      &  a_domain,
                         const EBIndexSpace *a_ebisPtr = Chombo_EBIS::instance() );
///
/**
 */
extern void
getNWOEBVTORelaxLoadsAndBoxes(Vector<unsigned long long>               &    a_loads,
                              Vector<Box>                              &    a_boxes,
                              const RefCountedPtr<BaseDomainBCFactory> &    a_domBC,
                              const RefCountedPtr<BaseEBBCFactory>     &    a_ebBC,
                              const DisjointBoxLayout                  &    a_dblOrig,
                              const ProblemDomain                      &    a_domain,
                              const EBIndexSpace *a_ebisPtr = Chombo_EBIS::instance() );


#include "NamespaceFooter.H"
#endif
