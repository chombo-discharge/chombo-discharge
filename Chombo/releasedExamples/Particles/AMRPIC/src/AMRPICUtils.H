#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRPICUTILS_H_
#define _AMRPICUTILS_H_

#include <cmath>

#include "AMR.H"
#include "AMRLevelPICFactory.H"
#include "memtrack.H"
#include "parstream.H"
#include "ZeldovichIBC.H"

#include "NamespaceHeader.H"

///
/**
 */
extern void
getProblemDomain(ProblemDomain& a_domain);

///
/**
 */
extern ZeldovichIBC*
getParticleIBC();

///
/**
 */
extern void
getAMRLevelPICFactory(RefCountedPtr<AMRLevelPICFactory>&  a_fact,
		      ZeldovichIBC* a_ibc);

///
/**
 */
extern void
setupAMRHierarchy(Vector<Vector<Box> >& a_hierarchy);

///
/**
 */
extern void
defineAMR(AMR&                                 a_amr,
          RefCountedPtr<AMRLevelPICFactory>&   a_fact,
          const ProblemDomain&                 a_prob_domain,
          const Vector<int>&                   a_refRat);

///
/**
 */
extern void
setupAMRForFixedHierarchyRun(AMR& a_amr,
			     const Vector<Vector<Box> >& a_hierarchy);

/****/

#include "NamespaceFooter.H"

#endif
