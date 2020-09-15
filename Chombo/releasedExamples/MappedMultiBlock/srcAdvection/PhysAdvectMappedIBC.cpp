#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PhysAdvectMappedIBC.H"
#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
PhysAdvectMappedIBC::PhysAdvectMappedIBC() : PhysMappedIBC()
{
}

PhysAdvectMappedIBC::~PhysAdvectMappedIBC()
{
}

void
PhysAdvectMappedIBC::setBdrySlopes(FArrayBox&       a_dW,
                                   const FArrayBox& a_W,
                                   const int&       a_dir,
                                   const Real&      a_time)
{
  MayDay::Error("PhysAdvectMappedIBC::setBdrySlopes not defined");
}

#include "NamespaceFooter.H"
