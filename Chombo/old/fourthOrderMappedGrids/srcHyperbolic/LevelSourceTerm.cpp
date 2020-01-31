#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelSourceTerm.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////

LevelSourceTerm::LevelSourceTerm()
{
  m_isDefined = false;
  m_coordSysPtr = NULL;
  m_molPhysics = NULL;
}

//////////////////////////////////////////////////////////////////////////////

LevelSourceTerm::~LevelSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void LevelSourceTerm::define(MultiBlockCoordSys* const a_coordSysPtr,
                             const MOLPhysics* const a_molPhysics,
                             const DisjointBoxLayout& a_grids)
{
  m_coordSysPtr = a_coordSysPtr;
  m_molPhysics = a_molPhysics;
  m_grids = a_grids;
  m_isDefined = true;
}

//////////////////////////////////////////////////////////////////////////////

#include "NamespaceFooter.H"
