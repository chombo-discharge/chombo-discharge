#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELSOURCETERM_H_
#define _LEVELSOURCETERM_H_

#include "MultiBlockCoordSys.H"
#include "MOLPhysics.H"
#include "LevelData.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

class LevelSourceTerm
{
public:

  /// Default constructor
  LevelSourceTerm();

  /// Destructor
  virtual ~LevelSourceTerm();

  /// Define the object so that time stepping can begin (actual constructor)
  virtual void define(MultiBlockCoordSys* const a_coordSysPtr,
                      const MOLPhysics* const a_molPhysics,
                      const DisjointBoxLayout& a_grids,
                      const RealVect         & a_dx)
  {
    define(a_coordSysPtr, a_molPhysics, a_grids);
  }

  /// Define the object so that time stepping can begin (actual constructor)
  virtual void define(MultiBlockCoordSys* const a_coordSysPtr,
                      const MOLPhysics* const a_molPhysics,
                      const DisjointBoxLayout& a_grids);

  /// Object factory for this class.  Pure virtual function.
  /**
   */
  virtual LevelSourceTerm* new_sourceTerm() const = 0;

  /// Increment a_rhs by source term, using a_U.  Pure = 0; virtual function.
  virtual void addSourceTerm(LevelData<FArrayBox>&   a_rhs,
                             LevelData<FArrayBox>&   a_U)
  {
    MayDay::Error("not implemented");
  }

  virtual void addSourceTerm(LevelData<FArrayBox>&   a_rhs,
                             LevelData<FArrayBox>&   a_U,
                             const RealVect      &   a_dx)
  {
    addSourceTerm(a_rhs, a_U);
  }

protected:

  bool m_isDefined;

  MultiBlockCoordSys* m_coordSysPtr;

  const MOLPhysics* m_molPhysics;

  DisjointBoxLayout m_grids;

private:

  // Disallowed for all the usual reasons
  void operator=(const LevelSourceTerm&);
  LevelSourceTerm(const LevelSourceTerm&);
};

#include "NamespaceFooter.H"

#endif
