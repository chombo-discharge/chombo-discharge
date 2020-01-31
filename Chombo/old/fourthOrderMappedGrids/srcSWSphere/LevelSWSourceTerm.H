#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELSWSOURCETERM_H_
#define _LEVELSWSOURCETERM_H_

#include "LevelSourceTerm.H"
#include "MOLShallowWaterPhysics.H"

#include "NamespaceHeader.H"

class LevelSWSourceTerm : public LevelSourceTerm
{
public:

  /// Default constructor
  LevelSWSourceTerm();

  /// Destructor
  virtual ~LevelSWSourceTerm();

  // For now, we'll assume constant bottom topography.
  // But eventually, we'll want to include nonzero Psi_B.
  // Then we'll have analytic z(alpha, beta) as ScalarFunction
  // or analytic dz/dalpha, dz/dbeta as VectorFunction.
  // (alpha, beta in range -pi/4:pi/4)

  /// Define the object so that time stepping can begin (actual constructor)
  virtual void define(MultiBlockCoordSys* const a_mblexPtr,
                      const MOLPhysics* const a_molPhysics,
                      const DisjointBoxLayout& a_grids);

  /// Object factory for this class.  Pure virtual function.
  /**
   */
  virtual LevelSourceTerm* new_sourceTerm() const;

  /// Increment a_rhs by source term, using a_U.  Pure virtual function.
  virtual void addSourceTerm(LevelData<FArrayBox>&   a_rhs,
                             LevelData<FArrayBox>&   a_U);

protected:

private:

  // Disallowed for all the usual reasons
  void operator=(const LevelSWSourceTerm&);
  LevelSWSourceTerm(const LevelSWSourceTerm&);
};

#include "NamespaceFooter.H"

#endif
