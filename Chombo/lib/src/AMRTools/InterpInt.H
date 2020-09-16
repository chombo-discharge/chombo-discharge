#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INTERPINT_H_
#define _INTERPINT_H_

#include <iostream>
#include <fstream>
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// fills fine LevelData<BaseFab<int> > from a coarse one using piecewise constant interpolation

/**
 */

class InterpInt
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  InterpInt()
    {
      m_isDefined = false;
    }

  ///
  virtual  ~InterpInt(){}

  ///
  InterpInt(const DisjointBoxLayout& a_gridsFine,
            const DisjointBoxLayout& a_gridsCoar,
            int a_ref_ratio,
            int a_num_comps);


  void define(const DisjointBoxLayout& a_gridsFine,
              const DisjointBoxLayout& a_gridsCoar,
              int a_refRat,
              int a_nComp);


  ///
  void
  fillFineFromCoarPWConst(      LevelData<BaseFab<int> >& a_fine,
                          const LevelData<BaseFab<int> >& a_coar,
                          int a_src_comp,
                          int a_dest_comp,
                          int a_num_comp);

protected:

  bool m_isDefined;
  // refinement ratio
  int m_refRat;
  int m_nComp;
  
  DisjointBoxLayout        m_dblCoFi;
  DisjointBoxLayout        m_dblCoar;
  DisjointBoxLayout        m_dblFine;
  LevelData<BaseFab<int> > m_bufCoFi;
};

#include "NamespaceFooter.H"
#endif
