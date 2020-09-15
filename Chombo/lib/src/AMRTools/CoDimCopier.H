#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CODIMCOPIER_H_
#define _CODIMCOPIER_H_

#include "DisjointBoxLayout.H"
#include "Pool.H"
#include "Vector.H"
#include "ProblemDomain.H"
#include "Copier.H"
#include "LevelData.H"
#include "NamespaceHeader.H"


/// An even strangerer (than CornerCopier) thing to copy from ghost cells to corner ghost cells.
/**
   This is a specialized Copier to fill corner ghost cells of Codim n from codim (n-1) ghost
   cells of adjacent boxes.  Calling exchange with this copier will
   fill _only_ the corner ghost cells.  All other exchange operations
   should be performed by calling a separate exchange with a normal
   Copier. 
   For this to work properly, one needs to call it DIM-1 times, starting with a coDim of 1, 
   and ending with coDim DIM-1. To automate this, we provide a CoDimCopierManager class, which 
   defines the appropriate CoDimCopiers and also managed the series of CoDim exchanges.
*/

class CoDimCopier: public Copier
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  CoDimCopier()
  {
    m_ghost = -IntVect::Unit;
    m_coDim = -1;
  }

  /// Basic defining constructor
  /**  Note that this constructor doesn't take a ghost cell argument,
       so it creates a useless CoDimCopier (no intersections are computed).
       Provided for compatibility with Copier base class.
  */
  CoDimCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
               bool a_exchange = false);

  /// this constructor contains support for periodic BC's
  /**  Note that this constructor doesn't take a ghost cell argument,
       so it creates a useless CoDimCopier (no intersections are computed).
       Provided for compatibility with Copier base class.
  */
  CoDimCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               bool  a_exchange  = false);

  /// Copier to copy into coDim corner ghost cells of a_dest from coDim-1 ghost cells of a_level
  /** Assumption is that a_dest and a_level have the same number of
      ghost cells. coDim is the co-dimensionality of the copier intersections.
  */
  CoDimCopier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const IntVect& a_ghost,
               int   a_coDim,
               bool  a_exchange = false);

  /// Copier to copy data into the coDim corner ghost cells of a_dest from the coDim-1 ghost cells of a_level with support for periodic BCs
  /** Assumption is that a_dest and a_level have the same number of
      ghost cells. coDim is the co-dimensionality of the copier intersections
  */
  CoDimCopier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               const IntVect& a_ghost,
               int   a_coDim,
               bool  a_exchange = false);

  ///
  virtual ~CoDimCopier();

  /// define to match constructor (also doesn't do anything)
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      bool  a_exchange = false);

  /// define to match constructor (also doesn't do anything)
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      bool  a_exchange = false);

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_ghost,
                      int   a_coDim,
                      bool  a_exchange = false);

  /// contains support for periodic BCs
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_ghost,
                      int   a_coDim,
                      bool  a_exchange = false);

  ///
  virtual void clear();

  const IntVect& ghost()
  {
    return m_ghost;
  }

  const int coDim() const
  {
    return m_coDim;
  }

protected:

private:

  IntVect m_ghost;
  int m_coDim;
};

std::ostream& operator<<(std::ostream& os, const CoDimCopier& copier);


#include "NamespaceFooter.H"
#endif
