#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REMOTECOPIER_H_
#define _REMOTECOPIER_H_

#include "Copier.H"
#include <map>
#include "NamespaceHeader.H"

class RemoteCopier: public Copier
{
  public:

  RemoteCopier();
  RemoteCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest, const IntVect& a_ghostVect = IntVect::Zero, bool a_exchange = false);

  virtual ~RemoteCopier();

  void define(const DisjointBoxLayout& a_level,
              const BoxLayout& a_dest,
              const IntVect& a_destGhost,
              bool  a_exchange = false);

  void define(const DisjointBoxLayout& a_level,
              const BoxLayout& a_dest,
              bool a_exchange = false)
  {
    define(a_level, a_dest, IntVect::Zero, a_exchange);
  }

  //! This allows one to define a transfer of information from boxes that are contained by
  //! those in a given disjoint box layout.
  void define(const DisjointBoxLayout& a_from,
              const DisjointBoxLayout& a_to,
              const DisjointBoxLayout& a_sourceBoxes,
              const BoxLayout& a_destBoxes,
              const std::map<Box, Box>& sourceToDestMapping,
              const IntVect& a_destGhost,
              bool a_exchange = false);

  //! This allows one to define a transfer of information from boxes that are contained by
  //! those in a given disjoint box layout.
  void define(const DisjointBoxLayout& a_from,
              const DisjointBoxLayout& a_to,
              const DisjointBoxLayout& a_sourceBoxes,
              const BoxLayout& a_destBoxes,
              const std::map<Box, Box>& sourceToDestMapping,
              bool a_exchange = false)
  {
    define(a_from, a_to, a_sourceBoxes, a_destBoxes, sourceToDestMapping, IntVect::Zero, a_exchange);
  }

  private:

  void define(const DisjointBoxLayout& a_level,
              const BoxLayout& a_dest,
              const ProblemDomain& a_domain,
              bool  a_exchange = false)
  {
  }

  void define(const BoxLayout& a_level,
              const BoxLayout& a_dest,
              const ProblemDomain& a_domain,
              const IntVect& a_destGhost,
              bool  a_exchange = false)
  {
  }

  void exchangeDefine(const DisjointBoxLayout& a_grids,
                      const IntVect& a_ghost)
  {
  }

  void ghostDefine(const DisjointBoxLayout& a_src,
                   const DisjointBoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_srcGhost)
  {
  }

  void reverse()
  {
  }

};

#include "NamespaceFooter.H"
#endif
