#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFINDEXSPACE_H_
#define _MFINDEXSPACE_H_

#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

class MFIndexSpace
{
public:

  MFIndexSpace();

  ~MFIndexSpace();

  ///
  /**
     a_geoservers contains a GeometryService per fluid
     phase.

      If a_ncellMax is set, that is the max width of
      an internal grid.  Otherwise use defaults
      of (16 in 3D, 64 in 2d)
   */
  void
  define(const Box                     & a_domain,
         const RealVect                & a_origin,
         const Real                    & a_dx,
         const Vector<GeometryService*>& a_geoservers,
         int                             a_nCellMax                              = -1,
         int                             a_maxCoarsenings                        = -1,
         bool                            a_fixOnlyFirstPhaseRegNextToMultiValued = false);

  ///
  void fillEBISLayout(Vector<EBISLayout>& a_ebis,
                      const DisjointBoxLayout& a_grids,
                      const Box& a_domain,
                      const int & nghost) const;

  void fillEBISLayout(EBISLayout& a_ebis,  int phase,
                      const DisjointBoxLayout& a_grids,
                      const Box& a_domain,
                      const int & nghost) const;

  int numPhases() const
  {
    return m_ebis.size();
  }

  const EBIndexSpace*  EBIS(int a_phase) const
  {
    return m_ebis[a_phase];
  }

  ///
  /**
      returns an IntVectSet that is the union of all the irregular
      cells in the entire MFIndexSpace.  Useful for initializing
      a multifluid calculation.
  */
  IntVectSet interfaceRegion(int a_depth=0) const;

private:

  Vector<EBIndexSpace*> m_ebis;

  /// disallowed...
  MFIndexSpace (const MFIndexSpace&);
  MFIndexSpace& operator = (const MFIndexSpace& );

};

#include "NamespaceFooter.H"
#endif
