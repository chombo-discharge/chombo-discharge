#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 31 Mar 2008

#ifndef _CYLINDRICALSPOKESSECTORMAP_H_
#define _CYLINDRICALSPOKESSECTORMAP_H_

#include "CylindricalSectorMap.H"

#include "NamespaceHeader.H"

/// Class to specify mapping for a cube.
class CylindricalSpokesSectorMap : public CylindricalSectorMap
{
public:

  /// null constructor.  creates identity mapping.
  CylindricalSpokesSectorMap(/// physical coordinates of center of cylinder
                             const RealVect&   a_physCenter,
                             /// length of central cube in each dimension
                             Real              a_bxWidth,
                             /// outer radius of cylinder
                             Real              a_outerRadius,
                             /// which block:  1, 2, 3, 4 for +x, +y, -x, -y
                             int               a_block,
                             /// CELL-centered domain of sector in index space
                             const Box&        a_bxCells);

  ///
  /** Get physical coordinates from mapped coordinates on a box.
      Added by petermc, 27 Mar 2008.
   */
  virtual void setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                  FArrayBox&         a_physFab,
                                  /// box on which to set physical coordinates
                                  const Box&         a_bx,
                                  /// mapped coordinates, SpaceDim components
                                  const FArrayBox&   a_mapFab) const;

  ///
  /** Get mapped (Cartesian) coordinates from physical coordinates on a box.
   */
  virtual void setMapFromPhysical(/// mapped coordinates, SpaceDim components
                                  FArrayBox&         a_mapFab,
                                  /// box on which to set physical coordinates
                                  const Box&         a_bx,
                                  /// physical coordinates, SpaceDim components
                                  const FArrayBox&   a_physFab) const;

  ///
  /** Get Jacobian.
   */
  virtual void setJacobian(/// jacobian on a_bx
                           FArrayBox&         a_jacobianFab,
                           /// box on which to set jacobian
                           const Box&         a_bx,
                           /// mapped coordinates, SpaceDim components
                           const FArrayBox&   a_mapFab) const;

protected:

};

#include "NamespaceFooter.H"

#endif // include guard
