#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 28 Mar 2008

#ifndef _RECTANGULAREQUALANGLESMAP_H_
#define _RECTANGULAREQUALANGLESMAP_H_

#include "RectangularMap.H"

#include "NamespaceHeader.H"

/// Class to specify mapping for a cube.
class RectangularEqualAnglesMap : public RectangularMap
{
public:

  /// null constructor.  creates identity mapping.
  RectangularEqualAnglesMap(/// physical coordinates of low corner
                            const RealVect&   a_cornerLo,
                            /// physical coordinates of high corner
                            const RealVect&   a_cornerHi,
                            /// CELL-centered domain in index space
                            const Box&        a_bxCells);
  // : RectangularMap(a_cornerLo, a_cornerHi, a_bxCells);

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
