#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 21 May 2008

#ifndef _SPHERICALSECTORMAP_H_
#define _SPHERICALSECTORMAP_H_

#include "BlockMap.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

/// Class to specify mapping for a cylindrical sector.
class SphericalSectorMap : public BlockMap
{
public:

  /// null constructor.  creates identity mapping.
  SphericalSectorMap(/// physical coordinates of center of cylinder
                     const RealVect&   a_physCenter,
                     /// length of central cube in each dimension
                     Real              a_bxWidth,
                     /// outer radius of cylinder
                     Real              a_outerRadius,
                     /// which block:  0:5 for +x, +y, -x, -y, +z, -z
                     int               a_block,
                     /// CELL-centered domain of sector in index space
                     const Box&        a_bxCells);

  ///
  /** Get physical coordinates from mapped coordinates on a box.
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
  /** Is physical point in this block?  Default is true.
   */
  virtual bool inThisBlock(const RealVect&   a_physPoint) const;

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

  /** physical coordinates of center of sphere
   */
  RealVect m_physCenter;

  /** length of central cube in each dimension
   */
  Real m_bxWidth;

  /** outer radius of sphere
   */
  Real m_outerRadius;

  /** block index:
       0   1   2   3   4   5
      +x  +y  -x  -y  +z  -z
   */
  int m_block;

  /** cartesian coordinates of lower corner of domain
   */
  RealVect m_cartesianOrigin;

  /** cartesian length of domain
   */
  RealVect m_cartesianLength;
};

#include "NamespaceFooter.H"

#endif // include guard
