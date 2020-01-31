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

#ifndef _RECTANGULARMAP_H_
#define _RECTANGULARMAP_H_

#include "BlockMap.H"
#include "FArrayBox.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

/// Class to specify mapping for a cube.
class RectangularMap : public BlockMap
{
public:

  /// null constructor.  creates identity mapping.
  RectangularMap(/// physical coordinates of low corner
                 const RealVect&   a_cornerLo,
                 /// physical coordinates of high corner
                 const RealVect&   a_cornerHi,
                 /// CELL-centered domain in index space
                 const Box&        a_bxCells); // : BlockMap();

  ///
  /** Get physical coordinates from mapped coordinates on a box.
      Added by petermc, 27 Mar 2008.
   */
  virtual void setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                  FArrayBox&         a_physFab,
                                  /// box on which to set physical coordinates
                                  const Box&         a_bx,
                                  /// mapped coordinates, SpaceDim components
                                  const FArrayBox&   a_mapFab) const = 0;

  ///
  /** Get mapped (Cartesian) coordinates from physical coordinates on a box.
   */
  virtual void setMapFromPhysical(/// mapped coordinates, SpaceDim components
                                  FArrayBox&         a_mapFab,
                                  /// box on which to set physical coordinates
                                  const Box&         a_bx,
                                  /// physical coordinates, SpaceDim components
                                  const FArrayBox&   a_physFab) const = 0;

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
                           const FArrayBox&   a_mapFab) const = 0;

  //===========================================================================
  // I/O functions.  Provides a means for the base class to write information
  // which then can be used by the base class to construct the appropriate
  // derived classes.

  /// unique identifier for I/O
  //  virtual std::string name() const;

  //  virtual void write(std::ostream& os) const;

  //  virtual RectangularMap* read(std::istream& is) const;

  //  static RectangularMap* factory(std::istream& is);

  //  static std::map<std::string, RectangularMap*> MapConstructors;

  //  static RectangularMap Identity;

protected:

  /** physical corners of rectangular domain
   */
  RealVect m_physCornerLo, m_physCornerHi;

  /** physical length of rectangular domain
   */
  RealVect m_bxWidth;

  /** cells of the domain
   */
  Box m_bxCells;

  /** cartesian coordinates of lower corner of domain
   */
  RealVect m_cartesianOrigin;

  /** cartesian length of domain
   */
  RealVect m_cartesianLength;
};

#include "NamespaceFooter.H"

#endif // include guard
