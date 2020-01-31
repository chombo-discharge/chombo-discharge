#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// bvs  06/06/06

#ifndef _BLOCKMAP_H_
#define _BLOCKMAP_H_


#include <iostream>
#include "RealVect.H"
#include "Box.H"
#include "Misc.H"
#include "SPACE.H"
#include <map>
#include <string>

#include "NamespaceHeader.H"

class FArrayBox;

// ---------------------------------------------------------

/// Class to specify mapping for a given block.
class BlockMap
{
public:

  /// null constructor.  creates identity mapping.
  BlockMap()
  {
    m_refinement = 1;
  }

  virtual ~BlockMap()
  {
  }

  /// Compute mapped coordinates
  /**
     Given a CELL-centered box, compute and return the
     physical coordinates of the corners of all the cells in a NODE-centered
     FArrayBox.

     The definition of this is supplied by the specific mapped grid object
     derived from this base class.

     the default implementation is the identity mapping.
  */
  virtual void getCornerCoordinates(/// NODE-centered, physical coordinates of NODEs surrounding the given CELLs
                                    FArrayBox&   a_cornerCoords,
                                    /// CELL-centered box
                                    const Box&   a_box) const ;

  ///
  /**
     Given a CELL-centered box, compute and return the
     physical coordinates of the centers of all the cells in a CELL-centered
     FArrayBox.

     The definition of this is supplied by the specific mapped grid object
     derived from this base class.

     the default implementation is the identity mapping.
  */
  virtual void getCellCenterCoordinates(/// CELL-centered, physical coordinates of centers of CELLs
                                        FArrayBox&   a_cellCenterCoords,
                                        /// CELL-centered box
                                        const Box&   a_box) const ;

  ///
  virtual void cartesianToReal(RealVect& a_point) const ;

  ///
  virtual void realToCartesian(RealVect& a_point) const ;

  ///
  /** Get physical coordinates from mapped coordinates on a box.
      Added by petermc, 27 Mar 2008.
   */
  virtual void setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                  FArrayBox&         a_physFab,
                                  /// box on which to set physical coordinates
                                  const Box&         a_bx,
                                  /// cartesian coordinates, SpaceDim components
                                  const FArrayBox&   a_mapFab) const ;

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
  /** Get Jacobian.  Default is zero.
   */
  virtual void setJacobian(/// jacobian on a_bx
                           FArrayBox&         a_jacobianFab,
                           /// box on which to set jacobian
                           const Box&         a_bx,
                           /// mapped coordinates, SpaceDim components
                           const FArrayBox&   a_mapFab) const;

  ///
  /** Is physical point in this block?  Default is true.
   */
  virtual bool inThisBlock(const RealVect&   a_physPoint) const;

  ///
  /** Set refinement ratio.  Default is 1.
   */
  void setRefinement(int a_refinement)
  {
    m_refinement = a_refinement;
  }

  //===========================================================================
  // I/O functions.  Provides a means for the base class to write information
  // which then can be used by the base class to construct the appropriate
  // derived classes.

  /// unique identifier for I/O
  virtual std::string name() const;

  virtual void write(std::ostream& os) const;

  virtual BlockMap* read(std::istream& is) const;

  static BlockMap* factory(std::istream& is);

  //  static std::map<std::string, BlockMap*> MapConstructors;

  static BlockMap Identity;

  /// refinement ratio
  int m_refinement;
};

#include "NamespaceFooter.H"

#endif // include guard
