#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BINITEM_H_
#define _BINITEM_H_

#include <iostream>

#include "RealVect.H"
#include "NamespaceHeader.H"

/// Virtual base class for particle data.
/**
  BinItem is a base class to provide the necessary interface
  with BinFab for particle-in-bin/cell problems. Since the only attribute
  of a particle needed by the BinFab is position, it is the only attribute
  present in this class.  All other attributes (charge, velocity, etc)
  must appear in a derived class.
*/

class BinItem
{
public:
  /// Null constructor, copy constructor and operator= can be compiler defined.
  BinItem();

  /// Constructs a \p BinItem at \a a_position
  BinItem(const RealVect& a_position);

  /// Destructor
  virtual ~BinItem();

  /// Same as the constructor
  virtual void define(const RealVect& a_position);

  /// Change the position
  void setPosition(const RealVect& a_position);

  /// Change one component of the position
  void setPosition(const Real a_position, const int a_dimension);

  /// Return a reference to the position
  RealVect& position();

  /// Return a constant reference to the position
  const RealVect& position() const;

  /// Gets one component of the position
  Real position(const int a_dir) const;

  /// Equality operator. Checks whether two BinItems have equal data.
  virtual bool operator==(const BinItem& a_p) const;

  /// Same as the above, but defined for a pointer
  virtual bool operator==(const BinItem* a_p) const;

  /// Inequality operator.
  virtual bool operator!=(const BinItem& a_p) const;

  // Linearization functions

  /// Number of bytes used by \b linearIn()/linearOut()
  /** Returns the size, in number of bytes, of a flat linear
      representation of the data in this object.
  */
  virtual int size() const;

  /// Write a serialized (packed) representation into \a a_buf
  /** Write a linear representation of the internal data. Assumes that
      sufficient memory for the buffer has already been allocated by
      the caller.
  */
  virtual void linearOut(void* a_buf) const;

  /// Extract a serialized (packed) representation from \a a_buf
  /** Given a linear representation of the class data previously made
      using \b linearOut() in \a a_buf, set the data for this class.
  */
  virtual void linearIn(void* a_buf);

  ///
  int preAllocatable()
  {
    return 0;
  }

  // Text output functions

  /// Write a human-readable version to stdout (for debugging)
  virtual void print() const;
  friend std::ostream& operator<<( std::ostream& ostr, const BinItem& a_item );

protected:
  /// this is the position in space of this item
  RealVect m_position;
};

std::ostream& operator<<( std::ostream& ostr, const BinItem& a_item );

#include "NamespaceFooter.H"
#endif
