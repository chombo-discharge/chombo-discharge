#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOUNDARYITERATOR_H_
#define _BOUNDARYITERATOR_H_

#include <cstdlib>

#include "Box.H"
#include "REAL.H"
#include "SPACE.H"
#include "IntVect.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

///iterates through the IntVects on the surface of a Box
/**
   BoundaryIterator iterates through the IntVects on the surface of a box.  The actual
   sqeuence of IntVects is implementation-specific.
   Typical usage:

   Box b;
   ...
   BoundaryIterator bit (b);
   for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
     (do operations involving iv)
    }
*/
class BoundaryIterator
{
public:
  ///
  /**
     Default constructor.  This constructs an invalid iterator.
     The user must call define before using.
  */
  BoundaryIterator();

  ///
  /**
     Constructs a BoundaryIterator and associates it with a Box.
     Arguments:
     a_bx (not modified) the Box to iterate over.
  */
  BoundaryIterator(const Box& a_bx);

  void setBox(const Box& a_bx);

  ///
  /**
     Associates a Box with this BoundaryIterator.
     Arguments:
     a_bx (not modified) the Box to iterate over.
  */
  void define(const Box& a_bx);

 
  ///
  ~BoundaryIterator ()
  {
  }

  ///
  /**
     Sets this BoundaryIterator to the first IntVect in its Box.  The
     definition of the "first" IntVect is
     implementation-dependent.
  */
  void begin();

  ///
  /**
     Sets this BoundaryIterator to the first IntVect in its Box.  The
     definition of the "first" IntVect is
     implementation-dependent.
  */
  void reset();

  ///
  /**
     Modifies this BoundaryIterator to set it to the next location in its
     Box.  The definition of the "next location" of a Box is
     implementation-dependent.
  */
  void operator ++ ();

  void next();

  ///
  /**
     Returns the value of the InVect for the current location of this
     BoundaryIterator.
  */
  const IntVect& operator () () const;

  ///
  /**
     Returns true if this BoundaryIterator's location is within its Box.
  */
  bool ok() const;

  ///
  /**
     Can request the current boundary box that this iterator is working on
  */
  const Box& box() const;

  ///
  /** skip to the next boundary box in this iterator
   */
  void nextBox();

protected:
  BoxIterator m_current;
  Box m_box[CH_SPACEDIM*2];
  int m_sdir;
};

inline
BoundaryIterator::BoundaryIterator()
{
  m_sdir=2*CH_SPACEDIM;
}

inline
BoundaryIterator::BoundaryIterator(const Box& a_bx)
{
  define(a_bx);
}

inline 
const Box& BoundaryIterator::box() const
{
  CH_assert(ok());
  return m_box[m_sdir];
}

inline 
void BoundaryIterator::nextBox()
{
  m_sdir++;
  if(ok())
    m_current.define(m_box[m_sdir]);
}

inline
void BoundaryIterator::begin()
{
  m_sdir=0;
  m_current.define(m_box[0]);
}

inline
void BoundaryIterator::reset()
{
  begin();
}

inline
void BoundaryIterator::operator ++ ()
{
  next();
}

inline
void BoundaryIterator::next()
{
  m_current.next();
  if(!m_current.ok())
    {
      m_sdir++;
      if(m_sdir > 2*CH_SPACEDIM) return;
      else
	{
	  m_current.define(m_box[m_sdir]);
	  m_current.begin();
	}
    }
}

inline
const IntVect& BoundaryIterator::operator () () const
{

  return m_current();
}

inline
bool BoundaryIterator::ok() const
{
  return m_sdir<2*CH_SPACEDIM;
}

#include "NamespaceFooter.H"
#endif
