#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBSTENCILITERATOR_H_
#define _MBSTENCILITERATOR_H_

#include <cstdlib>

#include "MBStencil.H"

#include "NamespaceHeader.H"

///iterates through the MBStencilElements of a MBStencil.
/**
   MBStencilIterator iterates through the MBStencilElements of a MBStencil.

   MBStencil stencil;
   ...
   MBStencilIterator sit(stencil);
   for (sit.begin(); sit.ok(); ++sit)
    {
      const MBStencilElement& stencilElement = sit();
      (do operations involving stencilElement)
    }
*/
class MBStencilIterator
{
public:
  ///
  /**
     Default constructor.  This constructs an invalid iterator.
     The user must call define before using.
  */
  MBStencilIterator();

  ///
  /**
     Constructs a MBStencilIterator and associates it with a Box.
     Arguments:
     a_bx (not modified) the Box to iterate over.
  */
  MBStencilIterator(const MBStencil& a_stencil);

  void setStencil(const MBStencil& a_stencil);

  ///
  /**
     Associates a MBStencil with this MBStencilIterator.
     Arguments:
     a_stencil (not modified) the MBStencil to iterate over.
  */
  void define(const MBStencil& a_stencil);

  ///
  /**
     Copy constructor.
     Arguments:
     a_iterIn (not modified) the MBStencilIterator to copy.
  */
  MBStencilIterator(const MBStencilIterator& a_iterIn);

  ///
  ~MBStencilIterator ()
  {
  }

  ///
  /**
     Sets this MBStencilIterator to first MBStencilElement in its MBStencil.
  */
  void begin();

  ///
  /**
     Sets this MBStencilIterator to first MBStencilElement in its MBStencil.
  */
  void reset();

  ///
  /**
     Modifies this MBStencilIterator to set it to
     the next MBStencilElement in its MBStencil.
  */
  void operator ++ ();

  void next();

  ///
  /**
     Returns the value of the MBStencilElement for
     the current location of this MBStencilIterator.
  */
  const MBStencilElement& operator () () const;

  ///
  /**
     Returns true if this MBStencilIterator's location is within its MBStencil.
  */
  bool ok();

protected:

  MBStencil m_stencil;
  int m_indCurrent;
  int m_indLo;
  int m_indHi;
};

inline
MBStencilIterator::MBStencilIterator()
{
  m_indCurrent = -1;
  m_indLo = 0;
  m_indHi = -1;
}

inline
MBStencilIterator::MBStencilIterator(const MBStencil& a_stencil)
{
  define(a_stencil);
}

inline
MBStencilIterator::MBStencilIterator(const MBStencilIterator& a_iterIn)
{
  m_stencil = a_iterIn.m_stencil;
  m_indCurrent = a_iterIn.m_indCurrent;
  m_indLo =   a_iterIn.m_indLo;
  m_indHi =   a_iterIn.m_indHi;
}

inline
void MBStencilIterator::begin()
{
  if (m_indLo <= m_indHi) m_indCurrent = m_indLo;
}

inline
void MBStencilIterator::reset()
{
  begin();
}

inline
void MBStencilIterator::operator ++ ()
{
  next();
}

inline
void MBStencilIterator::next()
{
  m_indCurrent++;
}

inline
const MBStencilElement& MBStencilIterator::operator () () const
{
  CH_assert(m_indCurrent <= m_indHi);
  CH_assert(m_indCurrent >= m_indLo);
  return m_stencil[m_indCurrent];
}

inline
bool MBStencilIterator::ok()
{
  return (m_indCurrent <= m_indHi);
}

#include "NamespaceFooter.H"
#endif
