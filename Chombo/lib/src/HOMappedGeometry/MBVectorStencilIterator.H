#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBVECTORSTENCILITERATOR_H_
#define _MBVECTORSTENCILITERATOR_H_

#include <cstdlib>

#include "MBVectorStencil.H"

#include "NamespaceHeader.H"

///iterates through the MBVectorStencilElements of a MBVectorStencil.
/**
   MBVectorStencilIterator iterates through the MBVectorStencilElements of a MBVectorStencil.

   MBVectorStencil stencil;
   ...
   MBVectorStencilIterator sit(stencil);
   for (sit.begin(); sit.ok(); ++sit)
    {
      const MBVectorStencilElement& stencilElement = sit();
      (do operations involving stencilElement)
    }
*/
class MBVectorStencilIterator
{
public:
  ///
  /**
     Default constructor.  This constructs an invalid iterator.
     The user must call define before using.
  */
  MBVectorStencilIterator();

  ///
  /**
     Constructs a MBVectorStencilIterator and associates it with a Box.
     Arguments:
     a_bx (not modified) the Box to iterate over.
  */
  MBVectorStencilIterator(const MBVectorStencil& a_stencil);

  void setStencil(const MBVectorStencil& a_stencil);

  ///
  /**
     Associates a MBVectorStencil with this MBVectorStencilIterator.
     Arguments:
     a_stencil (not modified) the MBVectorStencil to iterate over.
  */
  void define(const MBVectorStencil& a_stencil);

  ///
  /**
     Copy constructor.
     Arguments:
     a_iterIn (not modified) the MBVectorStencilIterator to copy.
  */
  MBVectorStencilIterator(const MBVectorStencilIterator& a_iterIn);

  ///
  ~MBVectorStencilIterator ()
  {
  }

  ///
  /**
     Sets this MBVectorStencilIterator to first MBVectorStencilElement in its MBVectorStencil.
  */
  void begin();

  ///
  /**
     Sets this MBVectorStencilIterator to first MBVectorStencilElement in its MBVectorStencil.
  */
  void reset();

  ///
  /**
     Modifies this MBVectorStencilIterator to set it to
     the next MBVectorStencilElement in its MBVectorStencil.
  */
  void operator ++ ();

  void next();

  ///
  /**
     Returns the value of the MBVectorStencilElement for
     the current location of this MBVectorStencilIterator.
  */
  const MBVectorStencilElement& operator () () const;

  ///
  /**
     Returns true if this MBVectorStencilIterator's location is within its MBVectorStencil.
  */
  bool ok();

protected:

  MBVectorStencil m_stencil;
  int m_indCurrent;
  int m_indLo;
  int m_indHi;
};

inline
MBVectorStencilIterator::MBVectorStencilIterator()
{
  m_indCurrent = -1;
  m_indLo = 0;
  m_indHi = -1;
}

inline
MBVectorStencilIterator::MBVectorStencilIterator(const MBVectorStencil& a_stencil)
{
  define(a_stencil);
}

inline
MBVectorStencilIterator::MBVectorStencilIterator(const MBVectorStencilIterator& a_iterIn)
{
  m_stencil = a_iterIn.m_stencil;
  m_indCurrent = a_iterIn.m_indCurrent;
  m_indLo =   a_iterIn.m_indLo;
  m_indHi =   a_iterIn.m_indHi;
}

inline
void MBVectorStencilIterator::begin()
{
  if (m_indLo <= m_indHi) m_indCurrent = m_indLo;
}

inline
void MBVectorStencilIterator::reset()
{
  begin();
}

inline
void MBVectorStencilIterator::operator ++ ()
{
  next();
}

inline
void MBVectorStencilIterator::next()
{
  m_indCurrent++;
}

inline
const MBVectorStencilElement& MBVectorStencilIterator::operator () () const
{
  CH_assert(m_indCurrent <= m_indHi);
  CH_assert(m_indCurrent >= m_indLo);
  return m_stencil[m_indCurrent];
}

inline
bool MBVectorStencilIterator::ok()
{
  return (m_indCurrent <= m_indHi);
}

#include "NamespaceFooter.H"
#endif
