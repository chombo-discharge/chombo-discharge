#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBSTENCIL_H_
#define _MBSTENCIL_H_

#include "MBStencilElement.H"
#include "RefCountedPtr.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe a multi-block interpolation stencil.
class MBStencil
{
public:

  /// null constructor leaves object in undefined state.
  MBStencil();

  ///
  /**
     Constructor.
  */
  MBStencil(RefCountedPtr< Vector<MBStencilElement> >  a_vec);

  /// destructor
  virtual ~MBStencil();

  ///
  /*
    Same as constructor.
   */
  virtual void define(RefCountedPtr< Vector<MBStencilElement> >  a_vec);

  const MBStencilElement& operator[](int  a_ind) const
  { return (*m_vec)[a_ind]; }

  MBStencilElement& operator[](int  a_ind)
  { return (*m_vec)[a_ind]; }

  int size() const;

protected:

  /// stencil elements
  RefCountedPtr< Vector<MBStencilElement> > m_vec;
};

#include "NamespaceFooter.H"

#endif // include guard
