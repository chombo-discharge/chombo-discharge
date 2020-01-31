#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBVECTORSTENCIL_H_
#define _MBVECTORSTENCIL_H_

#include "MBVectorStencilElement.H"
#include "RefCountedPtr.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe a multi-block interpolation stencil.
class MBVectorStencil
{
public:

  /// null constructor leaves object in undefined state.
  MBVectorStencil();

  ///
  /**
     Constructor.
  */
  MBVectorStencil(RefCountedPtr< Vector<MBVectorStencilElement> >  a_vec);

  /// destructor
  virtual ~MBVectorStencil();

  ///
  /*
    Same as constructor.
   */
  virtual void define(RefCountedPtr< Vector<MBVectorStencilElement> >  a_vec);

  const MBVectorStencilElement& operator[](int  a_ind) const
  { return (*m_vec)[a_ind]; }

  MBVectorStencilElement& operator[](int  a_ind)
  { return (*m_vec)[a_ind]; }

  int size() const;

protected:

  /// stencil elements
  RefCountedPtr< Vector<MBVectorStencilElement> > m_vec;
};

#include "NamespaceFooter.H"

#endif // include guard
