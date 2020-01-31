#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBVectorStencil.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
MBVectorStencil::MBVectorStencil()
{
  // RefCountedPtr does not do this.
  // m_vec = NULL;
}

// ---------------------------------------------------------
MBVectorStencil::MBVectorStencil(RefCountedPtr< Vector<MBVectorStencilElement> >  a_vec)
{
  define(a_vec);
}

// ---------------------------------------------------------
MBVectorStencil::~MBVectorStencil()
{
}

// ---------------------------------------------------------
void
MBVectorStencil::define(RefCountedPtr< Vector<MBVectorStencilElement> >  a_vec)
{
  m_vec = a_vec;
}

// ---------------------------------------------------------
int
MBVectorStencil::size() const
{
  if ( m_vec.isNull() )
    {
      return 0;
    }
  else
    {
      return m_vec->size();
    }
}

#include "NamespaceFooter.H"
