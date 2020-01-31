#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBStencil.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
MBStencil::MBStencil()
{
  // RefCountedPtr does not do this.
  // m_vec = NULL;
}

// ---------------------------------------------------------
MBStencil::MBStencil(RefCountedPtr< Vector<MBStencilElement> >  a_vec)
{
  define(a_vec);
}

// ---------------------------------------------------------
MBStencil::~MBStencil()
{
}

// ---------------------------------------------------------
void
MBStencil::define(RefCountedPtr< Vector<MBStencilElement> >  a_vec)
{
  m_vec = a_vec;
}

// ---------------------------------------------------------
int
MBStencil::size() const
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
