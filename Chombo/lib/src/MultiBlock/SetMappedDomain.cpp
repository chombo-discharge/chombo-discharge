#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SetMappedDomain.H"
// #include "UGIO.H"
// #include "CH_HDF5.H"

#include "NamespaceHeader.H"


//=====================================================================

SetMappedDomain::SetMappedDomain ()
{
  m_isDefined = false;
}


// ---------------------------------------------------------
SetMappedDomain::~SetMappedDomain ()
{
  if (isDefined())
    {
      for (int iblock = 0; iblock < m_nblocks; iblock++)
        {
          delete m_blockMaps[iblock];
        }
    }
}


// ---------------------------------------------------------
void SetMappedDomain::defineMappedDomain()
{
  m_mappedDomain.define(m_mappedBlocks);

  m_isDefined = true;
}


// ---------------------------------------------------------
bool SetMappedDomain::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
const MappedDomain& SetMappedDomain::mappedDomain() const
{
  CH_assert(isDefined());
  return m_mappedDomain;
}

#include "NamespaceFooter.H"
