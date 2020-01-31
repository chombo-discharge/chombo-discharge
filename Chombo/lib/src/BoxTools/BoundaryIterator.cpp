#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoundaryIterator.H"
#include "NamespaceHeader.H"

void BoundaryIterator::define (const Box& a_bx)
{
  m_sdir=0;
  IntVect lo=a_bx.smallEnd();
  IntVect hi=a_bx.bigEnd();
  for(int i=0; i<CH_SPACEDIM; i++)
    {
      int displace = hi[i]-lo[i];
      hi[i]=lo[i];
      m_box[i]= Box(lo,hi);
      m_box[i+CH_SPACEDIM]=m_box[i];
      m_box[i+CH_SPACEDIM].shift(i,displace);
      lo[i]++;
      hi[i]+=displace;
      hi[i]--;
    }
    
  m_current.define(m_box[0]);
}

void BoundaryIterator::setBox(const Box& a_bx)
{
  define(a_bx);
}

#include "NamespaceFooter.H"
