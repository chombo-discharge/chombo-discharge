#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CODIMCOPYMANAGERI_H_
#define _CODIMCOPYMANAGERI_H_

#include "CoDimCopier.H"
#include "CornerCopier.H"

#include "NamespaceHeader.H"

/// defining constructor -- calls matching define function
template<class T>
CoDimCopyManager<T>::CoDimCopyManager(const DisjointBoxLayout& a_grids, const BoxLayout& a_dest, const IntVect& a_ghostVect, bool a_exchange, int a_maxCoDim)
{
  define(a_grids, a_dest, a_ghostVect, a_exchange, a_maxCoDim);
}


/// destructor
template<class T>
CoDimCopyManager<T>::~CoDimCopyManager()
{
  /// clean up memory
  for (int n=0; n<m_copierVect.size(); n++)
    {
      if (m_copierVect[n] != NULL)
        {
          delete m_copierVect[n];
          m_copierVect[n] = NULL;
        }
    }
}
  

/// define
template<class T>
void
CoDimCopyManager<T>::define(const DisjointBoxLayout& a_grids, const BoxLayout& a_dest, const IntVect& a_ghostVect, bool a_exchange, int a_maxCoDim)
{
  m_ghostVect = a_ghostVect;
  m_srcGrids = a_grids;
  m_destGrids = a_dest;

  CH_assert(m_copierVect.size() == 0);
  m_copierVect.resize(a_maxCoDim-1, NULL);
  for (int n=0; n<a_maxCoDim-1; n++)
    {
      int coDim = n+2;      
      // define a CoDimCopier here
      CoDimCopier* copyPtr = new CoDimCopier(a_grids, a_dest, a_grids.physDomain(),
                                             a_ghostVect, coDim, a_exchange);
      m_copierVect[n] = static_cast<Copier*>(copyPtr);
    } // end loop over codims
}


/// orchestrates CoDim exchanges
template<class T>
void
CoDimCopyManager<T>::manageExchanges(LevelData<T>& a_data) const
{
  for (int n=0; n<m_copierVect.size(); n++)
    {
      a_data.exchange(*m_copierVect[n]);
    }
}


/// orchestrates CoDim copies
template<class T>
void
CoDimCopyManager<T>::manageCopies(const LevelData<T>& a_data, BoxLayoutData<T>& a_dest) const
{
  for (int n=0; n<m_copierVect.size(); n++)
    {
      a_data.copyTo(a_dest, *m_copierVect[n]);
    }
}


#include "NamespaceFooter.H"

#endif


