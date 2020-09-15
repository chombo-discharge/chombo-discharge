#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CODIMCOPYMANAGER_H_
#define _CODIMCOPYMANAGER_H_

#include "DisjointBoxLayout.H"
#include "Pool.H"
#include "Vector.H"
#include "ProblemDomain.H"
#include "Copier.H"
#include "NamespaceHeader.H"


/// Helper class to manage CoDim exchange operations
/**
   CoDimCopiers are designed to be exercised in order, starting with CoDim=2
   and increasing to SpaceDim. The CoDimCopyManager class manages these operations.
*/
template <class T> class CoDimCopyManager
{
public:
  /// null constructor
  CoDimCopyManager()
  {
  }
  
  /// defining constructor -- calls matching define function
  /** by default, does all coDim copies from CoDim=2->SpaceDim. If you know you only need coDim ghost cells up to a Co-dimension less than SpaceDim, 
      that can be specified, which will make things cheaper.
  */
  CoDimCopyManager(const DisjointBoxLayout& a_grids, const BoxLayout& a_dest, const IntVect& a_ghostVect, bool a_exchange=false, int a_maxCoDim=SpaceDim);
  
  /// destructor
  ~CoDimCopyManager();

  /// define
  void define(const DisjointBoxLayout& a_grids, const BoxLayout& a_dest, const IntVect& a_ghostVect, bool a_exchange=false, int a_maxCoDim=SpaceDim);

  /// orchestrates CoDim exchanges
  void manageExchanges(LevelData<T>& a_data) const;

  /// orchestrates CoDim copies
  void manageCopies(const LevelData<T>& a_data, BoxLayoutData<T>& a_dest) const;
  
  /// returns ghostVect used to define this object (to check for compatibility)
  IntVect ghostVect() const {return m_ghostVect;}

  /// returns reference to src layout used to define this object
  const DisjointBoxLayout& srcGrids() const {return m_srcGrids;}

  /// returns reference to dest layout used to define this object
  const BoxLayout& destGrids() const {return m_destGrids;}  

protected:
  /// Vector of CoDimCopiers
  Vector<Copier*> m_copierVect;

  /// keep ghostVect around just in case
  IntVect m_ghostVect;

  /// keep source and dest layouts around as well
  // src grids used to define this object
  DisjointBoxLayout m_srcGrids;

  // dest grids used to define this object
  BoxLayout m_destGrids;
  
  
private:

  
};

  
#include "NamespaceFooter.H"
#include "CoDimCopyManagerI.H"
  
#endif
