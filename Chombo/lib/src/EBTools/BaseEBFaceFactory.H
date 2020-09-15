#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _BASEEBFACEFACTORY_H_
#define _BASEEBFACEFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "BaseEBFaceFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"


/// Factory class to produce BaseEBFaceFABs.
/**
   Factory class to produce BaseEBFaceFABs.
   This is needed for LevelData construction.
 */
template <class T>
class BaseEBFaceFactory
  : public DataFactory<BaseEBFaceFAB<T> >
{
public:
  /// factory function.
  /**
      Creates a new baseebFacefab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual BaseEBFaceFAB<T>* create(const Box& a_box, int a_ncomps,
                                   const DataIndex& a_dit) const;

  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  BaseEBFaceFactory(const EBISLayout& a_ebisl, int a_idir);

  ///
  virtual ~BaseEBFaceFactory();

private:
  EBISLayout m_ebisl;
  int m_idir;
  //disallowed constructors and all that
  BaseEBFaceFactory(const BaseEBFaceFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  BaseEBFaceFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const BaseEBFaceFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#include "NamespaceFooter.H"
#include "BaseEBFaceFactoryI.H"

#endif
