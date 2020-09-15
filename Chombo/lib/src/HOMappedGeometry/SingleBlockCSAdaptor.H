#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SingleBlockCSAdapter.H
// Jeff J., 6 Jun 2011

#ifndef _SINGLEBLOCKCSADAPTOR_H_
#define _SINGLEBLOCKCSADAPTOR_H_

#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

//! \class SingleBlockCSAdapter
//! This class converts a NewCoordSys to a MultiBlockCoordSys so that
//! we can pay less attention to the differences in the two interfaces,
//! and write single- and multi-block algorithms in one framework.
class SingleBlockCSAdaptor: public MultiBlockCoordSys
{
public:

  //! Creates a MultiBlockCoordSys with exactly one block, having the
  //! given coordinate system. This class assumes control of \a a_coordSys.
  SingleBlockCSAdaptor(const NewCoordSys* const a_coordSys,
                       const ProblemDomain& a_problemDomain);

  //! Destructor.
  virtual ~SingleBlockCSAdaptor();

  // Overridden methods.

  void blockRemapping(RealVect& a_xi_valid,
                      int& a_n_valid,
                      const RealVect& a_xiSrc,
                      int a_nSrc) const;
};


//! \class SingleBlockCSAdaptorFactory
//! Factory associated with SingleBlockCSAdaptor.
class SingleBlockCSAdaptorFactory: public MultiBlockCoordSysFactory
{
  public:

  //! Assumes control of a_coordSysFactory.
  explicit SingleBlockCSAdaptorFactory(const NewCoordSysFactory* const a_coordSys);

  ~SingleBlockCSAdaptorFactory();
  MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                  const RealVect& a_dx) const;

  private:

  const NewCoordSysFactory* const m_coordSysFactory;
};

#include "NamespaceFooter.H"

#endif
