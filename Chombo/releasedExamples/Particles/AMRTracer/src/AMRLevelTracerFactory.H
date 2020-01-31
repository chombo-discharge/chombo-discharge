#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELTRACERFACTORY_H_
#define _AMRLEVELTRACERFACTORY_H_

#include "AMRLevelFactory.H"
#include "AdvectionFunctions.H"
#include "AMRLevelTracer.H"
#include "NamespaceHeader.H"

///
/**
 */
class AMRLevelTracerFactory: public AMRLevelFactory
{
public:

  /// Constructor
  AMRLevelTracerFactory(AdvectionVelocityFunction a_advFunc,
			const Real& a_cfl,
			const Real& a_domainLength,
			const int&  a_maxGridSize);


  /// Destructor
  virtual ~AMRLevelTracerFactory()
  {
  }

  /// Virtual constructor
  virtual AMRLevel* new_amrlevel() const;

protected:

  AdvectionVelocityFunction m_advFunc;
  Real                      m_cfl;
  Real                      m_domainLength;
  int                       m_maxGridSize;

private:

  /// Weak construction is disallowed
  AMRLevelTracerFactory()
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelTracerFactory&);
  AMRLevelTracerFactory(const AMRLevelTracerFactory&);
};

#include "NamespaceFooter.H"

#endif
