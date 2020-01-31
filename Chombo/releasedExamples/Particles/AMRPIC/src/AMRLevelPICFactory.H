#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELPICFACTORY_H_
#define _AMRLEVELPICFACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelPIC.H"
#include "NamespaceHeader.H"

///
/**
 */
class AMRLevelPICFactory: public AMRLevelFactory
{
public:

  AMRLevelPICFactory();

  /// Constructor
  AMRLevelPICFactory(const ZeldovichIBC* a_ibc,
		     const Real&         a_cfl,
		     const Real&         a_domainLength,
		     const int&          a_maxGridSize,
		     const InterpType&   a_poissonInterpType,
		     const Real&         a_refineThresh);


  /// Destructor
  virtual ~AMRLevelPICFactory();

  /// Virtual constructor
  virtual AMRLevel* new_amrlevel() const;

protected:

  const ZeldovichIBC*       m_ibc;
  Real                      m_cfl;
  Real                      m_domainLength;
  int                       m_maxGridSize;
  InterpType                m_poissonInterpType;
  Real                      m_refineThresh;

private:

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelPICFactory&);
  AMRLevelPICFactory(const AMRLevelPICFactory&);
};

#include "NamespaceFooter.H"

#endif
