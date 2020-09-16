#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELVOLTESTFACTORY_H_
#define _AMRLEVELVOLTESTFACTORY_H_

#include <string>

#include "AMRLevelFactory.H"
#include "AMRLevelVolTest.H"
#include "VolTestParameters.H"

#include "NamespaceHeader.H"

/*******************************************************************************
 */
///  Factory for creating AMRLevelVolTest class.
/*
 ******************************************************************************/

class AMRLevelVolTestFactory : public AMRLevelFactory
{
public:

  /// Constructor
  AMRLevelVolTestFactory(NewCoordSysFactory *const a_coordSysFact,
                         const VT::ICType          a_icType,
                         const VT::RefType         a_refType,
                         const int                 a_tagBufferSize,
                         const std::string&        a_plotPrefix);

  /// Destructor
  virtual ~AMRLevelVolTestFactory();

/*--------------------------------------------------------------------*
 * Member functions
 *--------------------------------------------------------------------*/

  /// Create a new AMR level
  virtual AMRLevel* new_amrlevel() const;

/*--------------------------------------------------------------------*
 * Data members
 *--------------------------------------------------------------------*/

private:

  MultiBlockCoordSysFactory *const m_coordSysFact;
  const VT::ICType                 m_icType;
  const VT::RefType                m_refType;
  const int                        m_tagBufferSize;
  const std::string                m_plotPrefix;
};

#include "NamespaceFooter.H"

#endif
