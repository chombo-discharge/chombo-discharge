#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevel.H"
#include "AMRLevelMappedCons.H"
#include "AMRLevelMappedConsFactory.H"

#include "NamespaceHeader.H"

AMRLevelMappedConsFactory::AMRLevelMappedConsFactory()
{
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////

// Virtual constructor
AMRLevel* AMRLevelMappedConsFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelMappedCons
  AMRLevelMappedCons* amrConsPtr = new AMRLevelMappedCons();

  // Set up new object
  transferSettings(amrConsPtr);
  // Return it
  return (static_cast <AMRLevel*> (amrConsPtr));
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::transferSettings(AMRLevelMappedCons* a_newPtr) const
{
  AMRLevelConsFactory::transferSettings(a_newPtr);
  a_newPtr->coordinateSystem(m_coordSysFactPtr);
  a_newPtr->m_plotfile_prefix = m_plotfile_prefix;
  // a_newPtr->m_writeMap = m_writeMap;
  a_newPtr->m_dtFromCells = m_dtFromCells;
  a_newPtr->useSourceTerm(m_useSourceTerm);
  a_newPtr->sourceTerm(m_sourceTerm);
}

//////////////////////////////////////////////////////////////////////////////

AMRLevelMappedConsFactory::~AMRLevelMappedConsFactory()
{
  if (m_sourceTerm != NULL)
    {
      delete m_sourceTerm;
      m_sourceTerm = NULL;
    }
  if (m_coordSysFactPtr != NULL)
    {
      delete m_coordSysFactPtr;
    }
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::dtFromCells(bool a_dtFromCells)
{
  m_dtFromCells = a_dtFromCells;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::coordinateSystemFactory(MultiBlockCoordSysFactory* a_coordSysFact)
{
  m_coordSysFactPtr = a_coordSysFact;
  m_coordSysFactSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::useSourceTerm(bool a_useSourceTerm)
{
  m_useSourceTerm = a_useSourceTerm;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::sourceTerm(const LevelSourceTerm* const a_sourceTerm)
{
  m_sourceTerm = a_sourceTerm->new_sourceTerm();
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::plotPrefix(const std::string& a_plotfile_prefix)
{
  CH_assert(isDefined());

  m_plotfile_prefix = a_plotfile_prefix;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::writeMap(bool a_writeMap)
{
  AMRLevelMappedCons::s_writeMap = a_writeMap;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::writeError(bool a_writeError)
{
  AMRLevelMappedCons::s_writeError = a_writeError;
}

//////////////////////////////////////////////////////////////////////////////

// Check that everything is defined
bool AMRLevelMappedConsFactory::isDefined() const
{
  return (AMRLevelConsFactory::isDefined() && m_coordSysFactSet);
}

//////////////////////////////////////////////////////////////////////////////

// Some default values
void AMRLevelMappedConsFactory::setDefaultValues()
{
  AMRLevelConsFactory::setDefaultValues();
  m_coordSysFactPtr = NULL;
  m_dtFromCells = false;
  m_useSourceTerm = false;
  m_sourceTerm = NULL;
}

#include "NamespaceFooter.H"
