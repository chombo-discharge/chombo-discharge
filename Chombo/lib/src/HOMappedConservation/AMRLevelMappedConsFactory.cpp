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
#include "SingleBlockCSAdaptor.H"

#include "NamespaceHeader.H"

AMRLevelMappedConsFactory::AMRLevelMappedConsFactory()
{
  setDefaultValues();
}

AMRLevelMappedConsFactory::
AMRLevelMappedConsFactory(NewCoordSysFactory *const a_coordSysFact,
                          AMRLevelMappedStabilityStrategy* a_stabilityStrategy,
                          AMRLevelMappedTaggingStrategy* a_taggingStategy)
{
  m_coordSysFact = new SingleBlockCSAdaptorFactory(a_coordSysFact);
  m_stabilityStrategy = a_stabilityStrategy;
  m_taggingStrategy = a_taggingStategy;
  setDefaultValues();
  m_singleBlock = true;
}

AMRLevelMappedConsFactory::
AMRLevelMappedConsFactory(MultiBlockCoordSysFactory *const a_coordSysFact,
                          AMRLevelMappedStabilityStrategy* a_stabilityStrategy,
                          AMRLevelMappedTaggingStrategy* a_taggingStategy)
{
  m_coordSysFact = a_coordSysFact;
  m_stabilityStrategy = a_stabilityStrategy;
  m_taggingStrategy = a_taggingStategy;
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////

// Virtual constructor
AMRLevel* AMRLevelMappedConsFactory::new_amrlevel() const
{
  CH_TIME("AMRLevelMappedConsFactory::new_amrlevel");
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelMappedCons
  AMRLevelMappedCons* amrConsPtr = new AMRLevelMappedCons(m_coordSysFact,
                                                          m_stabilityStrategy,
                                                          m_taggingStrategy);

  // Set up new object
  transferSettings(amrConsPtr);
  // Return it
  return (static_cast <AMRLevel*> (amrConsPtr));
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::transferSettings(AMRLevelCons* a_newPtr) const
{
  CH_TIME("AMRLevelMappedConsFactory::transferSettings");
  // Cast to what we should be dealing with
  AMRLevelMappedCons* amrMappedConsPtr =
    dynamic_cast<AMRLevelMappedCons*>(a_newPtr);
  CH_assert(amrMappedConsPtr != NULL);

  // Transfer what our base class knows about
  AMRLevelConsFactory::transferSettings(amrMappedConsPtr);

  // Now transfer what we know about
  amrMappedConsPtr->coordinateSystem(m_coordSysFact);
  amrMappedConsPtr->m_plotPrefix = m_plotPrefix;
  amrMappedConsPtr->useSourceTerm(m_useSourceTerm);
  amrMappedConsPtr->sourceTerm(m_sourceTerm);
  amrMappedConsPtr->useArtificialDissipation(m_useArtificialDissipation);
  amrMappedConsPtr->artificialDissipation(m_artificialDissipation);
  amrMappedConsPtr->writeJ(m_writeJ);
  amrMappedConsPtr->writeExact(m_writeExact);
  amrMappedConsPtr->writeError(m_writeError);
  amrMappedConsPtr->writeMap(m_writeMap);
}

//////////////////////////////////////////////////////////////////////////////

AMRLevelMappedConsFactory::~AMRLevelMappedConsFactory()
{
  delete m_sourceTerm;
  if (m_singleBlock)
    { // m_coordSysFact was set in constructor that takes NewCoordSysFactory
      delete m_coordSysFact;
    }
  delete m_stabilityStrategy;
  delete m_taggingStrategy;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::coordinateSystemFactory(MultiBlockCoordSysFactory* a_coordSysFact)
{
  m_coordSysFact = a_coordSysFact;
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

void AMRLevelMappedConsFactory::useArtificialDissipation(bool a_useArtificialDissipation)
{
  m_useArtificialDissipation = a_useArtificialDissipation;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::artificialDissipation(Real a_artificialDissipation)
{
  m_artificialDissipation = a_artificialDissipation;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::plotPrefix(const std::string& a_plotPrefix)
{
  CH_assert(isDefined());

  m_plotPrefix = a_plotPrefix;
}

//////////////////////////////////////////////////////////////////////////////

// Check that everything is defined
bool AMRLevelMappedConsFactory::isDefined() const
{
  return (AMRLevelConsFactory::isDefined() && (m_coordSysFact != NULL));
}

//////////////////////////////////////////////////////////////////////////////

// Some default values
void AMRLevelMappedConsFactory::setDefaultValues()
{
  CH_TIME("AMRLevelMappedConsFactory::setDefaultValues");
  AMRLevelConsFactory::setDefaultValues();
  m_useSourceTerm = false;
  m_sourceTerm = NULL;
  m_writeJ = true;
  m_writeExact = false;
  m_writeError = false;
  m_writeMap = true;
  m_useArtificialDissipation = false;
  m_artificialDissipation = 0.;
  m_singleBlock = false;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::
writeJ(bool a_writeJ)
{
  m_writeJ = a_writeJ;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::
writeExact(bool a_writeExact)
{
  m_writeExact = a_writeExact;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::
writeError(bool a_writeError)
{
  m_writeError = a_writeError;
}


//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::
writeMap(bool a_writeMap)
{
  m_writeMap = a_writeMap;
}


#include "NamespaceFooter.H"
