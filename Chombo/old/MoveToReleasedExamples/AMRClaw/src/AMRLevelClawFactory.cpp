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
#include "AMRLevelClaw.H"
#include "AMRLevelClawFactory.H"

AMRLevelClawFactory::AMRLevelClawFactory()
{
  setDefaultValues();
}

// Virtual constructor
AMRLevel* AMRLevelClawFactory::new_amrlevel() const
{
  // Make sure everything is defined
 CH_assert(isDefined());

  // Create a new AMRLevelClaw
  AMRLevelClaw* amrClawPtr = new AMRLevelClaw();

  // Set up new object
  amrClawPtr->CFL(m_cfl);
  amrClawPtr->domainLength(m_domainLength);
  amrClawPtr->refinementThreshold(m_refineThresh);
  amrClawPtr->tagBufferSize(m_tagBufferSize);
  amrClawPtr->initialDtMultiplier(m_initialDtMultiplier);
  amrClawPtr->clawPatch(m_clawPatch);
  amrClawPtr->verbosity(m_verbosity);

  // Return it
  return (static_cast <AMRLevel*> (amrClawPtr));
}

AMRLevelClawFactory::~AMRLevelClawFactory()
{
}

// CFL number
void AMRLevelClawFactory::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
  m_cflSet = true;
}

void AMRLevelClawFactory::verbosity(const int& a_verbosity)
{
  m_verbosity = a_verbosity;
}

// Physical dimension of the longest side of the domain
void AMRLevelClawFactory::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
  m_domainLengthSet = true;
}

// Refinement threshold
void AMRLevelClawFactory::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
  m_refineThreshSet = true;
}

// Tag buffer size
void AMRLevelClawFactory::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
  m_tagBufferSizeSet = true;
}

// Initial dt multiplier
void AMRLevelClawFactory::initialDtMultiplier(Real a_initialDtMultiplier)
{
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_initialDtMultiplierSet = true;
}

// ClawPatch object (used as a factory)
void AMRLevelClawFactory::clawPatch(const ClawPatch& a_clawPatch)
{
  m_clawPatch.define(a_clawPatch);
  m_clawPatchSet = true;
}

void AMRLevelClawFactory::set_stateNames(std::vector<string> stateNames)
{
  m_stateNames = stateNames;
}

// Check that everything is defined
bool AMRLevelClawFactory::isDefined() const
{
  return (m_cflSet &&
          m_domainLengthSet &&
          m_refineThreshSet &&
          m_tagBufferSizeSet &&
          m_initialDtMultiplierSet &&
          m_clawPatchSet);
}

// Some default values
void AMRLevelClawFactory::setDefaultValues()
{
  CFL(0.8);
  domainLength(1.0);
  refinementThreshold(0.2);
  tagBufferSize(2);
  initialDtMultiplier(0.1);
  m_verbosity = 0;
}
