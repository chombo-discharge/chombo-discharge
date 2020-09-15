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
#include "AMRLevelCons.H"
#include "AMRLevelConsFactory.H"

#include "NamespaceHeader.H"

AMRLevelConsFactory::AMRLevelConsFactory()
{
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////

// Virtual constructor
AMRLevel* AMRLevelConsFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelCons
  AMRLevelCons* amrConsPtr = new AMRLevelCons();

  // Set up new object
  transferSettings(amrConsPtr);
  // Return it
  return (static_cast <AMRLevel*> (amrConsPtr));
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::transferSettings(AMRLevelCons* a_newPtr) const
{
  CH_assert(a_newPtr != NULL);

  a_newPtr->CFL(m_cfl);
  a_newPtr->spaceOrder(m_spaceOrder);
  a_newPtr->limitFaceValues(m_limitFaceValues);
  a_newPtr->highOrderLimiter(m_highOrderLimiter);
  a_newPtr->useFlattening(m_useFlattening);
  a_newPtr->initialAverage(m_initialAverage);
  a_newPtr->noPPM(m_noPPM);
  a_newPtr->doDeconvolution(m_doDeconvolution);
  a_newPtr->doFaceDeconvolution(m_doFaceDeconvolution);
  a_newPtr->useArtificialViscosity(m_useArtificialViscosity);
  a_newPtr->artificialViscosity(m_artificialViscosity);
  a_newPtr->useArtVisc(m_useArtVisc);
  a_newPtr->ratioArtVisc(m_ratioArtVisc);
  a_newPtr->forwardEuler(m_forwardEuler);
  a_newPtr->enforceMinVal(m_enforceMinVal, m_minVal);
  a_newPtr->domainLength(m_domainLength);
  a_newPtr->refinementThreshold(m_refineThresh);
  a_newPtr->refinementIsScaled(m_refinementIsScaled);
  a_newPtr->tagPressure(m_tagPressure);
  a_newPtr->tagVorticity(m_tagVorticity);
  a_newPtr->tagBufferSize(m_tagBufferSize);
  a_newPtr->initialDtMultiplier(m_initialDtMultiplier);
  a_newPtr->verbosity(m_verbosity);
  a_newPtr->molPhysics(m_molPhysics);
}

//////////////////////////////////////////////////////////////////////////////

AMRLevelConsFactory::~AMRLevelConsFactory()
{
  if (m_molPhysics != NULL)
    {
      delete m_molPhysics;
      m_molPhysics = NULL;
    }
}

//////////////////////////////////////////////////////////////////////////////

// CFL number
void AMRLevelConsFactory::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
  m_cflSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// spatial order of accuracy
void AMRLevelConsFactory::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}

//////////////////////////////////////////////////////////////////////////////

// if true, limit face values in advection
// spatial order of accuracy
void AMRLevelConsFactory::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::initialAverage(bool a_initialAverage)
{
  m_initialAverage = a_initialAverage;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////
void AMRLevelConsFactory::useArtVisc(bool a_useArtVisc)
{
  m_useArtVisc = a_useArtVisc;
}


//////////////////////////////////////////////////////////////////////////////
void AMRLevelConsFactory::ratioArtVisc(Real a_ratioArtVisc)
{
  m_ratioArtVisc = a_ratioArtVisc;
}


//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::forwardEuler(bool a_forwardEuler)
{
  m_forwardEuler = a_forwardEuler;
}

//////////////////////////////////////////////////////////////////////////////

/// sets whether to enforce a min value in advection, along with value
void AMRLevelConsFactory::enforceMinVal(bool a_enforceMinVal, Real a_minVal)
{
  m_enforceMinVal = a_enforceMinVal;
  if (m_enforceMinVal) m_minVal = a_minVal;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::verbosity(const int& a_verbosity)
{
  m_verbosity = a_verbosity;
}

//////////////////////////////////////////////////////////////////////////////

// Physical dimension of the longest side of the domain
void AMRLevelConsFactory::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
  m_domainLengthSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// Refinement threshold
void AMRLevelConsFactory::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
  m_refineThreshSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether refinement threshold is scaled with dx
void AMRLevelConsFactory::refinementIsScaled(bool a_refinementIsScaled)
{
  m_refinementIsScaled = a_refinementIsScaled;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether to tag on pressure instead of on density
void AMRLevelConsFactory::tagPressure(bool a_tagPressure)
{
  m_tagPressure = a_tagPressure;
}

//////////////////////////////////////////////////////////////////////////////

// Set whether to tag on vorticity instead of on density
void AMRLevelConsFactory::tagVorticity(bool a_tagVorticity)
{
  m_tagVorticity = a_tagVorticity;
}

//////////////////////////////////////////////////////////////////////////////

// Tag buffer size
void AMRLevelConsFactory::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
  m_tagBufferSizeSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelConsFactory::molPhysics(const MOLPhysics* const a_molPhysics)
{
  m_molPhysics = a_molPhysics->new_molPhysics();
}

//////////////////////////////////////////////////////////////////////////////

// Initial dt multiplier
void AMRLevelConsFactory::initialDtMultiplier(Real a_initialDtMultiplier)
{
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_initialDtMultiplierSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// Check that everything is defined
bool AMRLevelConsFactory::isDefined() const
{
  return (m_cflSet &&
          m_domainLengthSet &&
          m_refineThreshSet &&
          m_tagBufferSizeSet &&
          m_initialDtMultiplierSet);
}

//////////////////////////////////////////////////////////////////////////////

// Some default values
void AMRLevelConsFactory::setDefaultValues()
{
  CFL(0.8);
  spaceOrder(4);
  limitFaceValues(false);
  highOrderLimiter(true);
  initialAverage(false);
  useFlattening(false);
  useArtVisc(false);
  noPPM(false);
  doDeconvolution(true);
  doFaceDeconvolution(true);
  useArtificialViscosity(false);
  artificialViscosity(0.);
  ratioArtVisc(0.);
  forwardEuler(false);
  enforceMinVal(false, -1);
  domainLength(1.0);
  refinementThreshold(0.2);
  refinementIsScaled(false);
  tagPressure(false);
  tagVorticity(false);
  tagBufferSize(2);
  initialDtMultiplier(0.1);
  m_verbosity = 0;
  m_molPhysics = NULL;
}

#include "NamespaceFooter.H"
