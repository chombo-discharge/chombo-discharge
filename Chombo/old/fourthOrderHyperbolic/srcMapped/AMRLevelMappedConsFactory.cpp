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
  amrConsPtr->CFL(m_cfl);
  amrConsPtr->spaceOrder(m_spaceOrder);
  amrConsPtr->limitFaceValues(m_limitFaceValues);
  amrConsPtr->useFlattening(m_useFlattening);
  amrConsPtr->initialAverage(m_initialAverage);
  amrConsPtr->noPPM(m_noPPM);
  amrConsPtr->doDeconvolution(m_doDeconvolution);
  amrConsPtr->doFaceDeconvolution(m_doFaceDeconvolution);
  amrConsPtr->useArtificialViscosity(m_useArtificialViscosity);
  amrConsPtr->artificialViscosity(m_artificialViscosity);
  amrConsPtr->useArtVisc(m_useArtVisc);
  amrConsPtr->ratioArtVisc(m_ratioArtVisc);
  amrConsPtr->forwardEuler(m_forwardEuler);
  amrConsPtr->enforceMinVal(m_enforceMinVal, m_minVal);
  amrConsPtr->domainLength(m_domainLength);
  amrConsPtr->refinementThreshold(m_refineThresh);
  amrConsPtr->tagBufferSize(m_tagBufferSize);
  amrConsPtr->initialDtMultiplier(m_initialDtMultiplier);
  amrConsPtr->verbosity(m_verbosity);
  // amrConsPtr->IBC(m_cons_ibc);
  amrConsPtr->godunovPhysics(m_gdnvPhysics);
  amrConsPtr->coordinateSystem(m_coordSysFactPtr);
  amrConsPtr->m_plotfile_prefix = m_plotfile_prefix;
  amrConsPtr->m_dtFromCells = m_dtFromCells;
  // Return it
  return (static_cast <AMRLevel*> (amrConsPtr));
}

//////////////////////////////////////////////////////////////////////////////

AMRLevelMappedConsFactory::~AMRLevelMappedConsFactory()
{
  delete m_coordSysFactPtr;
  if (m_gdnvPhysics != NULL)
    {
      delete m_gdnvPhysics;
      m_gdnvPhysics = NULL;
    }
}

//////////////////////////////////////////////////////////////////////////////

// CFL number
void AMRLevelMappedConsFactory::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
  m_cflSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// spatial order of accuracy
void AMRLevelMappedConsFactory::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}

//////////////////////////////////////////////////////////////////////////////

// if true, limit face values in advection
// spatial order of accuracy
void AMRLevelMappedConsFactory::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::initialAverage(bool a_initialAverage)
{
  m_initialAverage = a_initialAverage;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::useFlattening(bool a_useFlattening)
{
  m_useFlattening = a_useFlattening;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::noPPM(bool a_noPPM)
{
  m_noPPM = a_noPPM;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::doDeconvolution(bool a_doDeconvolution)
{
  m_doDeconvolution = a_doDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::doFaceDeconvolution(bool a_doFaceDeconvolution)
{
  m_doFaceDeconvolution = a_doFaceDeconvolution;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::useArtificialViscosity(bool a_useArtificialViscosity)
{
  m_useArtificialViscosity = a_useArtificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::artificialViscosity(Real a_artificialViscosity)
{
  m_artificialViscosity = a_artificialViscosity;
}

//////////////////////////////////////////////////////////////////////////////
void AMRLevelMappedConsFactory::useArtVisc(bool a_useArtVisc)
{
  m_useArtVisc = a_useArtVisc;
}


//////////////////////////////////////////////////////////////////////////////
void AMRLevelMappedConsFactory::ratioArtVisc(Real a_ratioArtVisc)
{
  m_ratioArtVisc = a_ratioArtVisc;
}


//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::forwardEuler(bool a_forwardEuler)
{
  m_forwardEuler = a_forwardEuler;
}

//////////////////////////////////////////////////////////////////////////////

/// sets whether to enforce a min value in advection, along with value
void AMRLevelMappedConsFactory::enforceMinVal(bool a_enforceMinVal, Real a_minVal)
{
  m_enforceMinVal = a_enforceMinVal;
  if (m_enforceMinVal) m_minVal = a_minVal;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::verbosity(const int& a_verbosity)
{
  m_verbosity = a_verbosity;
}

//////////////////////////////////////////////////////////////////////////////

// Physical dimension of the longest side of the domain
void AMRLevelMappedConsFactory::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
  m_domainLengthSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// Refinement threshold
void AMRLevelMappedConsFactory::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
  m_refineThreshSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// Tag buffer size
void AMRLevelMappedConsFactory::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
  m_tagBufferSizeSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::godunovPhysics(const GodunovPhysics* const a_gdnvPhysics)
{
  m_gdnvPhysics = a_gdnvPhysics->new_godunovPhysics();
}

//////////////////////////////////////////////////////////////////////////////

// Initial dt multiplier
void AMRLevelMappedConsFactory::initialDtMultiplier(Real a_initialDtMultiplier)
{
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_initialDtMultiplierSet = true;
}

//////////////////////////////////////////////////////////////////////////////

// Initial dt multiplier
void AMRLevelMappedConsFactory::dtFromCells(bool a_dtFromCells)
{
  m_dtFromCells = a_dtFromCells;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::coordinateSystemFactory(CoordSysFactory<FArrayBox,FluxBox>* a_coordSysFact)
{
  m_coordSysFactPtr = a_coordSysFact;
  m_coordSysFactSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedConsFactory::plotPrefix(const std::string& a_plotfile_prefix)
{
  CH_assert(isDefined());

  m_plotfile_prefix = a_plotfile_prefix;
}

//////////////////////////////////////////////////////////////////////////////

// Check that everything is defined
bool AMRLevelMappedConsFactory::isDefined() const
{
  return (m_cflSet &&
          m_domainLengthSet &&
          m_refineThreshSet &&
          m_tagBufferSizeSet &&
          m_initialDtMultiplierSet &&
          m_coordSysFactSet);
}

//////////////////////////////////////////////////////////////////////////////

// Some default values
void AMRLevelMappedConsFactory::setDefaultValues()
{
  CFL(0.8);
  spaceOrder(4);
  limitFaceValues(false);
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
  tagBufferSize(2);
  initialDtMultiplier(0.1);
  m_verbosity = 0;
  m_gdnvPhysics = NULL;
  m_coordSysFactPtr = NULL;
  m_dtFromCells = false;
}
