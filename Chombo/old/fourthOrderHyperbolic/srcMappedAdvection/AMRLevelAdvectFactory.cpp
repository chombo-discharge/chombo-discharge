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
#include "AMRLevelAdvect.H"
#include "AMRLevelAdvectFactory.H"

AMRLevelAdvectFactory::AMRLevelAdvectFactory()
{
  setDefaultValues();
}

// Virtual constructor
AMRLevel* AMRLevelAdvectFactory::new_amrlevel() const
{
  // Make sure everything is defined
 CH_assert(isDefined());

  // Create a new AMRLevelAdvect
  AMRLevelAdvect* amrAdvectPtr = new AMRLevelAdvect();

  // Set up new object
  amrAdvectPtr->CFL(m_cfl);
  amrAdvectPtr->spaceOrder(m_spaceOrder);
  amrAdvectPtr->limitFaceValues(m_limitFaceValues);
  amrAdvectPtr->enforceMinVal(m_enforceMinVal, m_minVal);
  amrAdvectPtr->lowOrderFluxScheme(m_lowOrderFluxScheme);
  amrAdvectPtr->redistributeNegativeVal(m_redistributeNegativeVal,
                                        m_maxRedistributionPasses);
  amrAdvectPtr->useHyperviscosity(m_useHyperviscosity);
  amrAdvectPtr->hyperviscosity(m_hyperviscosity);
  amrAdvectPtr->domainLength(m_domainLength);
  amrAdvectPtr->refinementThreshold(m_refineThresh);
  amrAdvectPtr->tagBufferSize(m_tagBufferSize);
  amrAdvectPtr->initialDtMultiplier(m_initialDtMultiplier);
  amrAdvectPtr->verbosity(m_verbosity);
  amrAdvectPtr->IBC(m_advect_ibc);
  amrAdvectPtr->coordinateSystem(m_coordSysFactPtr);
  amrAdvectPtr->m_plotfile_prefix = m_plotfile_prefix;
  // Return it
  return (static_cast <AMRLevel*> (amrAdvectPtr));
}

AMRLevelAdvectFactory::~AMRLevelAdvectFactory()
{
   delete m_coordSysFactPtr;
}

// CFL number
void AMRLevelAdvectFactory::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
  m_cflSet = true;
}

// spatial order of accuracy
void AMRLevelAdvectFactory::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}

// if true, limit face values in advection
// spatial order of accuracy
void AMRLevelAdvectFactory::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

  /// sets whether to enforce a min value in advection, along with value
void AMRLevelAdvectFactory::enforceMinVal(bool a_enforceMinVal, Real a_minVal)
{
  m_enforceMinVal = a_enforceMinVal;
  if (m_enforceMinVal)
  {
     m_minVal = a_minVal;
  }
}

void AMRLevelAdvectFactory::lowOrderFluxScheme(std::string a_lowOrderFluxScheme)
{
   m_lowOrderFluxScheme = a_lowOrderFluxScheme;
}

void AMRLevelAdvectFactory::redistributeNegativeVal(
   bool a_redistributeNegativeVal, int a_maxRedistributionPasses)
{
   m_redistributeNegativeVal = a_redistributeNegativeVal;
   m_maxRedistributionPasses = a_maxRedistributionPasses;
}

void AMRLevelAdvectFactory::useHyperviscosity(bool a_useHyperviscosity)
{
   m_useHyperviscosity = a_useHyperviscosity;
}

void AMRLevelAdvectFactory::hyperviscosity(Real a_hyperviscosity)
{
   CH_assert(a_hyperviscosity>=0);
   m_hyperviscosity = a_hyperviscosity;
}

void AMRLevelAdvectFactory::verbosity(const int& a_verbosity)
{
  m_verbosity = a_verbosity;
}

// Physical dimension of the longest side of the domain
void AMRLevelAdvectFactory::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
  m_domainLengthSet = true;
}

// Refinement threshold
void AMRLevelAdvectFactory::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
  m_refineThreshSet = true;
}

// Tag buffer size
void AMRLevelAdvectFactory::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
  m_tagBufferSizeSet = true;
}

// Initial dt multiplier
void AMRLevelAdvectFactory::initialDtMultiplier(Real a_initialDtMultiplier)
{
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_initialDtMultiplierSet = true;
}

// Initial dt multiplier
void AMRLevelAdvectFactory::coordinateSystemFactory(CoordSysFactory<FArrayBox,FluxBox>* a_coordSysFact)
{
  m_coordSysFactPtr = a_coordSysFact;
  m_coordSysFactSet = true;
}

void AMRLevelAdvectFactory::plotPrefix(const std::string& a_plotfile_prefix)
{
  CH_assert(isDefined());

  m_plotfile_prefix = a_plotfile_prefix;
}

// Check that everything is defined
bool AMRLevelAdvectFactory::isDefined() const
{
  return (m_cflSet &&
          m_domainLengthSet &&
          m_refineThreshSet &&
          m_tagBufferSizeSet &&
          m_initialDtMultiplierSet &&
          m_coordSysFactSet);
}

// Some default values
void AMRLevelAdvectFactory::setDefaultValues()
{
  CFL(0.8);
  spaceOrder(4);
  limitFaceValues(false);
  m_enforceMinVal = 0;
  m_minVal = -100000000.0;
  m_lowOrderFluxScheme = "CTU";
  m_redistributeNegativeVal = false;
  m_maxRedistributionPasses = 0;
  useHyperviscosity(false);
  hyperviscosity(0.1);
  domainLength(1.0);
  refinementThreshold(0.2);
  tagBufferSize(2);
  initialDtMultiplier(0.1);
  m_verbosity = 0;
  m_coordSysFactPtr = NULL;
}
