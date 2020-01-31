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
#include "AMRLevelSelfGravity.H"
#include "AMRLevelSelfGravityFactory.H"

AMRLevelSelfGravityFactory::AMRLevelSelfGravityFactory()
{
  m_godunovPhysics = NULL;
  m_refCellTagger  = NULL;
  m_isDefined = false;
}

void AMRLevelSelfGravityFactory::define(const Real&                 a_cfl,
                             const Real&                 a_domainLength,
                             const int&                  a_verbosity,
                             const int&                  a_tagBufferSize,
                             const int&                  a_maxInitRefLevel,
                             const Real&                 a_initialDtMultiplier,
                             const GodunovPhysics* const a_godunovPhysics,
                             const int&                  a_normalPredOrder,
                             const bool&                 a_useFourthOrderSlopes,
                             const bool&                 a_usePrimLimiting,
                             const bool&                 a_useCharLimiting,
                             const bool&                 a_useFlattening,
                             const bool&                 a_useArtificialViscosity,
                             const Real&                 a_artificialViscosity,
                             const RefCellTagger* const  a_refCellTagger,
                             const bool&                 a_useDeltaPhiCorr,
                             const StencilType&          a_stencil)
{
  // Set the CFL number
  m_cfl = a_cfl;

  // Set the physical dimension of the longest side of the domain
  m_domainLength = a_domainLength;

  // Store the verbosity of the object
  m_verbosity = a_verbosity;

  // Set the tag buffer size
  m_tagBufferSize = a_tagBufferSize;

  // Store the initial dt multiplier
  m_initialDtMultiplier = a_initialDtMultiplier;

  if (m_godunovPhysics != NULL)
  {
    delete m_godunovPhysics;
    m_godunovPhysics = NULL;
  }
  m_godunovPhysics = a_godunovPhysics->new_godunovPhysics();

  m_normalPredOrder = a_normalPredOrder;

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_usePrimLimiting      = a_usePrimLimiting;
  m_useCharLimiting      = a_useCharLimiting;
  m_useFlattening        = a_useFlattening;

  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity    = a_artificialViscosity;

  // create RefCellTagger object for tagging cells to be refined
  if (m_refCellTagger != NULL)
  {
    delete m_refCellTagger;
    m_refCellTagger = NULL;
  }
  m_refCellTagger = a_refCellTagger->new_refCellTagger();

  m_maxInitRefLevel = a_maxInitRefLevel;

#ifdef GRAVITY
  // Set the stencil for differentiating the grav. potential:
  m_stencil = a_stencil;

  // whether or not the deltaPhi correction should be applie
  m_useDeltaPhiCorr = a_useDeltaPhiCorr;
#endif

  m_isDefined = true;
}

// Virtual constructor
AMRLevel* AMRLevelSelfGravityFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelSelfGravity
  AMRLevelSelfGravity* amrGodPtr = new AMRLevelSelfGravity();

  // Set up new object
  amrGodPtr->defineParams(m_cfl,
                          m_domainLength,
                          m_verbosity,
                          m_tagBufferSize,
                          m_maxInitRefLevel,
                          m_initialDtMultiplier,
                          m_godunovPhysics,
                          m_normalPredOrder,
                          m_useFourthOrderSlopes,
                          m_usePrimLimiting,
                          m_useCharLimiting,
                          m_useFlattening,
                          m_useArtificialViscosity,
                          m_artificialViscosity,
                          m_refCellTagger,
                          m_useDeltaPhiCorr,
                          m_stencil);

  // Return it
  return (static_cast <AMRLevel*>(amrGodPtr));
}

AMRLevelSelfGravityFactory::~AMRLevelSelfGravityFactory()
{
  if (m_godunovPhysics != NULL)
  {
    delete m_godunovPhysics;
  }

  if (m_refCellTagger != NULL)
  {
    delete m_refCellTagger;
  }
}

// Check that everything is defined
bool AMRLevelSelfGravityFactory::isDefined() const
{
  return m_isDefined;
}

