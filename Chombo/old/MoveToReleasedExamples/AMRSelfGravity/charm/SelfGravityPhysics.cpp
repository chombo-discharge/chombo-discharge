#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "parstream.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"

#include "SelfGravityPhysics.H"

#include "CONSTANTS_C.H"
#include "SelfGravityPhysicsF_F.H"

SelfGravityPhysics::SelfGravityPhysics() : GodunovPhysics()
{
  //m_source = NULL;
}

SelfGravityPhysics::~SelfGravityPhysics()
{
  //m_source = NULL;
}

GodunovPhysics* SelfGravityPhysics::new_godunovPhysics() const
{
  // Make the new object
  SelfGravityPhysics* retval = new SelfGravityPhysics();

  // Pass the current initial and boundary condition (IBC) object
  retval->setPhysIBC(m_bc);
  retval->setSmallPressure(m_smallPressure);

  // Return the new object
  return static_cast<GodunovPhysics*>(retval);
}

// Compute the maximum wave speed
Real SelfGravityPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
                                         const Box&     a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));

  Real speed;

  FORT_MAXWAVESPEEDF(CHF_REAL(speed),
                    CHF_CONST_FRA(a_U),
                    CHF_BOX(a_box));

  return speed;
}

// Compute the maximum wave speed
Real SelfGravityPhysics::getMaxWaveSpeedWithSource(const FArrayBox& a_U,
                                                   const FArrayBox& a_S,
                                                   const Real&      a_dx,
                                                   const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));

  Real speed = 0.0;

  if (!a_S.box().isEmpty())
  {
    FORT_MAXWAVESPEEDWSOURCEF(CHF_REAL(speed),
                              CHF_CONST_FRA(a_U),
                              CHF_CONST_FRA(a_S),
                              CHF_CONST_REAL(a_dx),
                              CHF_BOX(a_box));
  }
  else
  {
    speed = getMaxWaveSpeed(a_U,a_box);
  }

  return speed;
}

// Names of the conserved variables
Vector<string> SelfGravityPhysics::plotNames()
{
  Vector<string> retval;

  retval.push_back("density");
  retval.push_back("x-velocity");
  if (SpaceDim >= 2)
  {
    retval.push_back("y-velocity");

    if (SpaceDim >= 3)
    {
      retval.push_back("z-velocity");
    }
  }
  retval.push_back("pressure");
  retval.push_back("specific-entropy");

  return retval;
}


// Names of the conserved variables
Vector<string> SelfGravityPhysics::stateNames()
{
  Vector<string> retval;

  retval.push_back("density");
  retval.push_back("x-momentum");
  if (SpaceDim >= 2)
  {
    retval.push_back("y-momentum");

    if (SpaceDim >= 3)
    {
      retval.push_back("z-momentum");
    }
  }
  retval.push_back("energy-density");
  retval.push_back("entropy");

  return retval;
}

// Number of conserved variables
int SelfGravityPhysics::numConserved()
{
  return UNUM;
}

// Number of flux variables
int SelfGravityPhysics::numFluxes()
{
  return UFLUX;
}

// Number of primitive variables
int SelfGravityPhysics::numPrimitives()
{
  return WNUM;
}

//  Number of primitive variables for which slopes are computed
int SelfGravityPhysics::numSlopes()
{
  return WSLOPE;
}

bool SelfGravityPhysics::isDefined()
{
  return GodunovPhysics::isDefined();
}

// Compute the primitive variables from the conserved variables
void SelfGravityPhysics::consToPrim(FArrayBox&       a_W,
                                    const FArrayBox& a_U,
                                    const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  FORT_CONSTOPRIMF(CHF_FRA(a_W),
                  CHF_CONST_FRA(a_U),
                  CHF_BOX(a_box));
}

void SelfGravityPhysics::postNormalPred(FArrayBox&       a_dWMinus,
                                        FArrayBox&       a_dWPlus,
                                        const FArrayBox& a_W,
                                        const Real&      a_dt,
                                        const Real&      a_dx,
                                        const int&       a_dir,
                                        const Box&       a_box)
{
  // Bound "a_dWMinus" and "a_dWPlus" so that density and pressure
  // will never be less than "smallr" and "smallp", respectively
  FORT_POSTNORMALPREDF(CHF_FRA(a_dWMinus),
                       CHF_FRA(a_dWPlus),
                       CHF_CONST_FRA(a_W),
                       CHF_BOX(a_box));
}


// Compute a Riemann problem and generate fluxes at the faces
void SelfGravityPhysics::riemann(FArrayBox&       a_WGdnv,
                                 const FArrayBox& a_WLeft,
                                 const FArrayBox& a_WRight,
                                 const FArrayBox& a_W,
                                 const Real&      a_time,
                                 const int&       a_dir,
                                 const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_WGdnv.box().contains(a_box));

  // Get the numbers of relevant variables
  const int numPrim = numPrimitives();

  CH_assert(a_WGdnv .nComp() == numPrim);
  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft.box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));

  FORT_RIEMANNF(CHF_FRA(a_WGdnv),
                CHF_CONST_FRA(shiftWLeft),
                CHF_CONST_FRA(shiftWRight),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}


void SelfGravityPhysics::getFlux(FArrayBox&       a_flux,
                                 const FArrayBox& a_whalf,
                                 const int&       a_dir,
                                 const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_GETFLUXF(CHF_FRA(a_flux),
                CHF_CONST_FRA(a_whalf),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));
}


// Compute expansion amplitudes of dW in right eigenvectors.
void SelfGravityPhysics::charAnalysis(FArrayBox&       a_dW,
                                      const FArrayBox& a_W,
                                      const int&       a_dir,
                                      const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARANALYSISF(CHF_FRA(a_dW),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(a_box));
}

void SelfGravityPhysics::charSynthesis(FArrayBox&       a_dW,
                                       const FArrayBox& a_W,
                                       const int&       a_dir,
                                       const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARSYNTHESISF(CHF_FRA(a_dW),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));
}

void SelfGravityPhysics::charValues(FArrayBox&       a_lambda,
                                    const FArrayBox& a_W,
                                    const int&       a_dir,
                                    const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARVALUESF(CHF_FRA(a_lambda),
                   CHF_CONST_FRA(a_W),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(a_box));
}


// Compute increment of primitive variables.
void SelfGravityPhysics::quasilinearUpdate(FArrayBox&       a_dWdx,
                                           const FArrayBox& a_WHalf,
                                           const FArrayBox& a_W,
                                           const Real&      a_scale,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_dWdx.box().contains(a_box));

  FORT_GETADWDXF(CHF_FRA(a_dWdx),
                 CHF_CONST_FRA(a_WHalf),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_REAL(a_scale),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(a_box));
}

// synchronize entropy and pressure
void SelfGravityPhysics::synchronize(FArrayBox&        a_U,
                                     const FArrayBox&  a_UOld,
                                     const Box&        a_box)
{
  CH_assert(a_U.box().contains(a_box));

  FORT_SYNCHRONIZEF(CHF_FRA(a_U),CHF_BOX(a_box));
}

// Add to (increment) the source terms given the current state
void SelfGravityPhysics::incrementSource(FArrayBox& a_localS,
                                         const FArrayBox& a_S,
                                         const FArrayBox& a_W,
                                         const Real&      a_time,
                                         const Real&      a_dt,
                                         const Box&       a_box)
{
  // local copy
  a_localS.copy(a_S);

  // cast off constantness so that S can be upgraded
  FArrayBox& S = (FArrayBox&)a_S;
  S.setVal(zero);

  // gradient box
  const Box gbox = a_box & grow(a_W.box(),-1);

  FORT_COMPUTESOURCEF(CHF_FRA(S),
                      CHF_CONST_FRA(a_localS),
                      CHF_CONST_FRA(a_W),
                      CHF_BOX(gbox),
                      CHF_BOX(a_box));

  FORT_COMPUTELOCALSOURCEF(CHF_FRA(a_localS),
                           CHF_CONST_FRA(a_S),
                           CHF_CONST_FRA(a_W),
                           CHF_BOX(a_box));

  //m_source = &S;
}

// Set up the source term. Used in AMRLevel.advance()... to pass in
// information for the cooling term which is ultimately setup within
// PatchGodunov.updateState(), once the boundary conditions for the
// state have been accounted for.
void SelfGravityPhysics::setupSourceTerm(FArrayBox&       a_S,
                                         const FArrayBox&  a_force,
                                         const FArrayBox& a_U,
                                         const Real&      a_time,
                                         const Real&      a_dt,
                                         const Box&       a_box)
{
  //  a_S.setVal(zero);

  const int compVel0 = velocityInterval().begin();
  a_S.copy(a_force,0,compVel0,SpaceDim);
}

// apply gravity + exp source terms; if no heating, reset the
// temperature to at least a floor value
Real SelfGravityPhysics::applySource(FArrayBox&       a_U,
                                     const FArrayBox& a_UOld,
                                     const FArrayBox& a_source,
                                     const Real&      a_time,
                                     const Real&      a_dt,
                                     const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));

  FORT_APPLYSOURCEF(CHF_FRA(a_U),
                    CHF_CONST_FRA(a_UOld),
                    CHF_CONST_FRA(a_source),
                    CHF_CONST_REAL(a_dt),
                    CHF_BOX(a_box));

  FArrayBox W(a_box,numPrimitives());
  consToPrim(W,a_U,a_box);

  Real maxspeed;
  FORT_MAXWAVESPEEDWSOURCEF(CHF_REAL(maxspeed),
                            CHF_CONST_FRA(W),
                            CHF_CONST_FRA(a_source),
                            CHF_CONST_REAL(m_dx),
                            CHF_BOX(a_box));

  return (m_dx/maxspeed);
}

//
void SelfGravityPhysics::forceCorrection(FArrayBox&       a_U,
                                         const FArrayBox&  a_force,
                                         const Real&      a_time,
                                         const Real&      a_dt,
                                         const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));

  const FArrayBox& force = a_force;

  FORT_APPLYFORCECORRF(CHF_FRA(a_U),
                       CHF_CONST_FRA(force),
                       CHF_CONST_REAL(a_dt),
                       CHF_BOX(a_box));
}

//
void SelfGravityPhysics::setPressureToEntropy(FArrayBox&       a_U,
                                              const Box&       a_box)
{
  CH_assert(a_U.box().contains(a_box));

  FORT_RESETPRESSURE(CHF_FRA(a_U),CHF_BOX(a_box));
}

// Interval within the primitive variables corresponding to the velocities
Interval SelfGravityPhysics::velocityInterval()
{
#if CH_SPACEDIM==2
  Interval retval(WVELX,WVELY);
#elif CH_SPACEDIM==3
  Interval retval(WVELX,WVELZ);
#else
  bogus_spacedim();
#endif

  return retval;
}

// Used to limit the absolute value of a "pressure" difference
Real SelfGravityPhysics::smallPressure()
{
  CH_assert(isDefined());
  return m_smallPressure;
}

// Component index within the primitive variables of the pressure
int SelfGravityPhysics::densityIndex()
{
  CH_assert(isDefined());
  return WRHO;
}

// Component index within the primitive variables of the pressure
int SelfGravityPhysics::pressureIndex()
{
  CH_assert(isDefined());
  return WPRES;
}

// Component index within the primitive variables of the bulk modulus
int SelfGravityPhysics::bulkModulusIndex()
{
  CH_assert(isDefined());
  return WPRES;
}

// Component index within the primitive variables of the entropy
int SelfGravityPhysics::entropyIndex()
{
  CH_assert(isDefined());
  return WENTR;
}

// Used to limit the absolute value of a "pressure" difference
void SelfGravityPhysics::setSmallPressure(const Real& a_smallPressure)
{
  m_smallPressure = a_smallPressure;
}

