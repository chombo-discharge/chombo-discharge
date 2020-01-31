#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoHiSide.H"

#include "MOLAdvectionPhysics.H"
#include "AdvectOpF_F.H"

#include "NamespaceHeader.H"

MOLAdvectionPhysics::MOLAdvectionPhysics()
  :
  MOLPhysics()
{
}

MOLAdvectionPhysics::~MOLAdvectionPhysics()
{
}

// Compute the maximum wave speed
Real MOLAdvectionPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
                                          const Box&       a_box)
{
  MayDay::Error("MOLAdvectionPhysics::getMaxWaveSpeed does not exist");
  Real speed = 0.0;
  return speed;
}

// Return zero
void MOLAdvectionPhysics::soundSpeed(FArrayBox& a_speed,
                                     const FArrayBox& a_U,
                                     const Box&       a_box)
{
  // need this in calculation of new dt
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));
  CH_assert(a_speed.contains(a_box));
  a_speed.setVal(0., a_box, 0); // component 0
}

MOLPhysics* MOLAdvectionPhysics::new_molPhysics() const
{
  CH_assert(m_isBCSet);

  MOLPhysics* retval = static_cast<MOLPhysics*>
    (new MOLAdvectionPhysics());
  retval->setPhysIBC(m_bc);
  // added by petermc, 8 Jul 2009
  if (fourthOrderArtificialViscosityIsDefined())
    retval->setFourthOrderArtificialViscosityParameter(getFourthOrderArtificialViscosityParameter());

  return retval;
}

// Number of conserved variables
int MOLAdvectionPhysics::numConserved() const
{
  // only phi is conserved
  int numCons = 1;
  return numCons;
}

// Names of the conserved variables
Vector<string> MOLAdvectionPhysics::stateNames()
{
  CH_assert(isDefined());

  Vector<string> retval;

  retval.push_back("scalar");

  return retval;
}

// Number of flux variables
int MOLAdvectionPhysics::numFluxes() const
{
  CH_assert(isDefined());
  // In some computations there may be more fluxes than conserved variables
  int nflux = 1;
  return nflux;
}

void MOLAdvectionPhysics::getFlux(FArrayBox&       a_flux,
                                  const FArrayBox& a_whalf,
                                  const int&       a_dir,
                                  const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(m_isVelFabSet);
  // Set a_flux[0] := a_whalf[0] * advVel[a_dir], because we alias to 0.
  a_flux.copy(a_whalf); // a_flux[0] := a_whalf[0]
  a_flux.mult(*m_advVelFabPtr, a_dir, 0); // a_flux[0] *= advVel[a_dir]
}

// Number of primitive variables
int MOLAdvectionPhysics::numPrimitives() const
{
  int numPrim = 1;
  return numPrim;
}

// Compute a Riemann problem and generate fluxes at the faces.
void MOLAdvectionPhysics::riemann(FArrayBox&       a_WGdnv,
                                  const FArrayBox& a_WLeft,
                                  const FArrayBox& a_WRight,
                                  const FArrayBox& a_W,
                                  const Real&      a_time,
                                  const int&       a_dir,
                                  const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_WGdnv.box().contains(a_box));
  CH_assert(m_isVelFabSet);
  const FArrayBox& advVelFab = *m_advVelFabPtr;
  CH_assert(advVelFab.box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();
  CH_assert(a_WGdnv .nComp() == numPrim);
  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&) a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&) a_WRight;

  // Solution to the Riemann problem

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft .box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));

  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.

  FORT_SELECTUPWIND(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WLeft),
                    CHF_CONST_FRA(a_WRight),
                    CHF_CONST_FRA1(advVelFab, a_dir),
                    CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

// Compute the primitive variables from the conserved variables
void MOLAdvectionPhysics::consToPrim(FArrayBox&       a_W,
                                     const FArrayBox& a_U,
                                     const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  a_W.copy(a_U, a_box);
}

// Interval within the primitive variables corresponding to the velocities
Interval MOLAdvectionPhysics::velocityInterval() const
{
  // return bogus interval
  Interval retval(0, -1);

  return retval;
}

// Component index within the primitive variables of the pressure
int MOLAdvectionPhysics::pressureIndex()
{
  CH_assert(isDefined());

  return -1;
}

// Used to limit the absolute value of a "pressure" difference
Real MOLAdvectionPhysics::smallPressure()
{
  CH_assert(isDefined());

  return 0.;
}

// Component index within the primitive variables of the bulk modulus
int MOLAdvectionPhysics::bulkModulusIndex()
{
  CH_assert(isDefined());

  return -1;
}

void MOLAdvectionPhysics::setVelocityFab(FArrayBox* a_advVelFabPtr)
{
  m_advVelFabPtr = a_advVelFabPtr;
  m_isVelFabSet = true;
  m_isVelFluxSet = false;
}

void MOLAdvectionPhysics::setVelocityFlux(FluxBox* a_advVelFluxPtr)
{
  m_advVelFluxPtr = a_advVelFluxPtr;
  m_isVelFluxSet = true;
  m_isVelFabSet = false;
}

#ifdef CH_USE_HDF5
void MOLAdvectionPhysics::expressions(HDF5HeaderData& a_expressions) const
{
  a_expressions.m_string["scalar gamma"] = "1.4";

  a_expressions.m_string["vector velocity"] = "momentum/density";
}

#endif

#include "NamespaceFooter.H"
