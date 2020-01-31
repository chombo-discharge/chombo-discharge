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

#include "MOLShallowWaterPhysics.H"

#include "ShallowWaterPhysicsF_F.H"
#include "CubedSphere2DPanelCS.H"
#include "SWintegrator.H"

#include "NamespaceHeader.H"

MOLShallowWaterPhysics::
MOLShallowWaterPhysics(RefCountedPtr<ScalarFunction> a_height,
                       RefCountedPtr<VectorFunction> a_velocity) :
  MOLPhysics(),
  m_heightFunc(a_height),
  m_velocityFunc(a_velocity),
  m_coordSys(NULL)
{
  CH_assert(!a_height.isNull());
  CH_assert(!a_velocity.isNull());
}

MOLShallowWaterPhysics::~MOLShallowWaterPhysics()
{
}

void
MOLShallowWaterPhysics::
setCurrentCoordSys(const NewCoordSys* a_coordSys)
{
  m_coordSys = const_cast<NewFourthOrderCoordSys*>(dynamic_cast<const NewFourthOrderCoordSys*>(a_coordSys));
  CH_assert(m_coordSys != NULL);
}

// Compute the maximum wave speed
Real MOLShallowWaterPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
                                             const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));

  Real speed = 0.0;

  FORT_SWMAXWAVESPEEDF(CHF_REAL(speed),
                       CHF_CONST_FRA(a_U),
                       CHF_BOX(a_box));

  return speed;
}

// Compute the speed of sound
void MOLShallowWaterPhysics::soundSpeed(FArrayBox& a_speed,
                                      const FArrayBox& a_U,
                                      const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));
  CH_assert(a_speed.contains(a_box));

  FORT_SWSOUNDSPEEDF(CHF_CONST_FRA1(a_speed, 0),
                     CHF_CONST_FRA(a_U),
                     CHF_BOX(a_box));
}

MOLPhysics* MOLShallowWaterPhysics::new_molPhysics() const
{
  CH_assert(m_isBCSet);

  MOLPhysics* retval = dynamic_cast<MOLPhysics*>(new MOLShallowWaterPhysics(m_heightFunc, m_velocityFunc));
  retval->setPhysIBC(m_bc);
  // added by petermc, 8 Jul 2009
  if (fourthOrderArtificialViscosityIsDefined())
    retval->setFourthOrderArtificialViscosityParameter(getFourthOrderArtificialViscosityParameter());
  MOLShallowWaterPhysics* phys =
    static_cast<MOLShallowWaterPhysics*>(retval);
  phys->m_coordSys = m_coordSys;
  return retval;
}

// Number of conserved variables
int MOLShallowWaterPhysics::numConserved() const
{
  CH_assert(isDefined());

  return UNUM;
}

// Names of the conserved variables
Vector<string> MOLShallowWaterPhysics::stateNames()
{
  CH_assert(isDefined());

  Vector<string> retval;

  retval.push_back("height");

  retval.push_back("x-momentum");
  if (SpaceDim >= 2)
  {
    retval.push_back("y-momentum");

    if (SpaceDim >= 3)
    {
      retval.push_back("z-momentum");
    }
  }

  return retval;
}

// Number of flux variables
int MOLShallowWaterPhysics::numFluxes() const
{
  CH_assert(isDefined());

  // In some computations there may be more fluxes than conserved variables
  return UNUM;
}

int MOLShallowWaterPhysics::densityIndex()
{
  return WHGT;
}

void MOLShallowWaterPhysics::getFlux(FArrayBox&       a_flux,
                                     const FArrayBox& a_whalf,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  CH_assert(isDefined());

  const Box& bx = a_flux.box();
  CH_assert(m_coordSys != NULL);
  FArrayBox metric(bx, SpaceDim);
  ((CubedSphere2DPanelCS*) m_coordSys)->contravariantMetric(metric, a_dir);

  // petermc, 1 Sep 2011: Note that MOLAnalyticAdvectionPhysics looks
  // a lot more complicated here.
  FORT_SWGETFLUXF(CHF_FRA(a_flux),
                  CHF_CONST_FRA(a_whalf),
                  CHF_CONST_FRA(metric), // contravariant 2D metric
                  CHF_CONST_INT(a_dir),
                  CHF_BOX(a_box));

  // Multiply a_flux by J.
  FArrayBox JavgFab(bx, 1);
  m_coordSys->getAvgJ(JavgFab, bx);
  for (int comp = 0; comp < SpaceDim; comp++)
    {
      // a_flux[comp] *= JavgFab[0]
      a_flux.mult(JavgFab, 0, comp);
    }
}

// Number of primitive variables
int MOLShallowWaterPhysics::numPrimitives() const
{
  CH_assert(isDefined());

  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return WNUM;
}

// Compute a Riemann problem and generate fluxes at the faces.
void MOLShallowWaterPhysics::riemann(FArrayBox&       a_WGdnv,
                                   const FArrayBox& a_WLeft,
                                   const FArrayBox& a_WRight,
                                   const FArrayBox& a_W,
                                   //const FArrayBox& a_faceN,
                                   //const IntVect&   a_metricTermComponents,
                                   const Real&      a_time,
                                   const int&       a_dir,
                                   const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_WGdnv.box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();

  CH_assert(a_WGdnv .nComp() == numPrim);
  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);
/*
  FArrayBox unitNormalFab(a_box, SpaceDim*SpaceDim);
  getBasisVectors(unitNormalFab, a_faceN, a_metricTermComponents, a_dir);
*/
  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Solution to the Riemann problem

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft .box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));
/*
  // Transform velocity components of shiftWLeft and shiftWRight.
  Interval velInt = velocityInterval();
  FArrayBox velLeftFab(velInt, shiftWLeft);
  FArrayBox velRightFab(velInt, shiftWRight);
  forwardBasisTransform(velLeftFab, unitNormalFab);
  forwardBasisTransform(velRightFab, unitNormalFab);
*/
  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.
  // petermc, 1 Sep 2011: Note that MOLAnalyticAdvectionPhysics looks
  // a lot more complicated here.
  FORT_SWRIEMANNROEF(CHF_FRA(a_WGdnv),
                   CHF_CONST_FRA(shiftWLeft),
                   CHF_CONST_FRA(shiftWRight),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
  //  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
  //  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);
/*
  // Transform back velocity components in a_WGdnv.
  FArrayBox velNewFab(velInt, a_WGdnv);
  reverseBasisTransform(velNewFab, unitNormalFab);
*/
  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

// Compute the primitive variables from the conserved variables
void MOLShallowWaterPhysics::consToPrim(FArrayBox&       a_W,
                                        const FArrayBox& a_U,
                                        const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  FORT_SWCONSTOPRIMF(CHF_FRA(a_W),
                   CHF_CONST_FRA(a_U),
                   CHF_BOX(a_box));
}

// Interval within the primitive variables corresponding to the velocities
Interval MOLShallowWaterPhysics::velocityInterval() const
{
  CH_assert(isDefined());

#if CH_SPACEDIM==1
  Interval retval(WVELX,WVELX);
#elif CH_SPACEDIM==2
  Interval retval(WVELX,WVELY);
#elif CH_SPACEDIM==3
  Interval retval(WVELX,WVELZ);
#else
  bogus_spacedim();
#endif

  return retval;
}

// Component index within the primitive variables of the pressure
int MOLShallowWaterPhysics::pressureIndex()
{
  CH_assert(isDefined());

  return -1;
}

// Used to limit the absolute value of a "pressure" difference
Real MOLShallowWaterPhysics::smallPressure()
{
  CH_assert(isDefined());

  return 0.;
}

// Component index within the primitive variables of the bulk modulus
int MOLShallowWaterPhysics::bulkModulusIndex()
{
  CH_assert(isDefined());

  return (-1);
}

#ifdef CH_USE_HDF5
void MOLShallowWaterPhysics::expressions(HDF5HeaderData& a_expressions) const
{
  a_expressions.m_string["scalar gravity"] = "11489.5722";

  a_expressions.m_string["vector velocity"] = "momentum/height";
  a_expressions.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  a_expressions.m_string["scalar pressure"] = "gravity*height*height/2";
  a_expressions.m_string["scalar soundspeed"] = "sqrt(gravity*height)";
}

#endif

#include "NamespaceFooter.H"
