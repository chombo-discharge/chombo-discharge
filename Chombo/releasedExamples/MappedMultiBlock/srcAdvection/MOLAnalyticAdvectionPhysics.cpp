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

#include "MOLAnalyticAdvectionPhysics.H"
#include "AdvectOpF_F.H"
#include "BoxIterator.H"
#include "SetCentersF_F.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
MOLAnalyticAdvectionPhysics::
MOLAnalyticAdvectionPhysics(RefCountedPtr<VectorFunction> a_advectionVelocity):
  MOLMappedPhysics(),
  m_velocityFunc(a_advectionVelocity),
  m_coordSys(NULL),
  m_time(0.0)
{
  CH_assert(!a_advectionVelocity.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MOLAnalyticAdvectionPhysics::
~MOLAnalyticAdvectionPhysics()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticAdvectionPhysics::
setCurrentTime(Real a_time)
{
  m_time = a_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticAdvectionPhysics::
setCurrentCoordSys(const NewCoordSys* a_coordSys)
{
  m_coordSys = const_cast<NewFourthOrderCoordSys*>(dynamic_cast<const NewFourthOrderCoordSys*>(a_coordSys));
  CH_assert(m_coordSys != NULL);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Compute the maximum wave speed
Real
MOLAnalyticAdvectionPhysics::
getMaxWaveSpeed(const FArrayBox& a_U,
                const Box&       a_box)
{
  MayDay::Error("MOLAnalyticAdvectionPhysics::getMaxWaveSpeed does not exist");
  Real speed = 0.0;
  return speed;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticAdvectionPhysics::
soundSpeed(FArrayBox& a_speed,
           const FArrayBox& a_U,
           const Box&       a_box)
{
  // need this in calculation of new dt
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));
  CH_assert(a_speed.contains(a_box));
  a_speed.setVal(0., a_box, 0); // component 0
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MOLPhysics*
MOLAnalyticAdvectionPhysics::
new_molPhysics() const
{
  CH_assert(m_isBCSet);

  MOLPhysics* retval = dynamic_cast<MOLPhysics*>(new MOLAnalyticAdvectionPhysics(m_velocityFunc));
  retval->setPhysIBC(m_bc);
  // added by petermc, 8 Jul 2009
  if (fourthOrderArtificialViscosityIsDefined())
    retval->setFourthOrderArtificialViscosityParameter(getFourthOrderArtificialViscosityParameter());
  MOLAnalyticAdvectionPhysics* phys = static_cast<MOLAnalyticAdvectionPhysics*>(retval);
  phys->m_coordSys = m_coordSys;
  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
MOLAnalyticAdvectionPhysics::
numConserved() const
{
  // only phi is conserved
  return 1;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<string>
MOLAnalyticAdvectionPhysics::
stateNames()
{
  CH_assert(isDefined());
  Vector<string> retval;
  retval.push_back("scalar");
  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
MOLAnalyticAdvectionPhysics::
numFluxes() const
{
  CH_assert(isDefined());
  return 1; // Normal flux only!
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticAdvectionPhysics::
getFlux(FArrayBox&       a_flux,
        const FArrayBox& a_whalf,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_assert(isDefined());

  // Xi:  cartesian coordinates (face centers)
  FArrayBox Xi(a_box, SpaceDim);
  FORT_SETFACECENTERS(CHF_FRA(Xi),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));

  // X:  physical coordinates (face centers)
  FArrayBox X(a_box, SpaceDim);
  if (m_coordSys != NULL)
    m_coordSys->realCoord(X, Xi, a_box);
  else
    X.copy(Xi); // Identity mapping.

  // Evaluate the flux function F = V*U on faces.
  if (m_velocityFunc->isHomogeneous()) // Velocity is constant in space
  {
    RealVect U0 = (*m_velocityFunc)(RealVect(), m_time);
    if (m_coordSys != NULL)
    {
      for (BoxIterator bit(a_box); bit.ok(); ++bit)
      {
        IntVect i = bit();
        RealVect Xii;
        for (int d = 0; d < SpaceDim; ++d)
          Xii[d] = Xi(i, d);
        RealVect V = m_coordSys->vectorTransformationRealToMapped(Xii).transform(U0); // Mapped velocity at Xii.
        a_flux(i, 0) = a_whalf(i, 0) * V[a_dir];
      }
    }
    else
    {
      for (int d = 0; d < SpaceDim; ++d)
        a_flux.setVal(U0[a_dir], 0);
      a_flux.mult(a_whalf);
    }
  }
  else // Velocity varies in space.
  {
    for (BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      IntVect i = bit();
      RealVect x;
      for (int d = 0; d < SpaceDim; ++d)
        x[d] = X(i, d);
      RealVect Ui = (*m_velocityFunc)(x, m_time);
      if (m_coordSys != NULL)
      {
        RealVect Xii;
        for (int d = 0; d < SpaceDim; ++d)
          Xii[d] = Xi(i, d);
        RealVect Vi = m_coordSys->vectorTransformationRealToMapped(Xii).transform(Ui); // Mapped velocity at Xii.
        a_flux(i, 0) = a_whalf(i, 0) * Vi[a_dir];
      }
      else
        a_flux(i, 0) = a_whalf(i, 0) * Ui[a_dir];
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Number of primitive variables
int
MOLAnalyticAdvectionPhysics::
numPrimitives() const
{
  int numPrim = 1;
  return numPrim;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Compute a Riemann problem and generate fluxes at the faces.
void
MOLAnalyticAdvectionPhysics::
riemann(FArrayBox&       a_WGdnv,
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

  // Compute the advection velocity on cell centers.
  FArrayBox advVel(a_box, SpaceDim);

  // Xi:  cartesian coordinates (face centers)
  FArrayBox Xi(a_box, SpaceDim);
  FORT_SETFACECENTERS(CHF_FRA(Xi),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));

  // X:  physical coordinates (face centers)
  FArrayBox X(a_box, SpaceDim);
  if (m_coordSys != NULL)
    m_coordSys->realCoord(X, Xi, a_box);
  else
    X.copy(Xi); // Identity mapping.

  // Evaluate the velocity function on faces.
  if (m_velocityFunc->isHomogeneous())
  {
    RealVect U0 = (*m_velocityFunc)(RealVect(), m_time);
    if (m_coordSys != NULL)
    {
      for (BoxIterator bit(a_box); bit.ok(); ++bit)
      {
        IntVect i = bit();
        RealVect Xii;
        for (int d = 0; d < SpaceDim; ++d)
          Xii[d] = Xi(i, d);
        RealVect V = m_coordSys->vectorTransformationRealToMapped(Xii).transform(U0); // Mapped velocity at Xii.
        for (int d = 0; d < SpaceDim; ++d)
          advVel.setVal(V[d], d);
      }
    }
    else
    {
      for (int d = 0; d < SpaceDim; ++d)
        advVel.setVal(U0[d], d);
    }
  }
  else
  {
    for (BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      IntVect i = bit();
      RealVect x;
      for (int d = 0; d < SpaceDim; ++d)
        x[d] = X(i, d);
      RealVect Ui = (*m_velocityFunc)(x, m_time);
      if (m_coordSys != NULL)
      {
        RealVect Xii;
        for (int d = 0; d < SpaceDim; ++d)
          Xii[d] = Xi(i, d);
        RealVect Vi = m_coordSys->vectorTransformationRealToMapped(Xii).transform(Ui); // Mapped velocity at Xii.
        for (int d = 0; d < SpaceDim; ++d)
          advVel(i, d) = Vi[d];
      }
      else
      {
        for (int d = 0; d < SpaceDim; ++d)
          advVel(i, d) = Ui[d];
      }
    }
  }

  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.
  FORT_SELECTUPWIND(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WLeft),
                    CHF_CONST_FRA(a_WRight),
                    CHF_CONST_FRA1(advVel, a_dir),
                    CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
//  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
//  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Compute the primitive variables from the conserved variables
void MOLAnalyticAdvectionPhysics::consToPrim(FArrayBox&       a_W,
                                     const FArrayBox& a_U,
                                     const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  a_W.copy(a_U, a_box);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Interval within the primitive variables corresponding to the velocities
Interval MOLAnalyticAdvectionPhysics::velocityInterval() const
{
  // return bogus interval
  Interval retval(0, -1);

  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Component index within the primitive variables of the pressure
int MOLAnalyticAdvectionPhysics::pressureIndex()
{
  CH_assert(isDefined());

  return -1;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Used to limit the absolute value of a "pressure" difference
Real MOLAnalyticAdvectionPhysics::smallPressure()
{
  CH_assert(isDefined());

  return 0.;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Component index within the primitive variables of the bulk modulus
int MOLAnalyticAdvectionPhysics::bulkModulusIndex()
{
  CH_assert(isDefined());

  return -1;
}
//-----------------------------------------------------------------------

#ifdef CH_USE_HDF5
//-----------------------------------------------------------------------
void MOLAnalyticAdvectionPhysics::expressions(HDF5HeaderData& a_expressions) const
{
}
//-----------------------------------------------------------------------

#endif

#include "NamespaceFooter.H"
