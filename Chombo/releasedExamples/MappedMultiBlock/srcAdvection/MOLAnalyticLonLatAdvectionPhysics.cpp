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

#include "MOLAnalyticLonLatAdvectionPhysics.H"
#include "AdvectOpF_F.H"
#include "BoxIterator.H"
#include "SetCentersF_F.H"
#include "CubedSphere2DCS.H"
#include "CubedSphere2DPanelCS.H"
#include "MBMiscUtil.H"
#include "MayDay.H"
#include "DebugOut.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
MOLAnalyticLonLatAdvectionPhysics::
MOLAnalyticLonLatAdvectionPhysics(RefCountedPtr<VectorFunction> a_advectionVelocity):
  MOLMappedPhysics(),
  m_velocityFunc(a_advectionVelocity),
  m_streamFunc(),
  m_coordSys(NULL),
  m_time(0.0)
{
  CH_assert(!a_advectionVelocity.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MOLAnalyticLonLatAdvectionPhysics::
MOLAnalyticLonLatAdvectionPhysics(RefCountedPtr<ScalarFunction> a_streamFunc):
  MOLMappedPhysics(),
  m_velocityFunc(),
  m_streamFunc(a_streamFunc),
  m_coordSys(NULL),
  m_time(0.0)
{
  CH_assert(!a_streamFunc.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MOLAnalyticLonLatAdvectionPhysics::
~MOLAnalyticLonLatAdvectionPhysics()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticLonLatAdvectionPhysics::
setCurrentTime(Real a_time)
{
  m_time = a_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticLonLatAdvectionPhysics::
setCurrentCoordSys(const NewCoordSys* a_coordSys)
{
  m_coordSys = const_cast<NewFourthOrderCoordSys*>(dynamic_cast<const NewFourthOrderCoordSys*>(a_coordSys));
  CH_assert(m_coordSys != NULL);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Compute the maximum wave speed
Real
MOLAnalyticLonLatAdvectionPhysics::
getMaxWaveSpeed(const FArrayBox& a_U,
                const Box&       a_box)
{
  MayDay::Error("MOLAnalyticLonLatAdvectionPhysics::getMaxWaveSpeed does not exist");
  Real speed = 0.0;
  return speed;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticLonLatAdvectionPhysics::
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
MOLAnalyticLonLatAdvectionPhysics::new_molPhysics() const
{
  CH_assert(m_isBCSet);

  MOLPhysics* retval;
  if (!m_velocityFunc.isNull())
    retval = dynamic_cast<MOLPhysics*>(new MOLAnalyticLonLatAdvectionPhysics(m_velocityFunc));
  else
    retval = dynamic_cast<MOLPhysics*>(new MOLAnalyticLonLatAdvectionPhysics(m_streamFunc));
  retval->setPhysIBC(m_bc);
  // added by petermc, 8 Jul 2009
  if (fourthOrderArtificialViscosityIsDefined())
    retval->setFourthOrderArtificialViscosityParameter(getFourthOrderArtificialViscosityParameter());
  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
MOLAnalyticLonLatAdvectionPhysics::
numConserved() const
{
  // only phi is conserved
  return 1;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<string>
MOLAnalyticLonLatAdvectionPhysics::
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
MOLAnalyticLonLatAdvectionPhysics::
numFluxes() const
{
  CH_assert(isDefined());
  return 1; // Normal flux only!
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MOLAnalyticLonLatAdvectionPhysics::
getVelocity(FArrayBox&       a_velFab,
            int              a_dir,
            const Box&       a_box)
{
  CH_assert(isDefined());
  int faceDir = faceDimension(a_box); // in MBMiscUtil
  // xiFab:  cartesian coordinates (face centers)
  FArrayBox xiFab(a_box, SpaceDim);
  FORT_SETFACECENTERS(CHF_FRA(xiFab),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_INT(faceDir),
                      CHF_BOX(a_box));

  // Extract the coordinate system for this box.
  const CubedSphere2DPanelCS* panelCS = dynamic_cast<const CubedSphere2DPanelCS*>(m_coordSys);
  CH_assert(panelCS != NULL);

  // Lon/lat coordinates for face centers.
  FArrayBox lonlatFab(a_box, SpaceDim);
  panelCS->fabTransformEquiangularToLonLat(xiFab, lonlatFab);

  // Set lonlatVelocityFab to velocity in lon/lat coordinates.
  FArrayBox lonlatVelocityFab(a_box, SpaceDim);
  lonlatVelocityFab.setVal(0.0);
  if (m_velocityFunc->isHomogeneous())
    {
      RealVect U0 = (*m_velocityFunc)(RealVect(), m_time);
      // Convert to lon/lat coordinates.
      for (int idirVel = 0; idirVel < SpaceDim; ++idirVel)
        lonlatVelocityFab.setVal(U0[idirVel], idirVel);
    }
  else
    {
      for (BoxIterator bit(a_box); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect lonlat(D_DECL(lonlatFab(iv, 0), lonlatFab(iv, 1), 0));
          RealVect Ui = (*m_velocityFunc)(lonlat, m_time);
          for (int idirVel = 0; idirVel < SpaceDim; ++idirVel)
            lonlatVelocityFab(iv, idirVel) = Ui[idirVel];
        }
    }

  // Transform velocity representation from lon/lat
  // to equiangular coordinates.
  // FIXME: This transformation is poorly named, as the velocity vector is assumed
  // FIXME: to be in lon/lat coordinates (not lat/lon)!
  // FArrayBox eqVelocityFab(a_flux.box(), SpaceDim);
  panelCS->fabVectorTransformLatLonToEquiangular(xiFab,
                                                 lonlatVelocityFab,
                                                 a_velFab);
}


//-----------------------------------------------------------------------
void
MOLAnalyticLonLatAdvectionPhysics::
getFlux(FArrayBox&       a_flux,
        const FArrayBox& a_whalf,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_assert(isDefined());
  // Evaluate the velocity function on faces and use it to set the flux vector.
  // Note that we need the vector-valued flux, not just the normal flux. This
  // is because we transform the flux vector from lon/lat coordinates to
  // equiangular coordinates.
  if (!m_velocityFunc.isNull())
    {
      FArrayBox velocityFab(a_flux.box(), SpaceDim);
      getVelocity(velocityFab, a_dir, a_flux.box());

      // Record the normal flux:
      // a_flux[0] := a_whalf[0] * eqVelocityFab[a_dir]
      a_flux.setVal(0., 0);
      a_flux.copy(a_whalf, 0, 0);
      a_flux.mult(velocityFab, a_dir, 0);
    }
  else
    {
      MayDay::Error("MOLAnalyticLonLatAdvectionPhysics::getFlux: Implement stream function physics!");
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Number of primitive variables
int
MOLAnalyticLonLatAdvectionPhysics::
numPrimitives() const
{
  int numPrim = 1;
  return numPrim;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Compute a Riemann problem and generate fluxes at the faces.
void
MOLAnalyticLonLatAdvectionPhysics::
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

  // Evaluate the velocity function on faces.
  FArrayBox advVelFab(a_box, SpaceDim);
  getVelocity(advVelFab, a_dir, a_box);

  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.
  FORT_SELECTUPWIND(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WLeft),
                    CHF_CONST_FRA(a_WRight),
                    CHF_CONST_FRA1(advVelFab, a_dir),
                    CHF_BOX(a_box));

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Compute the primitive variables from the conserved variables
void MOLAnalyticLonLatAdvectionPhysics::consToPrim(FArrayBox&       a_W,
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
Interval MOLAnalyticLonLatAdvectionPhysics::velocityInterval() const
{
  // return bogus interval
  Interval retval(0, -1);

  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Component index within the primitive variables of the pressure
int MOLAnalyticLonLatAdvectionPhysics::pressureIndex()
{
  CH_assert(isDefined());

  return -1;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Used to limit the absolute value of a "pressure" difference
Real MOLAnalyticLonLatAdvectionPhysics::smallPressure()
{
  CH_assert(isDefined());

  return 0.;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Component index within the primitive variables of the bulk modulus
int MOLAnalyticLonLatAdvectionPhysics::bulkModulusIndex()
{
  CH_assert(isDefined());

  return -1;
}
//-----------------------------------------------------------------------

#ifdef CH_USE_HDF5
//-----------------------------------------------------------------------
void MOLAnalyticLonLatAdvectionPhysics::expressions(HDF5HeaderData& a_expressions) const
{
}
//-----------------------------------------------------------------------

#endif

#include "NamespaceFooter.H"
