#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MappedAdvectionPhysics.H"
#include "MappedAdvectionPhysicsF_F.H"
#include "FourthOrderCoordSys.H"

#include "LoHiSide.H"

MappedAdvectionPhysics::MappedAdvectionPhysics():MappedGodunovPhysics()
{
  // Bogus values
  m_numCons = -1;
  m_numPrim = -1;
  m_numFlux = -1;

  // Undefined until these are defined
  m_nCompDefined = false;
  //m_advVelPtr = NULL;
  m_NTvelPtr = NULL;
  m_cellNTvelPtr = NULL;
  m_JPtr = NULL;
  m_NPtr = NULL;
  m_NJInversePtr = NULL;

}

MappedAdvectionPhysics::~MappedAdvectionPhysics()
{
}

void MappedAdvectionPhysics::setNComp(int nComp)
{
  CH_assert(nComp > 0);

  // All these are the same variable for advection
  m_numCons = nComp;
  m_numPrim = nComp;
  m_numFlux = nComp;

  m_nCompDefined = true;
}

// Compute the maximum wave speed
Real MappedAdvectionPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
                                             const Box&       a_box)
{
  CH_assert(isDefined());

  // note that for this advection problem, the max wavespeed
  // is based entirely on the cell-centered velocity
  Real maxSpeed, speed = 0.0;
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox& vel = *m_cellNTvelPtr;
      speed = vel.norm(0,dir,1);
      if (maxSpeed < speed) maxSpeed = speed;
    }

  return maxSpeed;
}

// Factory method - this object is its own factory.  It returns a pointer
// to new MappedAdvectionPhysics object with the same definition as this object.
GodunovPhysics* MappedAdvectionPhysics::new_godunovPhysics() const
{
  // Make the new object
  MappedAdvectionPhysics* newPtr = new MappedAdvectionPhysics();

  // Define the same domain and dx
  newPtr->define(m_domain,m_dx);

  // Pass the current initial and boundary condition (IBC) object to this new
  // patch integrator so it can create its own IBC object and define it as it
  // needs to (i.e. new domain and grid spacing if necessary)
  newPtr->setPhysIBC(getPhysIBC());

  newPtr->m_numCons = m_numCons;
  newPtr->m_numPrim = m_numPrim;
  newPtr->m_numFlux = m_numFlux;

  newPtr->m_nCompDefined = m_nCompDefined;

  GodunovPhysics* retval = static_cast<GodunovPhysics*>(newPtr);

  // Return the new object
  return retval;
}

// Number of conserved variables
int MappedAdvectionPhysics::numConserved()
{
  return m_numCons;
}

// Names of the conserved variables
Vector<string> MappedAdvectionPhysics::stateNames()
{
  Vector<string> retval(m_numPrim, "scalar");

  return retval;
}

// Number of flux variables
int MappedAdvectionPhysics::numFluxes()
{
  // In some computations there may be more fluxes than conserved variables
  return m_numFlux;
}

void MappedAdvectionPhysics::getFlux(FArrayBox&       a_flux,
                                     const FArrayBox& a_WHalf,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_flux.nComp() == a_WHalf.nComp());

  //FArrayBox& faceMetrics = (*m_NPtr)[a_dir];
  //FArrayBox& advVelDir = (*m_advVelPtr)[a_dir];
  //FArrayBox tempFlux(a_box, a_flux.nComp());

  FArrayBox& NTvel = (*m_NTvelPtr)[a_dir];

  // zero out flux
  a_flux.setVal(0.0);

  // flux is N^T*advVel*WHalf
  a_flux.copy(a_WHalf, a_box);

  // multiply all flux components by NTvel
  for (int comp=0; comp<a_flux.nComp(); comp++)
    {
      a_flux.mult(NTvel, a_box, 0, comp, 1);
    } // end loop over flux components

}


// Is everthing defined
bool MappedAdvectionPhysics::isDefined() const
{
  return (m_isDefined && m_nCompDefined);
}

// Number of primitive variables
int MappedAdvectionPhysics::numPrimitives()
{
  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return m_numPrim;
}

void MappedAdvectionPhysics::charAnalysis(FArrayBox&       a_dW,
                                          const FArrayBox& a_W,
                                          const int&       a_dir,
                                          const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_dW.nComp() == a_W.nComp());

  // This is the identity mapping for advection - do nothing
}

void MappedAdvectionPhysics::charSynthesis(FArrayBox&       a_dW,
                                           const FArrayBox& a_W,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_dW.nComp() == a_W.nComp());

  // This is the identity mapping for advection - do nothing
}

void MappedAdvectionPhysics::charValues(FArrayBox&       a_lambda,
                                        const FArrayBox& a_W,
                                        const int&       a_dir,
                                        const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_lambda.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  CH_assert(m_cellNTvelPtr != NULL);
  FArrayBox& cellNTvel = *m_cellNTvelPtr;

  int nComp = a_lambda.nComp();
  int copyComp = 1;

  for (int comp = 0; comp < nComp; comp++)
  {
    a_lambda.copy(cellNTvel,a_box,a_dir,a_box,comp,copyComp);
  }
}

// Increment the source term using the primitive state (not needed)
void MappedAdvectionPhysics::incrementSource(FArrayBox&       a_S,
                                             const FArrayBox& a_W,
                                             const Box&       a_box)
{
  CH_assert (isDefined());

  CH_assert (a_S.box().contains(a_box));
  CH_assert (a_W.box().contains(a_box));

  CH_assert (a_S.nComp() == a_W.nComp());

  // The source term does not explicitly depend on the current primitive state
  // multiply by J...
  for (int comp=0; comp<a_S.nComp(); comp++)
    {
      a_S.mult((*m_JPtr), a_box, 0, comp, 1);
    }
}

// Compute a Riemann problem and generate fluxes at the faces
void MappedAdvectionPhysics::riemann(FArrayBox&       a_WGdnv,
                                     const FArrayBox& a_WLeft,
                                     const FArrayBox& a_WRight,
                                     const FArrayBox& a_W,
                                     const Real&      a_time,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_WGdnv.box().contains(a_box));

  CH_assert(m_NTvelPtr != NULL);
  FluxBox& NTvel = *m_NTvelPtr;
  // CH_assert(advVel.box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();

  CH_assert(a_WGdnv .nComp() == numPrim);
  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

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

  //FArrayBox& advVelDir = advVel[a_dir];
  //CH_assert(advVelDir.box().contains(a_box));
  //FArrayBox& NTvelface = NTvel[a_dir];

  FArrayBox& NTvelDir = NTvel[a_dir];
  CH_assert(NTvelDir.box().contains(a_box));

  // Riemann solver computes WGdnv all edges that are not on the physical
  // boundary.
  FORT_MAPPEDADVRIEMANNF(CHF_FRA(a_WGdnv),
                         CHF_CONST_FRA(shiftWLeft),
                         CHF_CONST_FRA(shiftWRight),
                         CHF_CONST_FRA1(NTvelDir,0),
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

void MappedAdvectionPhysics::postNormalPred(FArrayBox&       a_WMinus,
                                            FArrayBox&       a_WPlus,
                                            const FArrayBox& a_W,
                                            const Real&      a_dt,
                                            const Real&      a_dx,
                                            const int&       a_dir,
                                            const Box&       a_box)
{


#if 0
  // this breaks freestream preservation, and is probably not the right
  // thing to be doing.
  // (DFM -- 9/17/09) -- reverted to advective form for quasilinear update
  // so this is no longer required...

  CH_assert(isDefined());
  // add in mapping-dependent piece...
  FArrayBox src(a_box, 1);

  // (N^T)(vel) in the dir direction
  FArrayBox& normalVel = (*m_NTvelPtr)[a_dir];

  // reality check
  Box faceBox(a_box);
  faceBox.surroundingNodes(a_dir);
  CH_assert(normalVel.box().contains(faceBox));

  // grad(N^Tvel)
  FORT_MAPPEDADVSIMPLEGRAD(CHF_FRA(src), CHF_CONST_FRA(normalVel),
                           CHF_CONST_REAL(m_dx), CHF_CONST_INT(a_dir),
                           CHF_BOX(a_box));


  // now multiply by 0.5*dt*W/J and increment predicted values
  for (int comp=0; comp<a_WMinus.nComp(); comp++)
    {
      src.mult(a_W, a_box, comp, 0, 1);
      src.divide(*m_JPtr, a_box, 0,0,1);
      src *= -0.5*a_dt;

      a_WMinus.plus(src, 0, comp, 1);
      a_WPlus.plus(src, 0, comp, 1);
    }

#endif

}

void MappedAdvectionPhysics::quasilinearUpdate(FArrayBox&       a_AdWdx,
                                               const FArrayBox& a_WHalf,
                                               const FArrayBox& a_W,
                                               const Real&      a_scale,
                                               const int&       a_dir,
                                               const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_AdWdx.box().contains(a_box));
  CH_assert(a_W    .box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();

  CH_assert(a_AdWdx.nComp() == numPrim);
  CH_assert(a_WHalf.nComp() == numPrim);
  CH_assert(a_W    .nComp() == numPrim);

  const FArrayBox& cellNTvel = *m_cellNTvelPtr;

  const FArrayBox& faceNTvel = (*m_NTvelPtr)[a_dir];

  CH_assert(faceNTvel.nComp() == 1);

  // advective form uses cell velocity
  FORT_MAPPEDADVQUASILINEARUPDATE(CHF_FRA(a_AdWdx),
                                  CHF_CONST_FRA(a_WHalf),
                                  CHF_CONST_FRA1(cellNTvel,0),
                                  CHF_CONST_REAL(a_scale),
                                  CHF_CONST_REAL(m_dx),
                                  CHF_CONST_INT(a_dir),
                                  CHF_BOX(a_box));

#if 0
  // do this on WHalf after it's computed and returned, so don't
  // do it here

  // this is a bit of a kluge -- basic idea is that this is
  // where positivity is breaking. so, test to see if
  // computed a_AdWdx will break positivity. if it will, then
  // limit AdWdx to enforce positivity
  FORT_MAPPEDADVPOSITIVITYFIX(CHF_FRA(a_AdWdx),
                              CHF_CONST_FRA(a_W),
                              CHF_BOX(a_box));
#endif

}

// Compute the primitive variables from the conserved variables
// for mapped advection, this means dividing by J
void MappedAdvectionPhysics::consToPrim(FArrayBox&       a_W,
                                        const FArrayBox& a_U,
                                        const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  // convert from conserved to primitive by dividing by J
  a_W.copy(a_U, a_box);
  const FArrayBox& J = *m_JPtr;
  for (int comp=0; comp<a_W.nComp(); comp++)
    {
      a_W.divide(J, 0, comp, 1);
    }
}

// Interval within the primitive variables corresponding to the velocities
// since this doesn't apply to this set of equations, return a bogus interval
Interval MappedAdvectionPhysics::velocityInterval()
{
  MayDay::Error("MappedAdvectionPhysics::velocityInterval - not defined");

  Interval retval(-1,-1);
  return retval;
}

// Component index within the primitive variables of the pressure
// since this doesn't apply to this set of equations, return a bogus value
int MappedAdvectionPhysics::pressureIndex()
{
  MayDay::Error("MappedAdvectionPhysics::pressureIndex - not defined");

  return -1;
}

// Used to limit the absolute value of a "pressure" difference
// since this doesn't apply to this set of equations, return a bogus value
Real MappedAdvectionPhysics::smallPressure()
{
  MayDay::Error("MappedAdvectionPhysics::smallPressure - not defined");

  return -1.0;
}

// Component index within the primitive variables of the bulk modulus
// since this doesn't apply to this set of equations, return a bogus value
int MappedAdvectionPhysics::bulkModulusIndex()
{
  MayDay::Error("MappedAdvectionPhysics::bulkModulusIndex - not defined");

  return -1;
}

// Set face-centered advection velocity pointer
void
MappedAdvectionPhysics::setNTvelPtr(FluxBox* a_NTvelPtr)
{
  // make sure this has at least one component
  CH_assert(a_NTvelPtr->nComp() >= 1);
  m_NTvelPtr = a_NTvelPtr;
}

// Get face-centered advection velocity pointer
const FluxBox*
MappedAdvectionPhysics::NTvelPtr() const
{
  return m_NTvelPtr;
}


// Set cell-centered advection velocity pointer
void
MappedAdvectionPhysics::setCellNTvelPtr(FArrayBox* a_cellNTvelPtr)
{
  m_cellNTvelPtr = a_cellNTvelPtr;
}

// Get cell-centered advection velocity pointer
const FArrayBox*
MappedAdvectionPhysics::cellNTvelPtr() const
{
  return m_cellNTvelPtr;
}


void
MappedAdvectionPhysics::setMappingParameters(const FArrayBox* a_JPtr,
                                             const FluxBox* a_NPtr,
                                             const FluxBox* a_NJInversePtr)
{
  // need to do a const_cast in order to make this work.
  // should be OK, since these guys should be effectively const
  m_JPtr = const_cast<FArrayBox*>(a_JPtr);
  m_NPtr = const_cast<FluxBox*>(a_NPtr);
  m_NJInversePtr = const_cast<FluxBox*>(a_NJInversePtr);
}
