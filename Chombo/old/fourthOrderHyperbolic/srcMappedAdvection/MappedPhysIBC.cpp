#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include "MappedPhysIBC.H"


/// Constructor
/**
 */
MappedPhysIBC::MappedPhysIBC()
{
}

/// Destructor
/**
 */
MappedPhysIBC::~MappedPhysIBC()
{
  if (m_basicIBCptr != NULL)
    {
      delete m_basicIBCptr;
    }
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Note that in this context,
   a_dx is really dXi
*/
void
MappedPhysIBC::define(const ProblemDomain& a_domain,
                      const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a point to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/

PhysIBC*
MappedPhysIBC::new_physIBC()
{
  MappedPhysIBC* newIBC = new MappedPhysIBC;
  newIBC->define(m_domain, m_dx);

  newIBC->setBasicIBC(m_basicIBCptr);

  newIBC->setCoordSys(m_coordSysPtr);
  return static_cast<PhysIBC*>(newIBC);

}

/// Set up initial conditions
/**
 */
void
MappedPhysIBC::initialize(LevelData<FArrayBox>& a_U)
{
  Real time = 0.0;
  m_basicIBCptr->initialize(a_U, m_domain, *m_coordSysPtr,
                            m_dx, time);
}

/// Set boundary fluxes
/**
 */
void
MappedPhysIBC::primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time)
{
  if (!m_domain.isPeriodic(a_dir))
    {
      // don't need to do 4th-order BC's, so we can use point values here,
      // (I think)
      int lohisign;
      Box tmp = a_WGdnv.box();

      // Determine which side and thus shifting directions
      lohisign = sign(a_side);
      tmp.shiftHalf(a_dir,lohisign);

      // Is there a domain boundary next to this grid
      if (!m_domain.contains(tmp))
        {
          tmp &= m_domain;

          Box boundaryBox;

          // Find the strip of cells next to the domain boundary
          if (a_side == Side::Lo)
            {
              boundaryBox = bdryLo(tmp,a_dir);
            }
          else
            {
              boundaryBox = bdryHi(tmp,a_dir);
            }

          bool includeJ = false;
          // note that boundaryBox is face-centered!
          m_basicIBCptr->pointVal(a_WGdnv, m_domain, *m_coordSysPtr,
                                  boundaryBox, includeJ, m_dx, a_time);
        }
    }
}

/// Set boundary slopes
/**
   The boundary slopes in a_dW are already set to one sided difference
   approximations.  If this function doesn't change them they will be
   used for the slopes at the boundaries.
*/

void
MappedPhysIBC::setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // bail on this for now -- will implement later
  if (!m_domain.isPeriodic(a_dir))
    {
      MayDay::Error("MappedPhysIBC::setBdrySlopes not implemented yet");
    }


}

/// Adjust boundary fluxes to account for artificial viscosity
/**
 */

void
MappedPhysIBC::artViscBC(FArrayBox&       a_F,
                         const FArrayBox& a_U,
                         const FArrayBox& a_divVel,
                         const int&       a_dir,
                         const Real&      a_time)
{
  // bail on this for now -- will implement later
  if (!m_domain.isPeriodic(a_dir))
    {
      MayDay::Error("MappedPhysIBC::artViscBC not implemented yet");
    }

}

/// set mapped-grid basicIBC member
void
MappedPhysIBC::setBasicIBC(BasicIBC* a_basicIBCptr)
{
  m_basicIBCptr = a_basicIBCptr->new_basicIBC();
}


/// set coordinate-system member
void
MappedPhysIBC::setCoordSys(CoordSys<FArrayBox,FluxBox>* a_coordSysPtr)
{
  m_coordSysPtr = a_coordSysPtr;
}
