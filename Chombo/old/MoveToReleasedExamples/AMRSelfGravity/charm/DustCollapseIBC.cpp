#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cassert>
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>
using std::ifstream;
using std::ios;

#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "DustCollapseIBC.H"

#include "LGintegrator.H"
#include "CONSTANTS_C.H"

#include "SelfGravF_F.H"
#include "DustCollapseIBCF_F.H"
#include "BoundariesF_F.H"


// Null constructor
DustCollapseIBC::DustCollapseIBC()
{
  m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines
DustCollapseIBC::DustCollapseIBC(const Real&  a_gamma,
                                 const Real&  a_radius,
                                 const Real&  a_density,
                                 const Real&  a_artvisc,
                                 const Real&  a_rsTolerance,
                                 const Real&  a_maxRsIter,
                                 const Real&  a_maxMach)
{
  setFortranCommon(a_gamma,
                   a_radius,
                   a_density,
                   a_artvisc,
                   a_rsTolerance,
                   a_maxRsIter,
                   a_maxMach);
}

void DustCollapseIBC::setFortranCommon(const Real&  a_gamma,
                                       const Real&  a_radius,
                                       const Real&  a_density,
                                       const Real&  a_artvisc,
                                       const Real&  a_rsTolerance,
                                       const Real&  a_maxRsIter,
                                       const Real&  a_maxMach)
{
  CH_assert(m_isFortranCommonSet == false);

  m_radius  = a_radius;
  m_density = a_density;

  FORT_SETSELFGRAV(CHF_CONST_REAL(a_gamma),
                   CHF_CONST_REAL(a_artvisc),
                   CHF_CONST_REAL(a_rsTolerance),
                   CHF_CONST_REAL(a_maxRsIter),
                   CHF_CONST_REAL(a_maxMach));


  const Real p0= pow(m_density*m_radius,2);
  const Real small_r0 = 1.e-4 *m_density;
  const Real small_p0 = small *p0 *1.e1;
  const Real small_u0 = small;
  //  const Real small_s0 = small *p0 *pow(m_density,-a_gamma);
  const Real small_s0 = 1.e-2 *p0 *pow(m_density,-a_gamma);

  FORT_SETFLOORS(CHF_CONST_REAL(small_r0),
                 CHF_CONST_REAL(small_u0),
                 CHF_CONST_REAL(small_p0),
                 CHF_CONST_REAL(small_s0));

  m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_physIBC() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void DustCollapseIBC::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
PhysIBC* DustCollapseIBC::new_physIBC()
{
  DustCollapseIBC* retval = new DustCollapseIBC();

  if (m_isFortranCommonSet == true)
  {
    retval->setFortranCommonSet();
    retval->m_radius  = m_radius;
    retval->m_density = m_density;
  }
  else
  {
    MayDay::Error(" DustCollapseIBC::new_physIBC : Fortran Common not set");
  }

  return static_cast<PhysIBC*>(retval);
}

// Set boundary fluxes
void DustCollapseIBC::primBC(FArrayBox&            a_WGdnv,
                             const FArrayBox&      a_Wextrap,
                             const FArrayBox&      a_W,
                             const int&            a_dir,
                             const Side::LoHiSide& a_side,
                             const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(PhysIBC::m_isDefined == true);

  if ( !PhysIBC::m_domain.isPeriodic(a_dir) )
  {
    // In periodic case, this doesn't do anything
    int lohisign;
    Box WGdnvBox = a_WGdnv.box();
    Box tmp = WGdnvBox;

    // Determine which side and thus shifting directions
    lohisign = sign(a_side);
    tmp.shiftHalf(a_dir,lohisign);

    // Is there a domain boundary next to this grid
    if (!PhysIBC::m_domain.contains(tmp))
    {
      tmp &= PhysIBC::m_domain;

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

      // Set the boundary fluxes
      FORT_INFLOWF(CHF_FRA(a_WGdnv),
                   CHF_CONST_FRA(a_Wextrap),
                   CHF_CONST_FRA(a_W),
                   CHF_CONST_INT(lohisign),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(boundaryBox));
    }
  }
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void DustCollapseIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void DustCollapseIBC::setBdrySlopes(FArrayBox&       a_dW,
                                 const FArrayBox& a_W,
                                 const int&       a_dir,
                                 const Real&      a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(PhysIBC::m_isDefined == true);

  if ( !PhysIBC::m_domain.isPeriodic(a_dir) )
  {
    // In periodic case, this doesn't do anything
    Box loBox,hiBox,centerBox,domain;
    int hasLo,hasHi;
    Box slopeBox = a_dW.box()&PhysIBC::m_domain;

    // Generate the domain boundary boxes, loBox and hiBox, if there are
    // domain boundarys there
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
               slopeBox,PhysIBC::m_domain,a_dir);

    // Set the boundary slopes if necessary
    if ((hasLo != 0) || (hasHi != 0))
    {
      FORT_SLOPEBCSF(CHF_FRA(a_dW),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi));
    }
  }
}

// Set up initial conditions
void DustCollapseIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(PhysIBC::m_isDefined == true);
  CH_assert(m_radius !=zero);

  const Real dx = PhysIBC::m_dx;

  const DisjointBoxLayout& grids = a_U.getBoxes();

  // Iterator of all grids in this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];
    U.setVal(zero);

    // Box of current grid
    Box box = grids.get(dit());
    box &= PhysIBC::m_domain;

    // Set up initial condition in this grid
    FORT_DUSTCOLLAPSE(CHF_FRA(U),
                      CHF_CONST_REAL(m_radius),
                      CHF_CONST_REAL(m_density),
                      CHF_CONST_REAL(dx),
                      CHF_BOX(box));
  }
}
