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
#include "LoHiCenter.H"

#include "RampIBC.H"
#include "RampIBCF_F.H"

// Null constructor
RampIBC::RampIBC()
{
    m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines
RampIBC::RampIBC(Real&       a_smallPressure,
                 const Real& a_gamma,
                 const Real& a_alpha,
                 const Real& a_ms,
                 const Real& a_xcorner,
                 const Real& a_artvisc,
                 const Real& a_riemannSolTol,
                 const int& a_maxRsIter,
                 const Real& a_maxMach)
{
  setFortranCommon(a_smallPressure,
                   a_gamma,
                   a_alpha,
                   a_ms,
                   a_xcorner,
                   a_artvisc,
                   a_riemannSolTol,
                   a_maxRsIter,
                   a_maxMach);
}

RampIBC::~RampIBC()
{
}

// Sets parameters in a common block used by Fortran routines:
//   a_smallPressure - Lower limit for pressure (returned)
//   a_gamma         - Gamma for polytropic, gamma-law gas
//   a_alpha         - The angle from horizontal to the wedge edge
//   a_ms            - Mach shock number
//   a_xcorner       - Position of the corner/tip of the wedge
//   a_artvisc       - Artificial viscosity coefficient
void RampIBC::setFortranCommon(Real&       a_smallPressure,
                               const Real& a_gamma,
                               const Real& a_alpha,
                               const Real& a_ms,
                               const Real& a_xcorner,
                               const Real& a_artvisc,
                               const Real& a_riemannSolTol,
                               const int& a_maxRsIter,
                               const Real& a_maxMach)
{
  FORT_RAMPSETF(CHF_REAL(a_smallPressure),
                CHF_CONST_REAL(a_gamma),
                CHF_CONST_REAL(a_alpha),
                CHF_CONST_REAL(a_ms),
                CHF_CONST_REAL(a_xcorner),
                CHF_CONST_REAL(a_artvisc),
                CHF_CONST_REAL(a_riemannSolTol),
                CHF_CONST_INT(a_maxRsIter),
                CHF_CONST_REAL(a_maxMach));

  m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_physIBC() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void RampIBC::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* RampIBC::new_physIBC()
{
  RampIBC* retval = new RampIBC();

  if (m_isFortranCommonSet == true)
  {
    retval->setFortranCommonSet();
  }

  return static_cast<PhysIBC*>(retval);
}

// Set up initial conditions
void RampIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // Iterator of all grids in this level
  for (DataIterator dit = a_U.dataIterator();
       dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;

    // Set up initial condition in this grid
    FORT_RAMPINITF(CHF_FRA(U),
                  CHF_CONST_REAL(m_dx),
                  CHF_BOX(uBox));
  }
}

// Set boundary fluxes
void RampIBC::primBC(FArrayBox&            a_WGdnv,
                     const FArrayBox&      a_Wextrap,
                     const FArrayBox&      a_W,
                     const int&            a_dir,
                     const Side::LoHiSide& a_side,
                     const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  int isign = sign(a_side);
  Box gdnvBox = a_WGdnv.box();
  Box bdryBox;

  if (a_side == Side::Lo)
  {
    bdryBox = bdryLo(m_domain,a_dir,1);
  }

  if (a_side == Side::Hi)
  {
    bdryBox = bdryHi(m_domain,a_dir,1);
  }

  bdryBox &= gdnvBox;

  // Is there a domain boundary next to this grid
  if (!bdryBox.isEmpty())
  {
    FORT_RAMPBCF(CHF_FRA(a_WGdnv),
                 CHF_CONST_FRA(a_Wextrap),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_REAL(a_time),
                 CHF_CONST_INT(isign),
                 CHF_CONST_REAL(m_dx),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(bdryBox));
  }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void RampIBC::setBdrySlopes(FArrayBox&       a_dW,
                            const FArrayBox& a_W,
                            const int&       a_dir,
                            const Real&      a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box loBox,hiBox,centerBox,domain;
  int hasLo,hasHi;
  Box slopeBox = a_dW.box();
  slopeBox.grow(a_dir,1);

  // Generate the domain boundary boxes, loBox and hiBox, if there are
  // domain boundarys there
  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
             slopeBox,m_domain,a_dir);

  // Set the boundary slopes if necessary
  if ((hasLo != 0) || (hasHi != 0))
  {
    FORT_RAMPSLOPEBCSF(CHF_FRA(a_dW),
                       CHF_CONST_FRA(a_W),
                       CHF_CONST_REAL(m_dx),
                       CHF_CONST_INT(a_dir),
                       CHF_BOX(loBox),
                       CHF_CONST_INT(hasLo),
                       CHF_BOX(hiBox),
                       CHF_CONST_INT(hasHi));
  }
}

// Adjust boundary fluxes to account for artificial viscosity
void RampIBC::artViscBC(FArrayBox&       a_F,
                        const FArrayBox& a_U,
                        const FArrayBox& a_divVel,
                        const int&       a_dir,
                        const Real&      a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  FArrayBox &divVel = (FArrayBox &)a_divVel;

  Box fluxBox = divVel.box();
  fluxBox.enclosedCells(a_dir);
  fluxBox.grow(a_dir,1);

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                 fluxBox,m_domain,a_dir);

  if (hasLo == 1)
  {
    loBox.shiftHalf(a_dir,1);
    divVel.shiftHalf(a_dir,1);
    a_F.shiftHalf(a_dir,1);

    int loHiSide = -1;

    FORT_RAMPARTVISCF(CHF_FRA(a_F),
                      CHF_CONST_FRA(a_U),
                      CHF_CONST_FRA1(divVel,0),
                      CHF_CONST_INT(loHiSide),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(loBox));

    loBox.shiftHalf(a_dir,-1);
    divVel.shiftHalf(a_dir,-1);
    a_F.shiftHalf(a_dir,-1);
  }

  if (hasHi == 1)
  {
    hiBox.shiftHalf(a_dir,-1);
    divVel.shiftHalf(a_dir,-1);
    a_F.shiftHalf(a_dir,-1);

    int loHiSide = 1;

    FORT_RAMPARTVISCF(CHF_FRA(a_F),
                      CHF_CONST_FRA(a_U),
                      CHF_CONST_FRA1(divVel,0),
                      CHF_CONST_INT(loHiSide),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(hiBox));

    hiBox.shiftHalf(a_dir,1);
    divVel.shiftHalf(a_dir,1);
    a_F.shiftHalf(a_dir,1);
  }
}
