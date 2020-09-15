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

#include "GaussianPBC.H"
#include "GaussianPBCF_F.H"

#include "NamespaceHeader.H"

// Null constructor
GaussianPBC::GaussianPBC()
{
  m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines
GaussianPBC::GaussianPBC(Real&           a_smallPressure,
                         const Real&     a_gamma,
                         const Real&     a_ambientDensity,
                         const Real&     a_ambientPressure,
                         const Real&     a_deltaPressure,
                         const RealVect& a_center,
                         const Real&     a_size,
                         const RealVect& a_velocity,
                         const Real&     a_artvisc)
{
  setFortranCommon(a_smallPressure,
                   a_gamma,
                   a_ambientDensity,
                   a_ambientPressure,
                   a_deltaPressure,
                   a_center,
                   a_size,
                   a_velocity,
                   a_artvisc);
}

GaussianPBC::~GaussianPBC()
{
}

// Sets parameters in a common block used by Fortran routines:
//   a_smallPressure  - Lower limit for pressure (returned)
//   a_gamma          - Gamma for polytropic, gamma-law gas
//   a_ambientDensity - Ambient density, constant
//   a_ambientPressure- Ambient pressure to add to the pressure gaussian
//   a_deltaPressure  - Mean of the gaussian
//   a_center         - Center of the gaussian
//   a_size           - Standard deviation of the gaussian
//   a_velocity       - Initial velocity of the gas
//   a_artvisc        - Artificial viscosity coefficient
void GaussianPBC::setFortranCommon(Real&           a_smallPressure,
                                   const Real&     a_gamma,
                                   const Real&     a_ambientDensity,
                                   const Real&     a_ambientPressure,
                                   const Real&     a_deltaPressure,
                                   const RealVect& a_center,
                                   const Real&     a_size,
                                   const RealVect& a_velocity,
                                   const Real&     a_artvisc)
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_GAUSSIANPSETF(CHF_REAL(a_smallPressure),
                     CHF_CONST_REAL(a_gamma),
                     CHF_CONST_REAL(a_ambientDensity),
                     CHF_CONST_REAL(a_ambientPressure),
                     CHF_CONST_REAL(a_deltaPressure),
                     CHF_CONST_REALVECT(a_center),
                     CHF_CONST_REAL(a_size),
                     CHF_CONST_REALVECT(a_velocity),
                     CHF_CONST_REAL(a_artvisc));

  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* GaussianPBC::new_physIBC()
{
  GaussianPBC* retval = new GaussianPBC();
  retval->m_isFortranCommonSet = m_isFortranCommonSet;

  return static_cast<PhysIBC*>(retval);
}

// Set up initial conditions
void GaussianPBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;

    // Set up initial condition in this grid
    FORT_GAUSSIANPINITF(CHF_FRA(U),
                        CHF_CONST_REAL(m_dx),
                        CHF_BOX(uBox));
  }
}

// Set boundary fluxes
void GaussianPBC::primBC(FArrayBox&            a_WGdnv,
                         const FArrayBox&      a_Wextrap,
                         const FArrayBox&      a_W,
                         const int&            a_dir,
                         const Side::LoHiSide& a_side,
                         const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
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

      // Set the boundary fluxes
      FORT_SOLIDBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_Wextrap),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
    }
  }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void GaussianPBC::setBdrySlopes(FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Real&      a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
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

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void GaussianPBC::artViscBC(FArrayBox&       a_F,
                            const FArrayBox& a_U,
                            const FArrayBox& a_divVel,
                            const int&       a_dir,
                            const Real&      a_time)
{
}

#include "NamespaceFooter.H"
