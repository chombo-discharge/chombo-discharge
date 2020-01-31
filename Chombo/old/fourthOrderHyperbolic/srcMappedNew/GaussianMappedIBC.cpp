#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// needed for setBdrySlopes
// #include "LoHiCenter.H"

#include "GaussianMappedIBC.H"

#include "GaussianSmoothBCF_F.H"
#include "GaussianMappedIBCF_F.H"
#include "SetCellCentersF_F.H"
#include "SolidBCF_F.H"

// Null constructor
GaussianMappedIBC::GaussianMappedIBC()
{
  m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
GaussianMappedIBC::GaussianMappedIBC(Real&           a_smallPressure,
                                     const Real&     a_gamma,
                                     const Real&     a_ambientDensity,
                                     const Real&     a_deltaDensity,
                                     const int&      a_pressure,
                                     const RealVect& a_center,
                                     const RealVect& a_width, // for smoothing
                                     const Real&     a_size,
                                     const RealVect& a_velocity,
                                     const Real&     a_artvisc)
{
  setFortranCommon(a_smallPressure,
                   a_gamma,
                   a_ambientDensity,
                   a_deltaDensity,
                   a_pressure,
                   a_center,
                   a_width,
                   a_size,
                   a_velocity,
                   a_artvisc);
}

GaussianMappedIBC::~GaussianMappedIBC()
{
}

// Sets parameters in a common block used by Fortran routines:
//   a_smallPressure  - Lower limit for pressure (returned)
//   a_gamma          - Gamma for polytropic, gamma-law gas
//   a_ambientDensity - Ambient density add to the density gaussian
//   a_deltaDensity   - Mean of the gaussian
//   a_pressure       - If 0, use isentropic pressure
//                      if 1, use the constant pressure of 1.0
//   a_center         - Center of the gaussian
//   a_size           - Standard deviation of the gaussian
//   a_velocity       - Initial velocity of the gas
//   a_artvisc        - Artificial viscosity coefficient
// Mapping changes nothing.
void GaussianMappedIBC::setFortranCommon(Real&           a_smallPressure,
                                         const Real&     a_gamma,
                                         const Real&     a_ambientDensity,
                                         const Real&     a_deltaDensity,
                                         const int&      a_pressure,
                                         const RealVect& a_center,
                                         const RealVect& a_width,
                                         const Real&     a_size,
                                         const RealVect& a_velocity,
                                         const Real&     a_artvisc)
{
  CH_assert(m_isFortranCommonSet == false);
  // no change from Chombo/example/AMRGodunov/srcPolytropic/
  FORT_GAUSSIANSMOOTHSETF(CHF_REAL(a_smallPressure),
                    CHF_CONST_REAL(a_gamma),
                    CHF_CONST_REAL(a_ambientDensity),
                    CHF_CONST_REAL(a_deltaDensity),
                    CHF_CONST_INT(a_pressure),
                    CHF_CONST_REALVECT(a_center),
                    CHF_CONST_REALVECT(a_width), // for smoothing
                    CHF_CONST_REAL(a_size),
                    CHF_CONST_REALVECT(a_velocity),
                    CHF_CONST_REAL(a_artvisc));

  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* GaussianMappedIBC::new_physIBC()
{
  GaussianMappedIBC* retval = new GaussianMappedIBC();
  retval->m_isFortranCommonSet = m_isFortranCommonSet;
  if (m_gotTime) retval->setTime(m_time);
  if (m_gotCoordSys) retval->setCoordSys(m_coordSysPtr);

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void GaussianMappedIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_gotCoordSys);
  CH_assert(m_gotTime);
  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  int nComp = a_U.nComp();
  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
    {
      // Storage for current grid
      FArrayBox& UFab = a_U[dit()];

      // Box of current grid
      Box uBox = UFab.box();
      uBox &= m_domain;

      // For each point:
      // set RealVect Xi, which is just linear,
      // and then RealVect X( m_coordSysPtr->realCoord( Xi ) );
      // and then Real J( m_coordSysPtr->pointwiseJ( X ) );

      // Xi:  cartesian coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      FORT_SETCELLCENTERS(CHF_FRA(XiFab),
                          CHF_CONST_REAL(m_dx),
                          CHF_BOX(uBox));

      // X:  physical coordinates
      FArrayBox XFab(uBox, SpaceDim);
      // m_coordSysPtr->realCoord(XFab, XiFab);
      m_coordSysPtr->realCoord(XFab, XiFab, uBox);
      // Set up initial condition in this grid
      FORT_GAUSSIANMAPPEDINITF(CHF_FRA(UFab),
                               CHF_CONST_FRA(XFab),
                               CHF_BOX(uBox));
      // Get pointwise J (not cell-averaged <J>).
      FArrayBox JFab(uBox, 1);
      // m_coordSysPtr->pointwiseJ(JFab);
      m_coordSysPtr->pointwiseJ(JFab, XiFab, uBox);

      // Multiply every component of U by J.
      for (int comp = 0; comp < nComp; comp++)
        {
          UFab.mult(JFab, 0, comp);
        }
  }
}

// Set boundary fluxes
void GaussianMappedIBC::primBC(FArrayBox&            a_WGdnv,
                         const FArrayBox&      a_Wextrap,
                         const FArrayBox&      a_W,
                         const int&            a_dir,
                         const Side::LoHiSide& a_side,
                         const Real&           a_time)
{
  CH_assert(m_isDefined);
  CH_assert(m_isFortranCommonSet);
#if 0
  // petermc, 16 Apr 2009:  kludge  :(
  // These two assertions fail.
  CH_assert(m_gotCoordSys);
  CH_assert(m_gotTime);

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
#endif
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void GaussianMappedIBC::setBdrySlopes(FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Real&      a_time)
{
  MayDay::Error("GaussianMappedIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void GaussianMappedIBC::artViscBC(FArrayBox&       a_F,
                            const FArrayBox& a_U,
                            const FArrayBox& a_divVel,
                            const int&       a_dir,
                            const Real&      a_time)
{
}
