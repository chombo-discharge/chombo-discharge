#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TrigMappedIBC.H"
#include "TrigIBCF_F.H"
#include "TrigMappedIBCF_F.H"
#include "LoHiSide.H"
// #include "LoHiCenter.H"  // needed for setBdrySlopes
#include "SetCentersF_F.H"
#include "SolidBCF_F.H"

// Null constructor
TrigMappedIBC::TrigMappedIBC()
{
  m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines
TrigMappedIBC::TrigMappedIBC(Real&       a_smallPressure,
                             const Real& a_gamma,
                             const Real& a_ambientDensity,
                             const Real& a_artvisc)
{
  setFortranCommon(a_smallPressure,
                   a_gamma,
                   a_ambientDensity,
                   a_artvisc);
}

TrigMappedIBC::~TrigMappedIBC()
{
}

// Sets parameters in a common block used by Fortran routines:
//      a_smallPressure  - Lower limit for pressure (returned)
//      a_gamma          - Gamma for polytropic, gamma-law gas
//      a_ambientDensity - Ambient density
//      a_pressure       - If 0, use isentropic pressure
//                         if 1, use the constant pressure of 1.0
//      a_artvisc        - Artificial viscosity coefficient
void TrigMappedIBC::setFortranCommon(Real&       a_smallPressure,
                                     const Real& a_gamma,
                                     const Real& a_ambientDensity,
                                     const Real& a_artvisc)
{
  CH_assert(m_isFortranCommonSet == false);

  int fakePressure = 0;  // Unused
  FORT_TRIGSETF(CHF_REAL(a_smallPressure),
                CHF_CONST_REAL(a_gamma),
                CHF_CONST_REAL(a_ambientDensity),
                CHF_CONST_INT(fakePressure),
                CHF_CONST_REAL(a_artvisc));

  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* TrigMappedIBC::new_physIBC()
{
  TrigMappedIBC* retval = new TrigMappedIBC();
  retval->m_isFortranCommonSet = m_isFortranCommonSet;
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void TrigMappedIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
    {
      // Storage for current grid
      FArrayBox& UFab = a_U[dit()];

      // Box of current grid
      Box uBox = UFab.box();
      uBox &= m_domain;

      // Extract the coordinate system for this box.
      const NewCoordSys* coordSys = m_coordSysPtr->getCoordSys(uBox);
      CH_assert(coordSys != NULL);

      // Xi:  cartesian coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      FORT_SETCELLCENTERS(CHF_FRA(XiFab),
                          CHF_CONST_REAL(m_dx),
                          CHF_BOX(uBox));

      // X:  physical coordinates
      FArrayBox XFab(uBox, SpaceDim);
      coordSys->realCoord(XFab, XiFab, uBox);

      // Set up initial condition in this grid
      FORT_TRIGMAPPEDINITF(CHF_FRA(UFab),
                           CHF_CONST_FRA(XFab),
                           CHF_BOX(uBox));
    }
}

/*--------------------------------------------------------------------*/
//  Set boundary fluxes
/** \param[out] a_WGdnv Face-averaged state to resolve for boundary
 *                      (stored on faces)
 *  \param[in]  a_Wextrap
 *                      Extrapolated face-averaged state set by the
 *                      scheme.  Unless you know something better,
 *                      this is the most accurate estimate of the
 *                      state on the face (stored on faces)
 *  \param[in]  a_W     Cell-averaged in the cell containing the face
 *                      (stored on cells, unused)
 *  \param[in]  a_unitNormalBasisPtr
 *                      A normal basis which is required if a Riemann
 *                      solve is applied at the boundary
 *  \param[in]  a_velIntv
 *                      Interval of velocity components
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_side  Side of the box
 *  \param[in]  a_time  Solution time
 *//*-----------------------------------------------------------------*/

void TrigMappedIBC::primBC(FArrayBox&             a_WGdnv,
                           const FArrayBox&       a_Wextrap,
                           const FArrayBox&       a_W,
                           const FArrayBox *const a_unitNormalBasisPtr,
                           const Interval&        a_velIntv,
                           const int&             a_dir,
                           const Side::LoHiSide&  a_side,
                           const Real&            a_time)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_isFortranCommonSet == true);
  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    MayDay::Error("TrigMappedIBC::primBC only supports periodic BCs");
  }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void TrigMappedIBC::setBdrySlopes(FArrayBox&       a_dW,
                                  const FArrayBox& a_W,
                                  const int&       a_dir,
                                  const Real&      a_time)
{
  MayDay::Error("TrigMappedIBC::setBdrySlopes not defined");
}
