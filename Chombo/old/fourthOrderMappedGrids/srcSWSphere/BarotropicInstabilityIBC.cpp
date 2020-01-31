#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "CubedSphere2DPanelCS.H"
#include "SWintegrator.H"
#include "BarotropicInstabilityIBC.H"
#include "BarotropicInstabilityF_F.H"
#include "NamespaceHeader.H"

// Null constructor
BarotropicInstabilityIBC::BarotropicInstabilityIBC()
{
  CH_assert(false);
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
BarotropicInstabilityIBC::BarotropicInstabilityIBC(
                                                   Real a_phi0,
                                                   Real a_phi1,
                                                   Real a_phi2,
                                                   Real a_umax,
                                                   Real a_hmean,
                                                   Real a_hhat,
                                                   Real a_BIalpha,
                                                   Real a_BIbeta)
{
  m_phi0 = a_phi0;
  m_phi1 = a_phi1;
  m_phi2 = a_phi2;
  m_umax = a_umax;
  m_hmean = a_hmean;
  m_hhat = a_hhat;
  m_BIalpha = a_BIalpha;
  m_BIbeta = a_BIbeta;
}

BarotropicInstabilityIBC::~BarotropicInstabilityIBC()
{
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* BarotropicInstabilityIBC::new_physIBC()
{
  BarotropicInstabilityIBC* retval =
    new BarotropicInstabilityIBC(m_phi0, m_phi1, m_phi2, m_umax, m_hmean, m_hhat, m_BIalpha, m_BIbeta);
  // retval->m_isFortranCommonSet = m_isFortranCommonSet;
  if (m_isFortranCommonSet)
    {
      retval->setFortranCommon(m_gravity, m_omega, m_alpha);
    }
/*
  if (m_gotTime) retval->setTime(m_time);
*/
  if (m_gotCoordSys) retval->setCoordSys(m_coordSysPtr);
  return static_cast<PhysMappedIBC*>(retval);
}


void BarotropicInstabilityIBC::setFortranCommon(const Real& a_gravity,
                                                const Real& a_omega,
                                                const Real& a_alpha)
{
  PhysShallowWaterMappedIBC::setFortranCommon(a_gravity, a_omega, a_alpha);
  FORT_BAROTROPICINSTABILITYSETPARAMS(
                                      CHF_CONST_REAL(m_phi0),
                                      CHF_CONST_REAL(m_phi1),
                                      CHF_CONST_REAL(m_phi2),
                                      CHF_CONST_REAL(m_umax),
                                      CHF_CONST_REAL(m_hmean),
                                      CHF_CONST_REAL(m_hhat),
                                      CHF_CONST_REAL(m_BIalpha),
                                      CHF_CONST_REAL(m_BIbeta));
}


// Set up initial conditions
void BarotropicInstabilityIBC::initializeUnified(LevelData<FArrayBox>& a_U,
                                        bool a_multiplyJ)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_gotCoordSys);
  CH_assert(m_gotTime);

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  int nComp = a_U.nComp();
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();

  LevelData<FArrayBox> W(layout, WNUM, a_U.ghostVect());
  LevelData<FArrayBox> Wlonlat(layout, WNUM, a_U.ghostVect());

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];

      // Storage for current grid
      FArrayBox& UFab = a_U[dit];

      // Box of current grid
      Box uBox = UFab.box();
      // uBox &= m_domain;

      const CubedSphere2DPanelCS* coordSysBlockPtr =
        dynamic_cast<const CubedSphere2DPanelCS*>(
          m_coordSysPtr->getCoordSys(baseBox));

      CH_assert(coordSysBlockPtr);

      // For each point:
      // set RealVect Xi, which is just linear,
      // and then RealVect X( m_coordSysPtr->realCoord( Xi ) );
      // and then Real J( m_coordSysPtr->pointwiseJ( X ) );

      // Xi: mapped space coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      coordSysBlockPtr->getCellMappedCoordinates(XiFab, uBox);

      FArrayBox lonlatFab(uBox, SpaceDim);
      coordSysBlockPtr->fabTransformEquiangularToLonLat(XiFab, lonlatFab);

      // Create RLL velocity Fab
      FArrayBox vecRLLFab(uBox, SpaceDim);

      // Create mapped velocity Fab
      FArrayBox velFab(uBox, SpaceDim);

      FORT_BAROTROPICINSTABILITYVELOCITY(CHF_FRA(vecRLLFab),
                                         CHF_CONST_FRA(lonlatFab),
                                         CHF_BOX(uBox));

      // Convert velocities from longitude-latitude space to mapped space
      coordSysBlockPtr->fabVectorTransformLatLonToEquiangular(XiFab,
                                                              vecRLLFab,
                                                              velFab);

      FORT_BAROTROPICINSTABILITYHEIGHT(CHF_FRA1(UFab, 0),
                                       CHF_CONST_FRA(lonlatFab),
                                       CHF_BOX(uBox));

      UFab.copy(UFab, UHGT, UMOMX); // UFab[UMOMX] := UFab[0]
      UFab.copy(UFab, UHGT, UMOMY); // UFab[UMOMY] := UFab[0]
      UFab.mult(velFab, 0, UMOMX); // UFab[UMOMX] *= velFab[0]
      UFab.mult(velFab, 1, UMOMY); // UFab[UMOMY] *= velFab[1]

      // Multiply every component of U by J
      if (a_multiplyJ)
        {
          // Get pointwise J (not cell-averaged <J>).
          FArrayBox JFab(uBox, 1);
          coordSysBlockPtr->pointwiseJ(JFab, XiFab, uBox);
          for (int comp = 0; comp < nComp; comp++)
            {
              UFab.mult(JFab, 0, comp);
            }
        }
    }
    // dummy statement in order to get around gdb bug
    int dummy_unused = 0; dummy_unused = 0;
}

// Set up initial conditions
void BarotropicInstabilityIBC::initialize(LevelData<FArrayBox>& a_U)
{
  // Initialize the flow field
  initializeUnified(a_U, false);
}

// Set up initial conditions with J
void BarotropicInstabilityIBC::initializeWithJ(LevelData<FArrayBox>& a_U)
{
  // Initialize the flow field with J
  initializeUnified(a_U, true);
}

// Set boundary fluxes
void BarotropicInstabilityIBC::primBC(FArrayBox&            a_WGdnv,
                         const FArrayBox&      a_Wextrap,
                         const FArrayBox&      a_W,
                         const int&            a_dir,
                         const Side::LoHiSide& a_side,
                         const Real&           a_time)
{
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void BarotropicInstabilityIBC::setBdrySlopes(FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Real&      a_time)
{
  MayDay::Error("BarotropicInstabilityIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void BarotropicInstabilityIBC::artViscBC(FArrayBox&       a_F,
                            const FArrayBox& a_U,
                            const FArrayBox& a_divVel,
                            const int&       a_dir,
                            const Real&      a_time)
{
}

#include "NamespaceFooter.H"
