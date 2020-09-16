#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GaussianMappedIBC.H"
#include "GaussianSmoothBCF_F.H"
#include "GaussianMappedIBCF_F.H"
#include "LoHiSide.H"
// #include "LoHiCenter.H"  // needed for setBdrySlopes
#include "SetCentersF_F.H"
#include "SolidBCF_F.H"

// Null constructor
GaussianMappedIBC::GaussianMappedIBC()
  :
  m_isFortranCommonSet(false)
{
  m_haveExactSoln = true;
}

// Constructor which defines parameters used by Fortran routines.
GaussianMappedIBC::GaussianMappedIBC(Real&           a_smallPressure,
                                     const Real&     a_gamma,
                                     const Real&     a_ambientDensity,
                                     const Real&     a_deltaDensity,
                                     const int&      a_pressure,
                                     const RealVect& a_center,
                                     const Real&     a_size,
                                     const Real&     a_radmax,
                                     const RealVect& a_velocity,
                                     const Real&     a_artvisc)
{
  setFortranCommon(a_smallPressure,
                   a_gamma,
                   a_ambientDensity,
                   a_deltaDensity,
                   a_pressure,
                   a_center,
                   a_size,
                   a_radmax,
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
//   a_radmax         - Cutoff radius of the gaussian
//   a_velocity       - Initial velocity of the gas
//   a_artvisc        - Artificial viscosity coefficient
// Mapping changes nothing.
void GaussianMappedIBC::setFortranCommon(Real&           a_smallPressure,
                                         const Real&     a_gamma,
                                         const Real&     a_ambientDensity,
                                         const Real&     a_deltaDensity,
                                         const int&      a_pressure,
                                         const RealVect& a_center,
                                         const Real&     a_size,
                                         const Real&     a_radmax,
                                         const RealVect& a_velocity,
                                         const Real&     a_artvisc)
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_GAUSSIANSMOOTHSETF(CHF_REAL(a_smallPressure),
                    CHF_CONST_REAL(a_gamma),
                    CHF_CONST_REAL(a_ambientDensity),
                    CHF_CONST_REAL(a_deltaDensity),
                    CHF_CONST_INT(a_pressure),
                    CHF_CONST_REALVECT(a_center),
                    CHF_CONST_REAL(a_size),
                    CHF_CONST_REAL(a_radmax),
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
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void GaussianMappedIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_haveCoordSys);
  CH_assert(m_haveTime);

  // Iterator of all grids in this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // Storage for current grid
      FArrayBox& UFab = a_U[dit()];

      // Box of current grid
      Box uBox = UFab.box();
      uBox &= m_domain;

      // Extract the coordinate system for this box.
      Box box = a_U.disjointBoxLayout()[dit()];
      const NewCoordSys* coordSys = m_coordSysPtr->getCoordSys(box);
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
      FORT_GAUSSIANMAPPEDINITF(CHF_FRA(UFab),
                               CHF_CONST_FRA(XFab),
                               CHF_BOX(uBox),
                               CHF_CONST_REAL(m_time));

//**You can't do this because gradients of pointwise JU may violate free-stream
//**preservation.  So this just returns pointwise U.  In the calling program,
//**it is better to find <U> and then <JU> using the product rule.
//       // Get pointwise J (not cell-averaged <J>).
//       FArrayBox JFab(uBox, 1);
//       // m_coordSysPtr->pointwiseJ(JFab);
//       m_coordSysPtr->pointwiseJ(JFab, XiFab, uBox);

//       // Multiply every component of U by J.
//       for (int comp = 0; comp < nComp; comp++)
//         {
//           UFab.mult(JFab, 0, comp);
//         }
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

void GaussianMappedIBC::primBC(FArrayBox&             a_WGdnv,
                               const FArrayBox&       a_Wextrap,
                               const FArrayBox&       a_W,
                               const FArrayBox *const a_unitNormalBasisPtr,
                               const Interval&        a_velIntv,
                               const int&             a_dir,
                               const Side::LoHiSide&  a_side,
                               const Real&            a_time)
{
  CH_assert(m_isDefined);
  CH_assert(m_isFortranCommonSet);
  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    MayDay::Error("GaussianMappedIBC::primBC only supports periodic BCs");
  }
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
}
