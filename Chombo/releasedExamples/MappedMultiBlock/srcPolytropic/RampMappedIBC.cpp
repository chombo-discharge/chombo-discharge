#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SchwarzChristoffelRampCS.H"
#include "RampMappedIBC.H"
#include "RampMappedIBCF_F.H"
#include "LoHiSide.H"
#include "Box.H"
#include "CHArray.H"
// #include "LoHiCenter.H"  // needed for setBdrySlopes
#include "SetCentersF_F.H"
#include "MOLPhysicsMappedArtViscF_F.H"
#include "PolytropicPhysicsF_F.H"
#include "UnitNormalsF_F.H"


/*******************************************************************************
 *
 *  Class RampMappedIBC: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*
 *  Definitions of static member data
 *--------------------------------------------------------------------*/

bool RampMappedIBC::s_isFortranCommonSet = false;

/*--------------------------------------------------------------------*/
//  Null constructor
/*--------------------------------------------------------------------*/

RampMappedIBC::RampMappedIBC()
  :
  m_alpha(0.),
  m_X0Start(0.)
{
}

/*--------------------------------------------------------------------*/
//  Constructor which defines parameters used by Fortran routines
/**
 *  \param[out] a_smallPressure
 *                     Lower limit for pressure (returned)
 *  \param[in]  a_gamma
 *                     Gamma for polytropic, gamma-law gas
 *  \param[in]  a_alpha
 *                     Ramp angle in radians
 *  \param[in]  a_ms   Mach shock number
 *  \param[in]  a_xStart
 *                     Starting x-dir location of the shock in
 *                     physical space.
 *  \param[in]  a_artvisc
 *                     Artificial viscosity coefficient
 *//*-----------------------------------------------------------------*/

RampMappedIBC::RampMappedIBC(Real&      a_smallPressure,
                             const Real a_gamma,
                             const Real a_alpha,
                             const Real a_ms,
                             const Real a_X0Start,
                             const Real a_artvisc)
  :
  m_alpha(0.),
  m_X0Start(0.)
{
  setFortranCommon(a_smallPressure,
                   a_gamma,
                   a_alpha,
                   a_ms,
                   a_X0Start,
                   a_artvisc);
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/

RampMappedIBC::~RampMappedIBC()
{
}

/*--------------------------------------------------------------------*/
//  Sets parameters in a common block used by Fortran routines
/**
 *  \param[out] a_smallPressure
 *                     Lower limit for pressure (returned)
 *  \param[in]  a_gamma
 *                     Gamma for polytropic, gamma-law gas
 *  \param[in]  a_alpha
 *                     Ramp angle in radians
 *  \param[in]  a_ms   Mach shock number
 *  \param[in]  a_X0Start
 *                     Starting x-dir location of the shock in
 *                     physical space.
 *  \param[in]  a_artvisc
 *                     Artificial viscosity coefficient
 *//*-----------------------------------------------------------------*/

void
RampMappedIBC::setFortranCommon(Real&      a_smallPressure,
                                const Real a_gamma,
                                const Real a_alpha,
                                const Real a_ms,
                                const Real a_X0Start,
                                const Real a_artvisc)
{
  FORT_RAMPMAPPEDSETF(CHF_REAL(a_smallPressure),
                      CHF_CONST_REAL(a_gamma),
                      CHF_CONST_REAL(a_alpha),
                      CHF_CONST_REAL(a_ms),
                      CHF_CONST_REAL(a_X0Start),
                      CHF_CONST_REAL(a_artvisc));

  s_isFortranCommonSet = true;
  m_alpha = a_alpha;
  m_X0Start = a_X0Start;
}

/*--------------------------------------------------------------------*/
//  Factory method - this object is its own factory
/** Return a pointer to a new PhysIBC object
 *//*-----------------------------------------------------------------*/

PhysMappedIBC* RampMappedIBC::new_physIBC()
{
  //**FIXME Use a copy constructor here
  RampMappedIBC* retRamp = new RampMappedIBC();

  if (m_haveTime) retRamp->setTime(m_time);
  if (m_haveCoordSys) retRamp->setCoordSys(m_coordSysPtr);

  retRamp->m_alpha = m_alpha;
  retRamp->m_X0Start = m_X0Start;
  return static_cast<PhysMappedIBC*>(retRamp);
}

/*--------------------------------------------------------------------*/
//  Set up initial conditions
/** \param[out] a_U     Initial conserved state
 *//*-----------------------------------------------------------------*/

void
RampMappedIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(s_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);

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

      // Make sure the CoordSys has the same parameters
      const SchwarzChristoffelRampCS* SCRampCS =
        static_cast<const SchwarzChristoffelRampCS*>(coordSys);
      CH_assert(m_alpha == SCRampCS->rampAngle());

      // Make sure the shock starts in the domain
      CH_assert(SCRampCS->physicalDomainX0Min() < m_X0Start);

      // Xi:  cartesian coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      FORT_SETCELLCENTERS(CHF_FRA(XiFab),
                          CHF_CONST_REAL(m_dx),
                          CHF_BOX(uBox));

      // X:  physical coordinates
      FArrayBox XFab(uBox, SpaceDim);
      coordSys->realCoord(XFab, XiFab, uBox);

      // Set up initial condition in this grid
      FORT_RAMPMAPPEDINITF(CHF_FRA(UFab),
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

void RampMappedIBC::primBC(FArrayBox&             a_WGdnv,
                           const FArrayBox&       a_Wextrap,
                           const FArrayBox&       a_W,
                           const FArrayBox *const a_unitNormalBasisPtr,
                           const Interval&        a_velIntv,
                           const int&             a_dir,
                           const Side::LoHiSide&  a_side,
                           const Real&            a_time)
{
  CH_assert(s_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);

  // Neither the x or y direction can be periodic
  if ((a_dir == 0 || a_dir == 1) && m_domain.isPeriodic(a_dir))
  {
    MayDay::Error("RampMappedIBC::primBC: Neither the x or y boundaries can be "
                  "periodic");
  }

  Box boundaryBox;
  getBoundaryFaces(boundaryBox, a_WGdnv.box(), a_dir, a_side);

  if (! boundaryBox.isEmpty() )
    {
      // Extract the coordinate system for this box.
      Box interiorBox = adjCellBox(boundaryBox, a_dir, Side::flip(a_side), 1);
      const NewCoordSys* coordSys = m_coordSysPtr->getCoordSys(interiorBox);
      CH_assert(coordSys != NULL);

      // Xi:  cartesian coordinates (face centers)
      FArrayBox XiFab(boundaryBox, SpaceDim);
      FORT_SETFACECENTERS(CHF_FRA(XiFab),
                          CHF_CONST_REAL(m_dx),
                          CHF_CONST_INT(a_dir),
                          CHF_BOX(boundaryBox));

      // X:  physical coordinates (face centers)
      FArrayBox XFab(boundaryBox, SpaceDim);
      coordSys->realCoord(XFab, XiFab, boundaryBox);

      // Temporary storage
      FArrayBox WFaceTmp(boundaryBox, a_WGdnv.nComp());

      // Set the boundary fluxes
      int lohisign = sign(a_side);
      FORT_RAMPMAPPEDBCF(CHF_FRA(a_WGdnv),
                         CHF_CONST_FRA(a_Wextrap),
                         CHF_CONST_FRA(WFaceTmp),
                         CHF_CONST_FRA(XFab),
                         CHF_CONST_REAL(a_time),
                         CHF_CONST_INT(lohisign),
                         CHF_CONST_REAL(m_dx),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(boundaryBox));

//--Set a_WGdnv as the solution to the Riemann problem between a_Wextrap and
//--WFaceTmp wherever desired

      if (a_dir == 1 && a_side == Side::Hi)
        {
          // a_Wextrap is left and WFaceTmp is right
          const FArrayBox &unitNormalBasisFab = *a_unitNormalBasisPtr;
          CH_assert(unitNormalBasisFab.box().contains(boundaryBox));

          // Forward transform
          // We want to modify a_Wextrap so make a copy
          FArrayBox WLeft(boundaryBox, a_Wextrap.nComp());
          WLeft.copy(a_Wextrap);
          FArrayBox velLeftFab(a_velIntv, WLeft);  // Alias
          FORT_FORWARDTRANSFORMF(CHF_FRA(velLeftFab),
                                 CHF_CONST_FRA(unitNormalBasisFab),
                                 CHF_BOX(boundaryBox));
          FArrayBox velRightFab(a_velIntv, WFaceTmp);  // Alias
          FORT_FORWARDTRANSFORMF(CHF_FRA(velRightFab),
                                 CHF_CONST_FRA(unitNormalBasisFab),
                                 CHF_BOX(boundaryBox));
          // Riemann problem
          FORT_RIEMANNF(CHF_FRA(a_WGdnv),
                        CHF_CONST_FRA(WLeft),
                        CHF_CONST_FRA(WFaceTmp),
                        CHF_CONST_INT(a_dir),
                        CHF_BOX(boundaryBox));
          // Reverse transform
          FArrayBox velNewFab(a_velIntv, a_WGdnv);  // Alias
          FORT_REVERSETRANSFORMF(CHF_FRA(velNewFab),
                                 CHF_CONST_FRA(unitNormalBasisFab),
                                 CHF_BOX(boundaryBox));
        }
    }
}

/*--------------------------------------------------------------------*/
//  Set boundary slopes
/** The boundary slopes in a_dW are already set to one sided
 *  difference approximations.  If this function doesn't change them
 *  they will be used for the slopes at the boundaries.
 *  \param[out] a_dW
 *  \param[in]  a_W
 *  \param[in]  a_dir
 *  \param[in]  a_time
 *//*-----------------------------------------------------------------*/

void RampMappedIBC::setBdrySlopes(FArrayBox&       a_dW,
                                  const FArrayBox& a_W,
                                  const int&       a_dir,
                                  const Real&      a_time)
{
  CH_assert(s_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  MayDay::Error("RampMappedIBC::setBdrySlopes: is not expected for fourth-"
                  "order solutions");
}

/*--------------------------------------------------------------------*/
//  Adjust boundary fluxes to account for artificial viscosity
/** \param[in]  a_NtFdir
 *                      All components at indices in a_loFaceBox and
 *                      a_hiFaceBox must be set to zero
 *  \param[out] a_NtFdir
 *                      Fluxes due to artificial viscosity added to
 *                      'a_dir' faces along the domain boundary
 *  \param[in]  a_Nctg  N, with components stored contiguously, on the
 *                      'a_dir' faces
 *  \param[in]  a_U     Conservative state \<U\> in the cells
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_divVel
 *                      Physical space divergence of the velocity on
 *                      the 'a_dir' faces
 *  \param[in]  a_csq   Min of c^2/gamma (over +1 layer of cells in
 *                      all transverse directions), stored on cells
 *  \param[in]  a_dxFace
 *                      Distance between cell centers in direction
 *                      'a_dir' stored on faces
 *  \param[in]  a_momIntv
 *                      Interval of momentum components
 *  \param[in]  a_alpha Artificial viscosity coefficient
 *  \param[in]  a_beta  Artificial viscosity coefficient
 *  \param[in]  a_loFaceBox
 *                      Domain faces on the low side
 *  \param[in]  a_hasLo Low-side box exists
 *  \param[in]  a_hiFaceBox
 *                      Domain faces on the high side
 *  \param[in]  a_hasHi High-side box exists
 *  \param[in]  a_dir   Direction of the faces
 *//*-----------------------------------------------------------------*/

void RampMappedIBC::artViscBC(
  FArrayBox&                                   a_NtFdir,
  const CHArray<Real, SpaceDim+1, ArRangeCol>& a_Nctg,
  const FArrayBox&                             a_U,
  const FArrayBox&                             a_unitNormalBasis,
  const FArrayBox&                             a_divVel,
  const FArrayBox&                             a_csq,
  const FArrayBox&                             a_dxFace,
  const Interval&                              a_momIntv,
  const Real                                   a_alpha,
  const Real                                   a_beta,
  const Box&                                   a_loFaceBox,
  const int                                    a_hasLo,
  const Box&                                   a_hiFaceBox,
  const int                                    a_hasHi,
  const int                                    a_dir)
{
  CH_assert(s_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);

  // Just for shifting
  FArrayBox &csq = (FArrayBox &)a_csq;
  FArrayBox &U = (FArrayBox &)a_U;

  if (a_dir == 1)
    {
      if (a_hasLo)
        {
          // All cell centered quantities need to be shifted to the faces
          csq.shiftHalf(a_dir, -1);
          U.shiftHalf(a_dir, -1);

          CH_assert(a_NtFdir.contains(a_loFaceBox));
          CH_assert(U.contains(a_loFaceBox));
          CH_assert(a_unitNormalBasis.contains(a_loFaceBox));
          CH_assert(a_divVel.contains(a_loFaceBox));
          CH_assert(csq.contains(a_loFaceBox));
          CH_assert(a_dxFace.contains(a_loFaceBox));

          // Alias to the momemntum components
          FArrayBox NtFdirMom(a_momIntv, a_NtFdir);  // Alias
          // For U, copy the momentum components so it can be altered
          FArrayBox UMom(a_loFaceBox, SpaceDim);
          UMom.copy(a_U, a_momIntv.begin(), 0, SpaceDim);

          // From MOLPhysicsMappedArtViscF.ChF
          FORT_MAPPEDARTVISCINVISCIDWALL(
            CHF_FRA(NtFdirMom),
            CHF_FRA(UMom),
            CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, a_Nctg),
            CHF_CONST_FRA(a_unitNormalBasis),
            CHF_CONST_FRA1(a_divVel, 0),
            CHF_CONST_FRA1(csq, 0),
            CHF_CONST_FRA1(a_dxFace, 0),
            CHF_CONST_REAL(a_alpha),
            CHF_CONST_REAL(a_beta),
            CHF_BOX(a_loFaceBox),
            CHF_CONST_INT(a_dir));

          // Shift back csq
          csq.shiftHalf(a_dir, 1);
          U.shiftHalf(a_dir, 1);
        }
    }
}
