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
#include "RossbyHaurwitzIBC.H"
#include "NamespaceHeader.H"

// Null constructor
RossbyHaurwitzIBC::RossbyHaurwitzIBC()
{
  CH_assert(false);
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
RossbyHaurwitzIBC::RossbyHaurwitzIBC(
                                     Real a_w,
                                     Real a_K,
                                     Real a_h0,
                                     int a_R)
{
  m_w = a_w;
  m_K = a_K;
  m_h0 = a_h0;
  m_R = a_R;
}

RossbyHaurwitzIBC::~RossbyHaurwitzIBC()
{
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* RossbyHaurwitzIBC::new_physIBC()
{
  RossbyHaurwitzIBC* retval = new RossbyHaurwitzIBC(m_w, m_K, m_h0, m_R);
  // retval->m_isFortranCommonSet = m_isFortranCommonSet;
  if (m_isFortranCommonSet) retval->setFortranCommon(m_gravity, m_omega, m_alpha);
/*
  if (m_gotTime) retval->setTime(m_time);
*/
  if (m_gotCoordSys) retval->setCoordSys(m_coordSysPtr);
  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void RossbyHaurwitzIBC::initializeUnified(LevelData<FArrayBox>& a_U,
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

      // Evaluate flow field
      BoxIterator bit(uBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          // Get lonlat coordinate
          // RealVect xi;
          // xi[0] = XiFab(iv,0);
          // xi[1] = XiFab(iv,1);
          // RealVect lonlat;
          // coordSysBlockPtr->pointTransformEquiangularToLonLat(xi, lonlat);
          Real lon = lonlatFab(iv, 0); // lonlat[0];
          Real lat = lonlatFab(iv, 1); // lonlat[1];

          // Calculate fluid height, Real dH.
          Real cos2 = cos(lat)*cos(lat);
          Real cos2R = pow(cos2, m_R);
          Real cosR = pow(cos(lat), m_R);
          Real A = 0.5*m_w * (2*m_omega + m_w) * cos2 +
            0.25 * m_K*m_K * cos2R *
            ((m_R + 1) * cos2 + (m_R*(2*m_R-1)-2)*1. - (2*m_R*m_R)*(1./cos2));
          Real B = ( 2.*(m_omega + m_w)*m_K / ((m_R+1)*(m_R+2)*1.) ) *
            cosR * ((m_R*(m_R+2)+2)*1. - ((m_R+1)*(m_R+1))*cos2);
          Real C = 0.25*m_K*m_K * cos2R * ((m_R+1)*cos2 - (m_R+2)*1.);

          Real dH = m_h0 + (A + B*cos(m_R*lon) + C*cos((2*m_R)*lon)) / m_gravity;

          Real cosRm1 = pow(cos(lat), m_R-1);
          Real dUlon = m_w * cos(lat) +
            m_K * cosRm1 * (m_R*sin(lat)*sin(lat) - cos(lat)*cos(lat)) *
            cos(m_R*lon);
          Real dUlat = -m_R * m_K * cosRm1 * sin(lat) * sin(m_R*lon);
          // FIXME!
          // dH = m_h0;
          // dUlon = cos(lat);
          // dUlat = 0.;

          UFab(iv,0) = dH;

          // Calculate fluid velocities in longitude-latitude coordinates
          vecRLLFab(iv,0) = dUlon;
          vecRLLFab(iv,1) = dUlat;
        }

      // Convert velocities from longitude-latitude space to mapped space
      coordSysBlockPtr->fabVectorTransformLatLonToEquiangular(XiFab,
                                                              vecRLLFab,
                                                              velFab);

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
void RossbyHaurwitzIBC::initialize(LevelData<FArrayBox>& a_U)
{
  // Initialize the flow field
  initializeUnified(a_U, false);
}

// Set up initial conditions with J
void RossbyHaurwitzIBC::initializeWithJ(LevelData<FArrayBox>& a_U)
{
  // Initialize the flow field with J
  initializeUnified(a_U, true);
}

// Set boundary fluxes
void RossbyHaurwitzIBC::primBC(FArrayBox&            a_WGdnv,
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
void RossbyHaurwitzIBC::setBdrySlopes(FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Real&      a_time)
{
  MayDay::Error("RossbyHaurwitzIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void RossbyHaurwitzIBC::artViscBC(FArrayBox&       a_F,
                            const FArrayBox& a_U,
                            const FArrayBox& a_divVel,
                            const int&       a_dir,
                            const Real&      a_time)
{
}

#include "NamespaceFooter.H"
