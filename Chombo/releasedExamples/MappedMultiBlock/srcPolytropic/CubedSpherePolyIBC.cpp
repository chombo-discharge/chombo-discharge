#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSpherePolyIBC.H"
#include "CubedSpherePolyIBCF_F.H"
#include "LoHiSide.H"
#include "Box.H"
#include "RampMappedIBCF_F.H"
#include "CHArray.H"
#include "BoxIterator.H"
#include "CubedSphereShellPanelCS.H"
#include "LGintegrator.H"
#include "LoHiCenter.H"
#include "SetCentersF_F.H"
#include "MOLPhysicsMappedArtViscF_F.H"
#include "PolytropicPhysicsF_F.H"
#include "UnitNormalsF_F.H"


/*******************************************************************************
 *
 *  Class CubedSpherePolyIBC: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*
 *  Definitions of static member data
 *--------------------------------------------------------------------*/

bool CubedSpherePolyIBC::s_isFortranCommonSet = false;

/*--------------------------------------------------------------------*/
//  Null constructor
/*--------------------------------------------------------------------*/
// Return the advection velocity at each cell center
void 
CubedSpherePolyIBC::
getAdvVel(LevelData<FArrayBox>& a_advVel, Real a_time)
{
  CH_TIME("CubedSphereFlowIBC::advVel");
  CH_assert(m_haveCoordSys);

  // For this particular IBC, a_advVel will be independent of a_time.
  const DisjointBoxLayout& layout = a_advVel.disjointBoxLayout();
  DataIterator dit = a_advVel.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];

      // Storage for current grid
      FArrayBox& velFab = a_advVel[dit];

      // Box of current grid
      Box uBox = velFab.box();

      // removed by petermc, 9 Feb 2011
      // uBox &= m_domain;

      const CubedSphereShellPanelCS* coordSysBlockPtr =
        dynamic_cast<const CubedSphereShellPanelCS*>(
          m_coordSysPtr->getCoordSys(baseBox));

      if(coordSysBlockPtr == NULL)
        {
          MayDay::Error("this ibc was kinda hardwired for CubedSphereShellPanelCS");
        }
      // Xi: mapped space coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      coordSysBlockPtr->getCellMappedCoordinates(XiFab, uBox);

      pointAdvVel(velFab, XiFab, a_time, *coordSysBlockPtr);
    }

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_advVel);
}


// Return the advection velocity at each point
void CubedSpherePolyIBC::pointAdvVel(FArrayBox& a_velFab,
                                     const FArrayBox& a_xiLoc, Real a_time, 
                                     const CubedSphereShellPanelCS& coordSysBlock)
{
    // rll: longitude-latitude coordinates
    Box velBox = a_xiLoc.box();
    FArrayBox rllFab(velBox, SpaceDim);
    coordSysBlock.fabTransformEquiangularToLonLat(a_xiLoc, rllFab);

    // vecRLL: vector in longitude-latitude basis
    FArrayBox vecRLLFab(velBox, SpaceDim);

    // Constant in time
    FORT_CSPVECLONLATSOLIDBODYFLOW(CHF_FRA(vecRLLFab),
                                   CHF_CONST_FRA(rllFab), 
                                   CHF_CONST_REAL(m_magnitude), 
                                   CHF_CONST_REAL(m_angle), 
                                   CHF_CONST_REAL(m_shift));

    // Conver to velFab: vector in equiangular basis
    coordSysBlock.fabVectorTransformLatLonToEquiangular(
        a_xiLoc, vecRLLFab, a_velFab);
}


CubedSpherePolyIBC::CubedSpherePolyIBC()
{
}
/*--------------------------------------------------------------------*/
CubedSpherePolyIBC::
CubedSpherePolyIBC(const Real           a_gamma,
                   const Real           a_r0,
                   const Real           a_p0,
                   const Real&          a_radius,
                   const Vector<Real>&  a_longitude,
                   const Vector<Real>&  a_latitude,
                   const Real&          a_delta,
                   const int            a_power,
                   const Real& a_magnitude, 
                   const Real& a_angle, 
                   const Real& a_shift)                   
:m_gamma     (a_gamma),                      
 m_r0        (a_r0),                            
 m_p0        (a_p0),                            
 m_radius    (a_radius),                    
 m_longitude (a_longitude),              
 m_latitude  (a_latitude),                
 m_delta     (a_delta),                      
 m_power     (a_power),
 m_magnitude (a_magnitude),
 m_angle     (a_angle),
 m_shift     (a_shift)
{
  Real smallp;
  setFortranCommon(smallp,
                   a_gamma);

}
/*--------------------------------------------------------------------*/
CubedSpherePolyIBC::
~CubedSpherePolyIBC()
{
}
/*--------------------------------------------------------------------*/
void
CubedSpherePolyIBC::setFortranCommon(Real&      a_smallPressure,
                                     const Real a_gamma)
{
  FORT_CUBEDSPHERECOSBSETF(CHF_REAL(a_smallPressure),
                           CHF_CONST_REAL(a_gamma));

  s_isFortranCommonSet = true;
}
/*--------------------------------------------------------------------*/
PhysMappedIBC* 
CubedSpherePolyIBC::
new_physIBC()
{
  CubedSpherePolyIBC* 
    retRamp = new CubedSpherePolyIBC(m_gamma     ,
                                     m_r0        ,
                                     m_p0        ,
                                     m_radius    ,
                                     m_longitude ,
                                     m_latitude  ,
                                     m_delta     ,
                                     m_power,
                                     m_magnitude,
                                     m_angle,
                                     m_shift
                                     );

  if (m_haveTime) retRamp->setTime(m_time);
  if (m_haveCoordSys) retRamp->setCoordSys(m_coordSysPtr);

  return static_cast<PhysMappedIBC*>(retRamp);
}
/*--------------------------------------------------------------------*/
void
CubedSpherePolyIBC::
initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  LevelData<FArrayBox> advectiveVelocity(layout, SpaceDim, a_U.ghostVect());

  Real time = 0;
  if (m_haveTime) time = m_time;
  getAdvVel(advectiveVelocity, time);

  DataIterator dit = layout.dataIterator();

  const int nbells = m_longitude.size();
  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];

      // Storage for current grid
      FArrayBox& UFab = a_U[dit];

      // Box of current grid
      Box uBox = UFab.box();
      // removed by petermc, 9 Feb 2011
      // uBox &= m_domain;

      const CubedSphereShellPanelCS* coordSysBlockPtr =
        dynamic_cast<const CubedSphereShellPanelCS*>(
          m_coordSysPtr->getCoordSys(baseBox));

      CH_assert(coordSysBlockPtr);

      // For each point:
      // set RealVect Xi, which is just linear,
      // and then RealVect X( m_coordSysPtr->realCoord( Xi ) );
      // and then Real J( m_coordSysPtr->pointwiseJ( X ) );

      // Xi: mapped space coordinates
      FArrayBox XiFab(uBox, SpaceDim);

      coordSysBlockPtr->getCellMappedCoordinates(XiFab, uBox);

      // Evaluate cosine bell
      BoxIterator bit(uBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          // Get equiangular coordinates
          RealVect xi;
          xi[0] = XiFab(iv,0);
          xi[1] = XiFab(iv,1);
          RealVect lonlat;
          coordSysBlockPtr->pointTransformEquiangularToLonLat(xi, lonlat);
          Real dLon = lonlat[0];
          Real dLat = lonlat[1];

          Real dH = 0;
          for (int i = 0; i < nbells; i++)
            {
              Real dRi =
                acos(sin(m_latitude[i])*sin(dLat)
                     + cos(m_latitude[i])*cos(dLat)*cos(dLon - m_longitude[i]));
              if (dRi < m_radius)
                {
                  Real cosshape = 0.5 * (1.0 + cos(M_PI * dRi / m_radius));
                  dH = pow(cosshape, m_power);
                  dH =  m_delta * dH;
                  break;
                }
            }

          Real rho = m_r0 + dH;
          UFab(iv,URHO) = rho;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              UFab(iv, UMOMX+idir) = rho*advectiveVelocity[dit()](iv, idir);
            }
          UFab(iv, UENG) = m_p0/(m_gamma -1.0);
        }

      //now multiply by J
      FArrayBox JFab(uBox, 1);
      coordSysBlockPtr->pointwiseJ(JFab, XiFab, uBox);
      for (int comp = 0; comp < UFab.nComp(); comp++)
        {
          UFab.mult(JFab, 0, comp);
        }
    }
}
/*--------------------------------------------------------------------*/
void CubedSpherePolyIBC::primBC(FArrayBox&             a_WGdnv,
                                const FArrayBox&       a_Wextrap,
                                const FArrayBox&       a_W,
                                const FArrayBox *const a_unitNormalBasisPtr,
                                const Interval&        a_velIntv,
                                const int&             a_dir,
                                const Side::LoHiSide&  a_side,
                                const Real&            a_time)
{
  // Do nothing.
  // MayDay::Error("no supposed to get to primBC");
}
/*--------------------------------------------------------------------*/
void CubedSpherePolyIBC::setBdrySlopes(FArrayBox&       a_dW,
                                       const FArrayBox& a_W,
                                       const int&       a_dir,
                                       const Real&      a_time)
{
  CH_assert(s_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  MayDay::Error("CubedSpherePolyIBC::setBdrySlopes: is not expected for fourth-"
                "order solutions");
}
/*--------------------------------------------------------------------*/

void CubedSpherePolyIBC::artViscBC(
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
}
