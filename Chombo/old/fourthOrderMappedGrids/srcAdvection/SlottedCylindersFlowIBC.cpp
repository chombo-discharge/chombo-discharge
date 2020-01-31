#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SlottedCylindersFlowIBC.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

///////////////////////////////////////////////////////////////////////////////

// Null constructor
SlottedCylindersFlowIBC::SlottedCylindersFlowIBC()
{
  // CH_assert(false);
  // m_params_are_set = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
SlottedCylindersFlowIBC::SlottedCylindersFlowIBC(
                     const Real&     a_hmax,
                     const Real&     a_radius,
                     const Vector<Real>&  a_longitude,
                     const Vector<Real>&  a_latitude,
                     const Real&     a_background,
                     const Real&     a_period,
                     const Real&     a_kappa,
                     const Real&     a_evalTime) : DefFlowIBC(a_period, a_kappa, a_evalTime)
{
  m_hmax = a_hmax;
  m_radius = a_radius;
  m_longitude = Vector<Real>(a_longitude);
  m_latitude = Vector<Real>(a_latitude);
  m_background = a_background;
  // DefFlowIBC(a_period, a_kappa, a_evalTime);
}

SlottedCylindersFlowIBC::~SlottedCylindersFlowIBC()
{
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* SlottedCylindersFlowIBC::new_physIBC()
{
  SlottedCylindersFlowIBC* retval = new SlottedCylindersFlowIBC();
/*
  retval->m_isFortranCommonSet = m_isFortranCommonSet;
*/
  if (m_gotTime) retval->setTime(m_time);
  if (m_gotCoordSys) retval->setCoordSys(m_coordSysPtr);
  retval->m_hmax = m_hmax;
  retval->m_radius = m_radius;
  retval->m_longitude = Vector<Real>(m_longitude);
  retval->m_latitude = Vector<Real>(m_latitude);
  retval->m_background = m_background;
  retval->m_period = m_period;
  retval->m_kappa = m_kappa;
  retval->m_evalTime = m_evalTime;

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void SlottedCylindersFlowIBC::initializeUnified(LevelData<FArrayBox>& a_U,
                                            bool a_includeJ)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_gotCoordSys);
  CH_assert(m_gotTime);

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  int nComp = a_U.nComp();
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();

  const int ncylinders = m_longitude.size();
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

          // Pointwise height
          Real dH = m_background;
          for (int i = 0; i < ncylinders; i++)
            {
              Real dRi =
                acos(sin(m_latitude[i])*sin(dLat)
                     + cos(m_latitude[i])*cos(dLat)*cos(dLon - m_longitude[i]));
              if (dRi <= m_radius)
                { // The point is inside full cylinder i.
                  if (fabs(dLon - m_longitude[i]) >= m_radius / 6.0)
                    { // Longitude of the point is far from m_longitude[i].
                      dH = m_hmax;
                    }
                  else // fabs(dLon - m_longitude[i]) < m_radius / 6.0
                    { // Longitude of the point is close to m_longitude[i].
                      Real diffLat = dLat - m_latitude[i];
                      // signi is +1 for i even, -1 for i odd
                      int signi = (i % 2 == 0) ? 1 : -1;
                      if (signi * diffLat < - 5.0 / 12.0 * m_radius)
                        {
                          dH = m_hmax;
                        }
                    }
                  if (dH == m_background)
                    {
                      int dummy_unused; dummy_unused = 0;
                    }
                  break;
                }
            }
          UFab(iv,0) = dH;
        }

      if (a_includeJ)
        { // Multiply every component of U by J.
          // Get pointwise J (not cell-averaged <J>).
          FArrayBox JFab(uBox, 1);
          coordSysBlockPtr->pointwiseJ(JFab, XiFab, uBox);
          for (int comp = 0; comp < nComp; comp++)
            {
              UFab.mult(JFab, 0, comp);
            }
        }
    }
}

#include "NamespaceFooter.H"
