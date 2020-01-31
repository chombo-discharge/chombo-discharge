#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "LoHiSide.H"
#include "BoxIterator.H"
#include "FourthOrderUtil.H"
#include "CoarseAverage.H"
#include "CircleAdvectMultiMappedIBC.H"

#include "NamespaceHeader.H"

// Null constructor
CircleAdvectMultiMappedIBC::CircleAdvectMultiMappedIBC()
{
  m_params_are_set = false;
}

// Sets parameters used by initial conditions
//     a_r0            - Full width of Circle.
void CircleAdvectMultiMappedIBC::setParams(Real a_r0,
                                      RealVect a_center,
                                      RealVect a_x0)
{
  CH_assert(m_params_are_set == false);

  m_rr0 = a_r0;
  m_x0 = a_x0;
  m_center = a_center;

  m_params_are_set = true;
}

/// set uniform velocity field
void
CircleAdvectMultiMappedIBC::setUniformVel(const RealVect& a_vel)
{
  m_velType = UNIFORM;
  m_uniformVel = a_vel;
}

// set parameter for solid-body rotation
void
CircleAdvectMultiMappedIBC::setSolidBodyRotation(const RealVect& a_rotationCenter,
                                            const Real a_omega)
{
  m_velType = SOLIDBODY;
  m_rotationCenter = a_rotationCenter;
  m_omega = a_omega;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new BasicIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* CircleAdvectMultiMappedIBC::new_physIBC()
{
  CircleAdvectMultiMappedIBC* retval = new CircleAdvectMultiMappedIBC();
  retval->setParams(m_rr0,m_center, m_x0);
  //[NOTE: do this even though setParams() will set it true
  //       because it might not have been set in the source
  //       object (this).  Of course, that would be a bad idea
  //       because then any object created by the factor wont
  //       be usable.  Caveat usor.  -dbs]

  retval->m_velType = m_velType;
  retval->m_uniformVel = m_uniformVel;
  retval->m_omega = m_omega;
  retval->m_rotationCenter = m_rotationCenter;
  retval->m_params_are_set = m_params_are_set ;
  if (m_gotTime) retval->setTime(m_time);
  if (m_gotCoordSys) retval->setCoordSys(m_coordSysPtr);

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void CircleAdvectMultiMappedIBC::initialize(LevelData<FArrayBox>& a_phi)
{
  initializeUnified(a_phi, false);
}

void CircleAdvectMultiMappedIBC::initializeWithJ(LevelData<FArrayBox>& a_phi)
{
  initializeUnified(a_phi, true);
}

void CircleAdvectMultiMappedIBC::initializeUnified(LevelData<FArrayBox>& a_phi,
                                                   bool a_includeJ)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_gotCoordSys);
  CH_assert(m_params_are_set == true);
  CH_assert(m_gotTime);

  /*
    Initialize the advected quantity.
  */
  const DisjointBoxLayout& layout = a_phi.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    FArrayBox& phiFab = a_phi[dit];
    const Box& bx = phiFab.box();
    const Box& bxBase = layout[dit];
    pointVal(phiFab, bx, bxBase, a_includeJ);
  }

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_phi);
}


void
CircleAdvectMultiMappedIBC::advVel(LevelData<FArrayBox>& a_advVel,
                                   Real a_time)
{
  Real twoPi = 8.0*atan(1.0);
  const DisjointBoxLayout& layout = a_advVel.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& advVelFab = a_advVel[dit];
      const Box& bxBase = layout[dit];
      const Box& bx = advVelFab.box();

      int blockNumber = m_coordSysPtr->whichBlock(bxBase);
      const NewCoordSys* coordSysCurrentPtr =
        m_coordSysPtr->getCoordSys(blockNumber);

      FArrayBox XFab(bx, SpaceDim);
      if (m_velType != UNIFORM)
        {
          FArrayBox XiFab(bx, SpaceDim);
          coordSysCurrentPtr->getCenterMappedCoordinates(XiFab, bx);
          coordSysCurrentPtr->realCoord(XFab, XiFab, bx);
        }

      for (int dir=0; dir<SpaceDim; dir++)
        {
          if (m_velType == UNIFORM)
            {
              advVelFab.setVal(m_uniformVel[dir], dir);
            }
          else if (m_velType == SOLIDBODY)
            {
              BoxIterator bit(bx);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  if (dir == 0)
                    {
                      //xVel = -y
                      Real y = XFab(iv,1) - m_rotationCenter[1];
                      advVelFab(iv,dir) = -twoPi*m_omega*y;
                    }
                  else if (dir == 1)
                    {
                      // yVel = x
                      Real x = XFab(iv,0) - m_rotationCenter[0];
                      advVelFab(iv,dir) = twoPi*m_omega*x;
                    }
                  else
                    {
                      // zvel = 0 for now
                      advVelFab(iv,dir) = 0.0;
                    }
                } // end loop over faces
            } // end if solid body
          else
            {
              MayDay::Error("CircleBC::advVel -- bad velType");
            }
        } // end loop over face directions
    }

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_advVel);
}


/// fill ghost cell values at domain boundaries
void
CircleAdvectMultiMappedIBC::ghostCellBC(LevelData<FArrayBox>& a_phi)
{
  IntVect ghostVect = a_phi.ghostVect();
  //const DisjointBoxLayout& grids = a_phi.getBoxes();
  const Box& domainBox = m_domain.domainBox();
  const DisjointBoxLayout& layout = a_phi.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];
      if (!m_domain.contains(thisPhi.box()) )
        {
          const Box& bxBase = layout[dit];
          //const Box& gridBox = grids[dit];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // this is designed to ensure that
              // corner cells get filled
              IntVect tanGrow(ghostVect);
              tanGrow[dir] = 0;

              if (!m_domain.isPeriodic(dir))
                {
                  SideIterator sit;
                  for (sit.begin(); sit.ok(); ++sit)
                    {
                      Box ghostBox;
                      if (sit() == Side::Lo)
                        {
                          ghostBox = adjCellLo(domainBox,
                                               dir,
                                               ghostVect[dir]);
                        }
                      else
                        {
                          ghostBox = adjCellHi(domainBox,
                                               dir,
                                               ghostVect[dir]);
                        }

                      // ensure that corner cells get filled
                      ghostBox.grow(tanGrow);
                      ghostBox &= thisPhi.box();

                      if (!ghostBox.isEmpty())
                        {
                          // need to grow the ghost box by
                          // one for 4th-order averaging
                          ghostBox.grow(1);
                          FArrayBox ghostData(ghostBox,
                                              thisPhi.nComp());
                          // includeJ = true
                          pointVal(ghostData, ghostBox, bxBase, true);

                          fourthOrderAverageCell(ghostData);

                          // now copy into ghost cells
                          ghostBox.grow(-1);

                          thisPhi.copy(ghostData, ghostBox);
                        }  // end if there are domain ghost cells here
                    } // end loop over hi-lo
                } // end if not periodic in this direction
            } // end loop over directions
        } // end if phi sticks out of domain
    } // end loop over grids
}


/// compute exact solution
void
CircleAdvectMultiMappedIBC::exactSoln(
                                 LevelData<FArrayBox>& a_phi)
{
  //initialize(a_phi, m_domain);

  // do this by computing an exact finer solution and then averaging
  // down

  int nRef = 8;
  const DisjointBoxLayout& grids = a_phi.getBoxes();
  DisjointBoxLayout finerGrids;
  refine(finerGrids, grids, nRef);

  // Real fineDx = m_dx/nRef;
  ProblemDomain finerDomain(m_domain);
  finerDomain.refine(nRef);
  LevelData<FArrayBox> finerPhi(finerGrids,
                                a_phi.nComp(),
                                a_phi.ghostVect());

   // compute 4th-order solution on finer mesh:  need fineDx !
  initialize(finerPhi);

  CoarseAverage averager(finerGrids,
                         grids,
                         a_phi.nComp(),
                         nRef);

  averager.averageToCoarse(a_phi, finerPhi);

  return;
}

/// Set up initial conditions
void
CircleAdvectMultiMappedIBC::pointVal(FArrayBox& a_phi,
                                     const Box& a_box,
                                     const Box& a_bxBase,
                                     bool a_includeJ)
{
  int ncomp = a_phi.nComp();

  // Box of current grid
  //Box intersectBox(a_phi.box());
  int blockNumber = m_coordSysPtr->whichBlock(a_bxBase);
  const NewCoordSys* coordSysCurrentPtr =
    m_coordSysPtr->getCoordSys(blockNumber);

  // compute initial values

  // compute current center
  RealVect center = m_center;
  if (m_velType == UNIFORM)
    {
      center += m_time*m_uniformVel;
    }
  else if (m_velType == SOLIDBODY)
    {
      // amount of rotation
      Real twoPi = 8.0*atan(1.0);
      Real thetaOffset = twoPi*m_omega*m_time;
      RealVect localCenter = m_center - m_rotationCenter;
      Real radCenter = sqrt(D_TERM6(localCenter[0]*localCenter[0],
                                    +localCenter[1]*localCenter[1],
                                    +0, +0, +0, +0));
      Real oldTheta = acos(localCenter[0]/radCenter);
      Real newTheta = oldTheta + thetaOffset;
      D_TERM6(center[0] = radCenter*cos(newTheta);,
              center[1] = radCenter*sin(newTheta);,
              center[2] = center[2];,
              center[3] = center[3];,
              center[4] = center[4];,
              center[5] = center[5];)
        center += m_rotationCenter;
    }

  // 0.5 for cell-centered, 0 for node-centering
  RealVect meshOffset(0.5*RealVect::Unit);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (a_box.type(dir) == IndexType::NODE) meshOffset[dir] = 0.0;
    }

  FArrayBox XiFab(a_box, SpaceDim);
  FArrayBox XFab(a_box, SpaceDim);
  coordSysCurrentPtr->getCenterMappedCoordinates(XiFab, a_box);
  coordSysCurrentPtr->realCoord(XFab, XiFab, a_box);

  BoxIterator bit(a_box);
  for (bit.begin();bit.ok();++bit)
    {
      IntVect iv = bit();
      RealVect X;
      for (int idir = 0;idir < SpaceDim; idir++)
        {
          X[idir] = XFab(iv, idir);
        }
      Real dist = 0.;
      for (int idir = 0;idir < SpaceDim; idir++)
        {
          Real loc = m_x0[idir] + X[idir] - center[idir];
          if (m_domain.isPeriodic(idir)  )
            {
              loc = min(abs(loc), abs(loc + 1.0));
              loc = min(abs(loc), abs(loc - 1.0));
            }
          dist = dist + loc*loc;
        }

      Real J = 1.0;
      if (a_includeJ)
        {
          J = coordSysCurrentPtr->pointwiseJ( X );
        }
      for (int icomp = 0;icomp < ncomp;icomp++)
        {
          if (dist <= m_rr0)
            {
              a_phi(iv,icomp) = J;
            }
          else
            {
              a_phi(iv, icomp) = 0.0;
            }
        }
    }
}


// Set boundary fluxes
void CircleAdvectMultiMappedIBC::primBC(FArrayBox&            a_WGdnv,
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
void CircleAdvectMultiMappedIBC::setBdrySlopes(FArrayBox&       a_dW,
                                          const FArrayBox& a_W,
                                          const int&       a_dir,
                                          const Real&      a_time)
{
  MayDay::Error("CircleAdvectMultiMappedIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void CircleAdvectMultiMappedIBC::artViscBC(FArrayBox&       a_F,
                                      const FArrayBox& a_U,
                                      const FArrayBox& a_divVel,
                                      const int&       a_dir,
                                      const Real&      a_time)
{
}

#include "NamespaceFooter.H"
