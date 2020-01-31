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
#include "GaussianAdvectMultiMappedIBC.H"

#include "NamespaceHeader.H"

// Null constructor
GaussianAdvectMultiMappedIBC::GaussianAdvectMultiMappedIBC()
{
  m_params_are_set = false;
}

// Sets parameters used by initial conditions
//     a_r0            - Full width of Gaussian.
void GaussianAdvectMultiMappedIBC::setParams(Real a_r0,
                                             RealVect a_center,
                                             RealVect a_x0,
                                             Real a_initialTime)
{
  CH_assert(m_params_are_set == false);

  m_rr0 = a_r0;
  m_x0 = a_x0;
  m_center = a_center;

  m_rr0powInv2 = 1. / (m_rr0 * m_rr0);
  m_rr0powInvD = 1.;
  for (int idir = 0; idir < SpaceDim; idir++) m_rr0powInvD /= m_rr0;

  m_initialTime = a_initialTime;

  m_params_are_set = true;
}

/// set uniform velocity field
void
GaussianAdvectMultiMappedIBC::setUniformVel(const RealVect& a_vel)
{
  m_velType = UNIFORM;
  m_uniformVel = a_vel;
}

// set parameter for solid-body rotation
void
GaussianAdvectMultiMappedIBC::setSolidBodyRotation(const RealVect& a_rotationCenter,
                                            const Real a_omega)
{
  m_velType = SOLIDBODY;
  m_rotationCenter = a_rotationCenter;
  m_omega = a_omega;
}

// set parameters for translating-oscilating velocity
void
GaussianAdvectMultiMappedIBC::setTranslatingOscillation(
                                                   const RealVect& a_translation_vel,
                                                   const Real a_oscillation_amplitude)
{
  m_velType = TRANSOSC;
  m_uniformVel = a_translation_vel;
  m_oscAmp = a_oscillation_amplitude;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new BasicIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* GaussianAdvectMultiMappedIBC::new_physIBC()
{
  GaussianAdvectMultiMappedIBC* retval = new GaussianAdvectMultiMappedIBC();
  retval->setParams(m_rr0, m_center, m_x0, m_initialTime);
  //[NOTE: do this even though setParams() will set it true
  //       because it might not have been set in the source
  //       object (this).  Of course, that would be a bad idea
  //       because then any object created by the factor wont
  //       be usable.  Caveat usor.  -dbs]

  retval->m_velType = m_velType;
  retval->m_uniformVel = m_uniformVel;
  retval->m_omega = m_omega;
  retval->m_oscAmp = m_oscAmp;
  retval->m_rotationCenter = m_rotationCenter;
  retval->m_params_are_set = m_params_are_set ;
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void GaussianAdvectMultiMappedIBC::initialize(LevelData<FArrayBox>& a_phi)
{
  CH_assert(m_isDefined);
  CH_assert(m_haveCoordSys);
  CH_assert(m_params_are_set);
  CH_assert(m_haveTime);

  /*
    Initialize the advected quantity.
  */
  const DisjointBoxLayout& layout = a_phi.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    FArrayBox& phiFab = a_phi[dit];
    const Box& bxBase = layout[dit];
    const Box& bx = phiFab.box();
    pointVal(phiFab, bx, bxBase, false); // a_includeJ);
  }

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_phi);
}


void
GaussianAdvectMultiMappedIBC::advVel(FArrayBox& a_advVelFab,
                                     Real a_time)
{
  Real twoPi = 8.0*atan(1.0);

  const Box& bx = a_advVelFab.box();

  int blockNumber = m_coordSysPtr->whichBlockOverlap(bx);
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
          a_advVelFab.setVal(m_uniformVel[dir], dir);
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
                  a_advVelFab(iv,dir) = -twoPi*m_omega*y;
                }
              else if (dir == 1)
                {
                  // yVel = x
                  Real x = XFab(iv,0) - m_rotationCenter[0];
                  a_advVelFab(iv,dir) = twoPi*m_omega*x;
                }
              else
                {
                  // zvel = 0 for now
                  a_advVelFab(iv,dir) = 0.0;
                }
            } // end loop over faces
        } // end if solid body
      else if (m_velType == TRANSOSC)
        {
          Real speed = sqrt( m_uniformVel.dotProduct( m_uniformVel ) );
          Real L = speed / m_uniformVel.dotProduct( RealVect::Unit );
          BoxIterator bit(bx);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect X;
              for (int idir = 0;idir < SpaceDim; idir++)
                {
                  X[idir] = XFab(iv, idir);
                }
              Real zeta = X.dotProduct( m_uniformVel ) / ( L * speed );
              Real osc = m_oscAmp * cos( twoPi * zeta );
              D_TERM6( a_advVelFab(iv,0) = L*L*(m_uniformVel[0]-m_uniformVel[1]*osc);,
                       a_advVelFab(iv,1) = L*L*(m_uniformVel[1]+m_uniformVel[0]*osc);,
                       a_advVelFab(iv,2) = 0; ,
                       a_advVelFab(iv,3) = 0; ,
                       a_advVelFab(iv,4) = 0; ,
                       a_advVelFab(iv,5) = 0; );
            }
        }
      else
        {
          MayDay::Error("GaussianBC::advVel -- bad velType");
        }
    } // end loop over face directions

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_advVel);
}


/// fill ghost cell values at domain boundaries
void
GaussianAdvectMultiMappedIBC::ghostCellBC(LevelData<FArrayBox>& a_phi)
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
GaussianAdvectMultiMappedIBC::exactSoln(
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
GaussianAdvectMultiMappedIBC::pointVal(FArrayBox& a_phi,
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
      center += (m_initialTime + m_time)*m_uniformVel;
    }
  else if (m_velType == SOLIDBODY)
    {
      // amount of rotation
      Real twoPi = 8.0*atan(1.0);
      Real thetaOffset = twoPi*m_omega*(m_initialTime + m_time);
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
      for (int icomp = 0;icomp < ncomp;icomp++)
        {
          a_phi(iv, icomp) = exp(-dist * m_rr0powInv2) * m_rr0powInvD;
        }
    }
  if (a_includeJ)
    {
      FArrayBox JFab(a_box, 1);
      coordSysCurrentPtr->pointwiseJ(JFab, XiFab, a_box);
      for (int icomp = 0;icomp < ncomp;icomp++)
        {
          a_phi.mult(JFab, a_box, 0, icomp);
        }
    }
}

#include "NamespaceFooter.H"
