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
#include "GaussianIBC.H"

// Null constructor
GaussianIBC::GaussianIBC()
{
  m_params_are_set = false;
}
// Sets parameters used by initial conditions
//     a_r0            - Full width of Gaussian.
void GaussianIBC::setParams(Real a_r0, RealVect a_center, RealVect a_x0)
{
  CH_assert(m_params_are_set == false);

  m_rr0 = a_r0;
  m_x0 = a_x0;
  m_center = a_center;

  m_params_are_set = true;
}


/// set uniform velocity field
void
GaussianIBC::setUniformVel(const RealVect& a_vel)
{
  m_velType = uniform;
  m_uniformVel = a_vel;
}

// set parameter for solid-body rotation
void
GaussianIBC::setSolidBodyRotation(const RealVect& a_rotationCenter,
                                  const Real a_omega)
{
  m_velType = solidBody;
  m_rotationCenter = a_rotationCenter;
  m_omega = a_omega;
}

// set parameters for translating-oscilating velocity
void
GaussianIBC::setTranslatingOscillation(
                                       const RealVect& a_translation_vel,
                                       const Real a_oscillation_amplitude)
{
  m_velType = transOsc;
  m_uniformVel = a_translation_vel;
  m_oscAmp = a_oscillation_amplitude;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new BasicIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
BasicIBC* GaussianIBC::new_basicIBC()
{
  GaussianIBC* retval = new GaussianIBC();
  retval -> setParams(m_rr0,m_center, m_x0);
  //[NOTE: do this even though setParams() will set it true
  //       because it might not have been set in the source
  //       object (this).  Of course, that would be a bad idea
  //       because then any object created by the factor wont
  //       be usable.  Caveat usor.  -dbs]

  retval ->m_velType = m_velType;
  retval->m_uniformVel = m_uniformVel;
  retval->m_omega = m_omega;
  retval->m_oscAmp = m_oscAmp;
  retval->m_rotationCenter = m_rotationCenter;
  retval->m_params_are_set = m_params_are_set ;

  return static_cast<BasicIBC*>(retval);
}

// Set up initial conditions
void GaussianIBC::initialize(LevelData<FArrayBox>& a_phi,
                             const ProblemDomain& a_domain,
                             const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                             Real a_dx,
                             Real a_time)
{
  CH_assert(m_params_are_set == true);

  DataIterator dit = a_phi.boxLayout().dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    // includeJ = true
    pointVal(a_phi[dit], a_domain, a_coordSys,
             a_phi[dit].box(), true, a_dx, a_time);
  }

  // convert point values into 4th-order cell averages
  fourthOrderAverage(a_phi);
}


void
GaussianIBC::advVel(LevelData<FArrayBox>& a_advVel,
                    const ProblemDomain& a_domain,
                    const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                    Real a_dXi, Real a_time)
{
  Real twoPi = 8.0*atan(1.0);
  DataIterator dit = a_advVel.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
       FArrayBox X(a_advVel[dit].box(),SpaceDim);
       if (m_velType != uniform)
       {
          BoxIterator bit(X.box());
          for (bit.begin(); bit.ok(); ++bit)
          {
             IntVect iv = bit();
             RealVect Xi;
             for (int idir=0; idir<SpaceDim; idir++)
             {
                Xi[idir] = ( 0.5 + iv[idir] ) * a_dXi;
             }
             RealVect xtmp( a_coordSys.realCoord( Xi ) );
             for (int idir=0; idir<SpaceDim; idir++)
             {
                X(iv,idir) = xtmp[idir];
             }
          }
       }

       for (int dir=0; dir<SpaceDim; dir++)
        {
          if (m_velType == uniform)
            {
              a_advVel[dit].setVal(m_uniformVel[dir],dir);
            }
          else if (m_velType == solidBody)
            {
              FArrayBox& advVel = a_advVel[dit];

              BoxIterator bit(advVel.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  if (dir == 0)
                    {
                      //xVel = -y
                      Real y = X(iv,1) - m_rotationCenter[1];
                      advVel(iv,dir) = -twoPi*m_omega*y;
                    }
                  else if (dir == 1)
                    {
                      // yVel = x
                      Real x = X(iv,0) - m_rotationCenter[0];
                      advVel(iv,dir) = twoPi*m_omega*x;
                    }
                  else
                    {
                      // zvel = 0 for now
                      advVel(iv,dir) = 0.0;
                    }
                } // end loop over faces
            } // end if solid body

          else if (m_velType == transOsc)
            {
              FArrayBox& advVel = a_advVel[dit];
              Real speed = sqrt( m_uniformVel.dotProduct( m_uniformVel ) );
              Real L = speed / m_uniformVel.dotProduct( RealVect::Unit );

              BoxIterator bit(advVel.box());
              for (bit.begin();bit.ok();++bit)
                {
                  // get cell center
                  IntVect iv = bit();
                  RealVect xi;
                  for (int idir=0; idir<SpaceDim; idir++)
                    {
                      xi[idir] = ( 0.5 + iv[idir] ) * a_dXi;
                    }
                  RealVect x( a_coordSys.realCoord( xi ) );
#if 1
                  Real zeta = x.dotProduct( m_uniformVel ) / ( L * speed );
                  Real osc = m_oscAmp * cos( twoPi * zeta );
                  D_TERM( advVel(iv,0) = L*L*(m_uniformVel[0]-m_uniformVel[1]*osc);,
                          advVel(iv,1) = L*L*(m_uniformVel[1]+m_uniformVel[0]*osc);,
                          advVel(iv,2) = 0; );
#else
                  D_TERM( advVel(iv,0) = 1;,
                          advVel(iv,1) = 1+0.25*sin(twoPi*x[0])*sin(twoPi*x[1]);,
                          advVel(iv,2) = 0; );
#endif
                }

            }

          else
            {
              MayDay::Error("GaussianBC::advVel -- bad velType");
            }
        } // end loop over face directions
    }
  fourthOrderAverage(a_advVel);
}


/// fill ghost cell values at domain boundaries
void
GaussianIBC::ghostCellBC(LevelData<FArrayBox>& a_phi,
                         const ProblemDomain& a_domain,
                         const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                         Real a_dx,
                         Real a_time)
{
  IntVect ghostVect = a_phi.ghostVect();
  //const DisjointBoxLayout& grids = a_phi.getBoxes();
  const Box& domainBox = a_domain.domainBox();

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];
      if (!a_domain.contains(thisPhi.box()) )
        {
          //const Box& gridBox = grids[dit];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // this is designed to ensure that
              // corner cells get filled
              IntVect tanGrow(ghostVect);
              tanGrow[dir] = 0;

              if (!a_domain.isPeriodic(dir))
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
                          pointVal(ghostData,
                                   a_domain,
                                   a_coordSys,
                                   ghostBox,
                                   true,
                                   a_dx, a_time);

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
GaussianIBC::exactSoln(
                       LevelData<FArrayBox>& a_phi,
                       const ProblemDomain& a_domain,
                       const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                       Real a_dx,
                       Real a_time)
{
  DataIterator dit = a_phi.boxLayout().dataIterator();
  //initialize(a_phi, a_domain, a_dx, a_time);

  // do this by computing an exact finer solution and then averaging
  // down

  int nRef = 8;
  const DisjointBoxLayout& grids = a_phi.getBoxes();
  DisjointBoxLayout finerGrids;
  refine(finerGrids, grids, nRef);

  Real fineDx = a_dx/nRef;
  ProblemDomain finerDomain(a_domain);
  finerDomain.refine(nRef);
  LevelData<FArrayBox> finerPhi(finerGrids,
                                a_phi.nComp(),
                                a_phi.ghostVect());

   // compute 4th-order solution on finer mesh
  initialize(finerPhi,
             finerDomain,
             a_coordSys,
             fineDx,
             a_time);

  CoarseAverage averager(finerGrids,
                         grids,
                         a_phi.nComp(),
                         nRef);

  averager.averageToCoarse(a_phi, finerPhi);


  return;
}

/// Set up initial conditions
void
GaussianIBC::pointVal(FArrayBox& a_phi,
                      const ProblemDomain& a_domain,
                      const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                      const Box& a_box,
                      bool a_includeJ,
                      Real a_dXi,
                      Real a_time)
{
  int ncomp = a_phi.nComp();

  // Box of current grid
  //Box intersectBox(a_phi.box());

  // compute initial values

  // compute current center
  RealVect center = m_center;
  if (m_velType == uniform)
    {
      center += a_time*m_uniformVel;
    }
  else if (m_velType == solidBody)
    {
      // amount of rotation
      Real twoPi = 8.0*atan(1.0);
      Real thetaOffset = twoPi*m_omega*a_time;
      RealVect localCenter = m_center - m_rotationCenter;
      Real radCenter = sqrt(D_TERM(localCenter[0]*localCenter[0],
                                   +localCenter[1]*localCenter[1],
                                   +0));
      Real oldTheta = acos(localCenter[0]/radCenter);
      Real newTheta = oldTheta + thetaOffset;
      D_TERM(center[0] = radCenter*cos(newTheta);,
             center[1] = radCenter*sin(newTheta);,
             center[2] = center[2];)

        center += m_rotationCenter;
    }

  // 0.5 for cell-centered, 0 for node-centering
  RealVect meshOffset(0.5*RealVect::Unit);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (a_box.type(dir) == IndexType::NODE) meshOffset[dir] = 0.0;
    }

  BoxIterator bit(a_box);
  for (bit.begin();bit.ok();++bit)
    {
      IntVect iv = bit();
      Real dist = 0.;

      RealVect Xi;
      for (int idir=0; idir<SpaceDim; idir++)
      {
         Xi[idir] = ( meshOffset[idir] + iv[idir] ) * a_dXi;
      }
      RealVect X( a_coordSys.realCoord( Xi ) );

      for (int idir = 0;idir < SpaceDim; idir++)
        {
          Real loc = m_x0[idir] + X[idir] - center[idir];
          if (a_domain.isPeriodic(idir)  )
            {
              loc = min(abs(loc), abs(loc + 1.0));
              loc = min(abs(loc), abs(loc - 1.0));
            }

          dist = dist + loc*loc;
        }

      Real J = 1.0;
      if (a_includeJ)
        {
          J = a_coordSys.pointwiseJ( X );
        }
      for (int icomp = 0;icomp < ncomp;icomp++)
        {
          a_phi(iv,icomp) = J*exp(-dist/m_rr0/m_rr0)/pow(m_rr0,SpaceDim);;
        }
    }

}
