#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHPARTICLEI_H_
#define _PATCHPARTICLEI_H_

#include <iomanip>

#include "NamespaceHeader.H"

// Update the particle variables, return the maximum part. speed
// on this patch/grid
template <class P, class FTypeBox>
void PatchParticle::updateState(ListBox<P>&     a_listBox,
                                Real&           a_maxPartSpeed,
                                const FTypeBox& a_field,
                                const Real&     a_dt,
                                const Box&      a_box)
{
  CH_assert(isDefined());
  a_maxPartSpeed = KDK_Pusher(a_listBox,a_field,a_box);
}

// apply second order corrections to particle positions and velocities
template <class P, class FTypeBox>
void PatchParticle::correctState(ListBox<P>&     a_listBox,
                                 const FTypeBox& a_field,
                                 const Real&     a_dt,
                                 const Box&      a_box)
{
  CH_assert(isDefined());
  CH_assert(m_isPusherCoeffSet==true);

  // correct particle's x and v
  KDK_Corrector(a_listBox, a_field, a_box);
}

/// Kick, then drift using current velocity and acceleration
template <class P, class FTypeBox>
Real PatchParticle::KDK_Pusher(ListBox<P>&     a_listBox,
                               const FTypeBox& a_field,
                               const Box&      a_box)
{
  CH_assert(isDefined());
  CH_assert(m_isPusherCoeffSet == true);

  Real maxSpeed = 1.0e-12;

  // interpolate the force to the particle positions
  m_meshInterp->interpolate(a_listBox.listItems(), 
			    a_field,
			    m_poissonInterpType);

  for (ListIterator<P> li( a_listBox.listItems() ); li; ++li)
    {
      RealVect& x = li().position();
      RealVect& v = li().velocity();
      RealVect& a = li().acceleration();

      // kick and then drift particles
      v *= m_coef_exp_pred;
      v += m_coef_acc_pred * a;
      x += m_coef_vel      * v;

      for (int dir=0; dir<SpaceDim; ++dir)
        {
          // get directional velocity and acceleration
          Real vi = Abs(v[dir]);
          Real ai = Abs(a[dir]) * m_dx;

          if (ai > 1e-3*vi*vi) {
            vi = ai / (sqrt(vi*vi + 2.0*ai) - vi);
          }
          maxSpeed = Max(vi, maxSpeed);
        }

      RealVect rightEdge(D_DECL(10.0, 10.0, 10.0));
      RealVect leftEdge(D_DECL(0.0, 0.0, 0.0));
      RealVect Lbox(D_DECL(10.0, 10.0, 10.0));

      // enforce periodic boundary conditions here
      for (int idir = 0; idir < CH_SPACEDIM; idir++)
      {
        if (x[idir] > rightEdge[idir])
        {
         x[idir] = x[idir] - Lbox[idir];
        }
      
        else if (x[idir] < leftEdge[idir])
        {
         x[idir] = x[idir] + Lbox[idir];
        }
      }
    }

  return maxSpeed;
}


// KDK corrector scheme: uses coeffiecients set in KDK_PusherCoeff(...)
template <class P, class FTypeBox>
void PatchParticle::KDK_Corrector(ListBox<P>&     a_listBox,
                                  const FTypeBox& a_field,
                                  const Box&      a_box)
{
  CH_assert(isDefined());
  CH_assert(m_isPusherCoeffSet==true);

  for (ListIterator<P> li( a_listBox.listItems() ); li; ++li)
    {
      RealVect& v = li().velocity();
      RealVect& a = li().acceleration();

      // apply final kick
      v *= m_coef_exp;
      v += m_coef_acc * a;
    }
}

// Compute the maximum wave speed
template <class P>
Real PatchParticle::getMaxPartSpeed(const ListBox<P>& a_listBox,
                                    const Box&        a_box)
{
  CH_assert(isDefined());

  Real maxVel = 1.0e-12;

  for (ListIterator<P> li( a_listBox.listItems() ); li; ++li)
    {
      for (int dir=0; dir<SpaceDim; ++dir)
        {
          Real vi= Abs(li().velocity(dir));
          maxVel = Max(vi,maxVel);
        }
    }

  return maxVel;
}

// Distribute particle mass onto the grid. If the particles' clouds
// extend outside the current box, an error will be thrown.
template <class P>
void PatchParticle::deposit(FArrayBox&        a_rho,
                            const ListBox<P>& a_listBox,
                            const Box&        a_box)
{
  CH_assert(isDefined());

  CH_assert(m_dx == a_listBox.meshSpacing()[0]);
  CH_assert(a_rho.box().contains(a_box));

  m_meshInterp->deposit(a_listBox.listItems(), 
			a_rho,
			m_poissonInterpType);
}

#include "NamespaceFooter.H"

#endif
