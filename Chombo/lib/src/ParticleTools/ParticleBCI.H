#ifndef _PARTICLEBCI_H_
#define _PARTICLEBCI_H_

#include "CH_Timer.H"
#include "NamespaceHeader.H"

template <class P>
void ParticleBC::enforcePeriodic(ParticleData<P>& a_particleData,
				 const RealVect& a_leftEdge,
				 const RealVect& a_rightEdge)
{

  CH_TIME("ParticleBC::enforcePeriodic");

  RealVect Lbox = a_rightEdge - a_leftEdge;
  BoxLayout BL = a_particleData.getBoxes();
  DataIterator dit(BL);

  // for each ListBox on this proc
  for (dit.reset(); dit.ok(); ++dit)
  {
    // loop over all the particles in this ListBox
    for (ListIterator<P> lit(a_particleData[dit()].listItems()); lit; ++lit)
    {
      // and apply periodic BCs...

      P& this_particle = lit();
      RealVect& this_x = this_particle.position();
      
      for (int idir = 0; idir < CH_SPACEDIM; idir++)
      {
        if (this_x[idir] > a_rightEdge[idir])
        {
         this_x[idir] = this_x[idir] - Lbox[idir];
        }
      
        else if (this_x[idir] < a_leftEdge[idir])
        {
         this_x[idir] = this_x[idir] + Lbox[idir];
        }
      }

      this_particle.setPosition(this_x);
    }
  }
}

#include "NamespaceFooter.H"

#endif // include guard
