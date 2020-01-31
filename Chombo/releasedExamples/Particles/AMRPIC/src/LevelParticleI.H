#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELPARTICLEI_H_
#define _LEVELPARTICLEI_H_

#include "MPI_util.H"
#include "ParticleData.H"

#include "NamespaceHeader.H"

// Advance the particle solution by "a_dt": method defined in PatchParticle
template <class P>
Real LevelParticle::step(ParticleData<P>&                a_levP,
                         const BoxLayoutData<FArrayBox>& a_force,
                         const Real&                     a_time,
                         const Real&                     a_dt,
                         const bool                      a_globalDtNew)
{
  // timer
  CH_TIME("LevelParticle::step");

  // Make sure everything is defined
  CH_assert(m_defined);

  // Used by the particle pusher
  m_patchParticle->setPusherCoeff(a_time, a_dt);

  // Use to restrict maximum wave speed away from zero
  Real maxPartSpeed = 1.0e-12;

  // Beginning of loop through patches/grids.
  const BoxLayout& grids = a_levP.getBoxes();

  for (DataIterator di = grids.dataIterator(); di.ok(); ++di)
    {
      // The current box
      const Box& curBox = grids.get(di);

      // The current grid of Particle variables
      ListBox<P>& curP = a_levP[di];

      // acceleration data
      const FArrayBox& force = a_force[di];

      Real maxPartSpeedGrid;

      // Update the current particle variables (velocity, position)
      // and return the maximum speed for this grid.
      m_patchParticle->updateState(curP,
                                   maxPartSpeedGrid,
                                   force,
                                   a_dt,
                                   curBox);

      // Clamp away from zero
      maxPartSpeed = Max(maxPartSpeed, maxPartSpeedGrid);
    }

  // Find the minimum of dt's over this level
  Real dtNew = m_dx / maxPartSpeed;

#ifdef CH_MPI
  if (a_globalDtNew)
    {
      Real globalMin = 1.0e300; // an absurdly large value
      int srcProc = 0;
      Vector<Real> allDts(numProc());
      gather(allDts, dtNew, srcProc);
      
      if( procID() == srcProc)
	{
	  for(int ivec = 0; ivec < numProc(); ivec++)
	    {
	      globalMin = Min(globalMin, allDts[ivec]);
	    }	  
	}
      
      broadcast(globalMin, srcProc);
      dtNew = globalMin;
    }
#endif
  
  // Return the maximum stable time step
  return dtNew;
}

// Implement 2nd order correction: method defined in PatchParticle
template <class P>
void LevelParticle::correct(ParticleData<P>&                a_levP,
                            const BoxLayoutData<FArrayBox>& a_force,
                            const Real&                     a_time,
                            const Real&                     a_dt)
{
  // timer
  CH_TIME("LevelParticle::correct");

  // Make sure everything is defined
  CH_assert(m_defined);

  // Beginning of loop through patches/grids.
  const BoxLayout& grids = a_levP.getBoxes();

  for (DataIterator di(grids); di.ok(); ++di)
    {
      // The current box
      const Box& curBox = grids.get(di);

      // The current grid of Particle variables
      ListBox<P>& curP = a_levP[di];

      // acceleration data
      const FArrayBox& force = a_force[di];

      // Correct the current particle variables and return the maximum
      // speed for this grid.
      m_patchParticle->correctState(curP,
                                    force,
                                    a_dt,
                                    curBox);
    }

  // Used by the particle pusher
  m_patchParticle->resetPusherCoeff();
}

//
template<class P>
void LevelParticle::postCorrect(ParticleData<P>& a_levP)
{
  // timer
  CH_TIME("LevelParticle::postCorrect");

  // Used by the particle pusher
  m_patchParticle->resetPusherCoeff();
}


// Find the maximum wave speed on the current level
template <class P>
Real LevelParticle::getMaxPartSpeed(const ParticleData<P>& a_levP,
                                    const bool             a_globalMax)
{
  // Initial maximum part. speed
  Real maxSpeed = 1.e-12;

  // Loop over all grids to get the maximum wave speed
  const BoxLayout& grids = a_levP.getBoxes();

  for (DataIterator di(grids); di.ok(); ++di)
    {
      const Box& curBox = grids.get(di);
      const ListBox<P>& curP = a_levP[di];
      Real speedOverBox = m_patchParticle->getMaxPartSpeed(curP,curBox);

      // Compute a running maximum
      maxSpeed = Max(maxSpeed,speedOverBox);
    }

#ifdef CH_MPI
  if (a_globalMax)
    {
      Real globalMax = -1.0e300; // an absurdly large negative value
      int srcProc = 0;
      Vector<Real> allSpeeds(numProc());
      gather(allSpeeds, maxSpeed, srcProc);
      
      if( procID() == srcProc)
	{
	  for(int ivec = 0; ivec < numProc(); ivec++)
	    {
	      globalMax = Max(globalMax, allSpeeds[ivec]);
	    }	  
	}
      
      broadcast(globalMax, srcProc);
      maxSpeed = globalMax;
    }
#endif

  // Return the maximum wave speed
  return maxSpeed;
}

#include "NamespaceFooter.H"

#endif
