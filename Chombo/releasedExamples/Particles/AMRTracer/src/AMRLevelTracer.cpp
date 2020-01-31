#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include "parstream.H"

#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "computeSum.H"
#include "computeNorm.H"
#include "ParticleIO.H"
#include "AMRLevelTracer.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

/*******/
AMRLevelTracer::
~AMRLevelTracer()
{}

/********/
void
AMRLevelTracer::
define(AdvectionVelocityFunction a_advFunc,
       const Real&               a_cfl,
       const Real&               a_domainLength,
       const int&                a_maxBoxSize)
{
  m_isDefined = true;
  m_advFunc = a_advFunc;
  m_cfl = a_cfl;
  m_domainLength = a_domainLength;
  m_maxBoxSize = a_maxBoxSize;
}

/********/
void AMRLevelTracer::define(AMRLevel*            a_coarserLevelPtr,
			    const ProblemDomain& a_problemDomain,
			    int                  a_level,
			    int                  a_refRatio)
{
  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  if (a_coarserLevelPtr != NULL)
    {
      AMRLevelTracer* amrPartPtr = dynamic_cast<AMRLevelTracer*>(a_coarserLevelPtr);

      if (amrPartPtr != NULL)
        {
          define(amrPartPtr->m_advFunc,
                 amrPartPtr->m_cfl,
                 amrPartPtr->m_domainLength,
                 amrPartPtr->m_maxBoxSize);
        }
      else
        {
          MayDay::Error("AMRLevelTracer::define: a_coarserLevelPtr is not castable to AMRLevelTracer*");
        }
    }

  // Compute the grid spacing
  m_dx          = m_domainLength/a_problemDomain.domainBox().longside();
  m_meshSpacing = RealVect(D_DECL(m_dx, m_dx, m_dx));
  m_origin      = RealVect(D_DECL(0.0, 0.0, 0.0));

  m_numGhost = 1;
}

/*******/
Real
AMRLevelTracer::
advance()
{
  if (m_verbosity >= 2)
    {
      pout() << "AMRLevelTracer::advance " << m_level << endl;
    }

  // shift all the particles according to the velocity field.
  RealVect shift;
  RealVect velocity;

  m_maxVelocityLocal = 0.0;
  DataIterator dit(m_grids);
  for (dit.begin(); dit.ok(); ++dit)
  {
    for (ListIterator<Particle> li (m_PNew[dit].listItems()); li.ok(); ++li)
    {
      RealVect x=li().position();

      // grab velocity at Particle location
      Real velocityMag = 0.0;
      for (int dir = 0; dir < CH_SPACEDIM; ++dir)
      {
        velocity[dir] = m_advFunc(x, dir);
        velocityMag += velocity[dir]*velocity[dir];
      }

      velocityMag = sqrt(velocityMag);
      m_maxVelocityLocal = max(m_maxVelocityLocal, velocityMag);

      shift = velocity * m_dt;

      RealVect new_position = x + shift;

      // enforce periodic boundary conditions
      for (int d = 0; d < CH_SPACEDIM; d++)
      {
        if (new_position[d] < m_origin[d])
        {
          new_position[d] += m_domainLength;
        }

        if (new_position[d]  > m_domainLength + m_origin[d])
        {
          new_position[d] -= m_domainLength;
        }
      }

      li().setPosition(new_position);
    }
  }

  // now collect any invalid particles
  for (dit.reset(); dit.ok(); ++dit)
  {
    m_PNew[dit].getInvalidDestructive(m_PNew.outcast(), m_grids[dit]);
  }  

  // rebin
  m_PNew.gatherOutcast();
  m_PNew.remapOutcast();

  // Update the time and store the new timestep
  m_time += m_dt;
  return computeDt();
}

/*******/
void
AMRLevelTracer::
postTimeStep()
{
  // if this level just reached the same time as the next coarser level, transfer particles
  if (m_hasCoarser)
  {

    AMRLevelTracer* amrPartCoarserPtr = getCoarserLevel();

    const Real eps = 5.e-2;
    const Real crseTime = amrPartCoarserPtr->time();
    const bool stepsLeft = Abs(crseTime - m_time) > (eps*m_dt);
    if (!stepsLeft)
    {
      // First, collect all the particles that go from crse to fine:
      collectValidParticles(m_PNew.outcast(), 
                            amrPartCoarserPtr->m_PNew, 
                            m_PVR.mask(),
                            m_meshSpacing, 
                            amrPartCoarserPtr->m_ref_ratio);

    // now take the particles that have left the fine level and put them on the crse
    collectValidParticles(amrPartCoarserPtr->m_PNew.outcast(), m_PNew, 
                          m_PVR.mask(), m_meshSpacing, 1, true);

    // finally, remap outcasts on both levels
    m_PNew.remapOutcast();
    amrPartCoarserPtr->m_PNew.remapOutcast();

    } 
  } 

  if (m_level == 0)
  {
    int totalParticleCount = countItems();
    int totalOutcastCount = countOutcast();

    AMRLevelTracer* nextFinestLevelPtr = getFinerLevel();

    while (nextFinestLevelPtr != NULL)
    {
      totalParticleCount += nextFinestLevelPtr->countItems();
      totalOutcastCount += nextFinestLevelPtr->countOutcast();
      nextFinestLevelPtr = nextFinestLevelPtr->getFinerLevel();
    }

    pout() << totalParticleCount << endl;
    pout() << totalOutcastCount << endl;
    pout() << endl;
  }
}

/*******/
void
AMRLevelTracer::
tagCells(IntVectSet& a_tags)
{
  int numPointsLevel = 32 * pow(2, m_level);
  int loCorner = (numPointsLevel - 16) / 2;
  int hiCorner = loCorner + 16 - 1;

  IntVect lo(D_DECL(loCorner, loCorner, loCorner));
  IntVect hi(D_DECL(hiCorner, hiCorner, hiCorner));

  IntVectSet localTags;

  DataIterator dit(m_grids);
  for (dit.reset(); dit.ok(); ++dit)
  {
    const Box& b = m_grids[dit()];
    BoxIterator bit(b);
    for (bit.reset(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      bool is_contained = true;
      for (int dim = 0; dim < CH_SPACEDIM; ++dim)
      {
        is_contained = (is_contained and 
                        (lo[dim] <= iv[dim]) and
                        (hi[dim] > iv[dim]));
      }

      if (is_contained)
      {
        localTags |= iv;
      }

    }
  }

  a_tags = localTags;

}

/*******/
void
AMRLevelTracer::
tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

/// Set up data on this level before regridding
void
AMRLevelTracer::
preRegrid(int a_base_level,
          const Vector<Vector<Box> >& a_newGrids)
{

  // timer
  CH_TIME("AMRLevelTracer::preRegrid");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelTracer::preRegrid " << m_level << endl;
    }

  // get pointer to a_base_level
  AMRLevelTracer* baseLevelPtr = this;
  while (baseLevelPtr->m_level > a_base_level)
  {
    baseLevelPtr = baseLevelPtr->getCoarserLevel();
  }

  CH_assert(baseLevelPtr->m_level == a_base_level);

  DataIterator dit(m_grids);
  for (dit.reset(); dit.ok(); ++dit)
  {
    // Transfer all particles on this level to a_base_level
    for (ListIterator<Particle> li (m_PNew[dit].listItems()); li.ok(); )
    {
      baseLevelPtr->m_PNew.outcast().transfer(li);
    }

  }
    
  // This level should now be empty
  CH_assert(countItems() == 0);

  // remap outcast on a_base_level
  baseLevelPtr->m_PNew.remapOutcast();

}

/*******/
void
AMRLevelTracer::
regrid(const Vector<Box>& a_newGrids)
{
  
  // timer
  CH_TIME("AMRLevelTracer::regrid");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelTracer::regrid " << m_level << endl;
    }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);

  // Define old and new state data structures
  m_PNew.define(m_grids, m_problem_domain, m_maxBoxSize,
		m_meshSpacing, m_origin);
  m_POld.define(m_grids, m_problem_domain, m_maxBoxSize,
		m_meshSpacing, m_origin);

  // Set up data structures
  levelSetup();

}

/// postRegridding ops
void 
AMRLevelTracer::
postRegrid(int a_baseLevel)
{

  // timer
  CH_TIME("AMRLevelTracer::postRegrid");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelTracer::postRegrid " << m_level << endl;
    }

  // get pointer to a_base_level, and compute total refRatio between
  // this level and a_baseLevel.
  AMRLevelTracer* baseLevelPtr = this;
  int refRatio = 1;
  while (baseLevelPtr->m_level > a_baseLevel)
  {
    baseLevelPtr = baseLevelPtr->getCoarserLevel();
    refRatio *= baseLevelPtr->m_ref_ratio;
  }

  CH_assert(baseLevelPtr->m_level == a_baseLevel);

  // Transfer all particles from a_base_level to this level
  // if they contained within the valid region.
  collectValidParticles(m_PNew.outcast(),
                        baseLevelPtr->m_PNew, 
                        m_PVR.mask(),
                        m_meshSpacing, refRatio);

  // remap outcast on this level
  m_PNew.remapOutcast();

}

/*******/
void
AMRLevelTracer::
initialGrid(const Vector<Box>& a_newGrids)
{

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);

  // Define old and new state data structures
  m_PNew.define(m_grids, m_problem_domain, m_maxBoxSize,
		m_meshSpacing, m_origin);
  m_POld.define(m_grids, m_problem_domain, m_maxBoxSize,
		m_meshSpacing, m_origin);

  // Set up data structures
  levelSetup();
}

/*******/
void
AMRLevelTracer::
initialData()
{

  if (m_level == 0)
  {  

    // if we are on the first level, put a particle at each cell center.
    CH_XD::List<Particle> thisList;
    DataIterator dit(m_grids);
    for (dit.reset(); dit.ok(); ++dit)
    {
      const Box thisBox = m_grids.get(dit);
      BoxIterator bit(thisBox);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect position = ((RealVect)iv + 0.5)*m_dx;
        Particle particle(1.0, position);
        thisList.append(particle);
      }

      m_PNew[dit].addItemsDestructive(thisList);
    }
  }

  else
  {

    // particles have already been created, grab them from the coarser level
    AMRLevelTracer* amrPartCoarserPtr = getCoarserLevel();

    collectValidParticles(m_PNew.outcast(), 
                          amrPartCoarserPtr->m_PNew, 
                          m_PVR.mask(),
                          m_meshSpacing, 
                          amrPartCoarserPtr->refRatio());

    // put the particles in the proper bins
    m_PNew.remapOutcast();

  }
}

/*******/
void
AMRLevelTracer::
postInitialize()
{
  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelTracer::postInitialize " << m_level << endl;
    }

  if (m_level == 0)
  {
    pout() << countItems() << " " << m_level << endl;

    AMRLevelTracer* nextFinestLevelPtr = getFinerLevel();

    while (nextFinestLevelPtr != NULL)
    {
      pout() << nextFinestLevelPtr->countItems() << " " << nextFinestLevelPtr->m_level << endl;
      nextFinestLevelPtr = nextFinestLevelPtr->getFinerLevel();
    }
  }
}

#ifdef CH_USE_HDF5

/*******/
void
AMRLevelTracer::
writeCheckpointHeader(HDF5Handle& a_handle) const
{
  writePlotHeader(a_handle);
}

/*******/
void
AMRLevelTracer::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  writePlotLevel(a_handle);
}

/*******/
void
AMRLevelTracer::
readCheckpointHeader(HDF5Handle& a_handle)
{
  MayDay::Error("Restart not currently implemented.");
}

/*******/
void
AMRLevelTracer::
readCheckpointLevel(HDF5Handle& a_handle)
{
  MayDay::Error("Restart not currently implemented.");
}

/*******/
void
AMRLevelTracer::
writePlotHeader(HDF5Handle& a_handle) const
{

  // timer
  CH_TIME("AMRLevelTracer::writePlotHeader");

  if (m_verbosity >= 3)
  {
    pout() << "AMRLevelTracer::writePlotHeader" << endl;
  }

 //
  a_handle.setGroup("Charm_global");
  a_handle.setGroup("/");

  // Setup the component names
  HDF5HeaderData header;
  int numComps = 1 + 3*CH_SPACEDIM;

  Vector<string> vectNames;
  char field_name [50];
  char comp_name [50];
  char coords[3] = {'x', 'y', 'z'};

  for (int dir = 0; dir < CH_SPACEDIM; dir++)
  {
    sprintf(field_name, "particle_position_%c", coords[dir]);
    vectNames.push_back(field_name);
  }

  for (int dir = 0; dir < CH_SPACEDIM; dir++)
  {
    sprintf(field_name, "particle_velocity_%c", coords[dir]);
    vectNames.push_back(field_name);
  }

  for (int dir = 0; dir < CH_SPACEDIM; dir++)
  {
    sprintf(field_name, "particle_acceleration_%c", coords[dir]);
    vectNames.push_back(field_name);
  }

  vectNames.push_back("particle_mass");

  for (int i = 0; i < numComps; ++i)
  {
    sprintf(comp_name, "particle_component_%d", i);
    header.m_string[comp_name] = vectNames[i];
  }

  // also write out total number of particles
  int totalParticleCount = 0;
  if (m_level == 0)
  {
    totalParticleCount = countItems();

    AMRLevelTracer* nextFinestLevelPtr = getFinerLevel();

    while (nextFinestLevelPtr != NULL)
    {
      totalParticleCount += nextFinestLevelPtr->countItems();
      nextFinestLevelPtr = nextFinestLevelPtr->getFinerLevel();
    }
  }

  header.m_int["num_particles"] = totalParticleCount;

  // Write the header
  header.writeToFile(a_handle);

}

/*******/
void
AMRLevelTracer::
writePlotLevel(HDF5Handle& a_handle) const
{
  // timer
  CH_TIME("AMRLevelTracer::writePlotLevel");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelTracer::writePlotLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
           header.m_int ["is_periodic_0"] = 1;
         else
           header.m_int ["is_periodic_0"] = 0; ,

         if (m_problem_domain.isPeriodic(1))
           header.m_int ["is_periodic_1"] = 1;
         else
           header.m_int ["is_periodic_1"] = 0; ,

         if (m_problem_domain.isPeriodic(2))
           header.m_int ["is_periodic_2"] = 1;
         else
           header.m_int ["is_periodic_2"] = 0; );

  // Write the header for this level
  header.writeToFile(a_handle);

    // now write out our particles
  writeParticlesToHDF(a_handle, m_PNew, "particles");
}
#endif

/*******/
Real
AMRLevelTracer::
computeDt()
{

  // m_maxVelocityLocal stores the max velocity on each proc
  // gather and broadcast here

  int srcProc = 0;
  Vector<Real> allVelocities(numProc());
  gather(allVelocities, m_maxVelocityLocal, srcProc);

  Real maxVelocity = 0.0;
  if(procID() == srcProc)
  {
    for(int ivec = 0; ivec < numProc(); ivec++)
    {
      maxVelocity = max(maxVelocity, allVelocities[ivec]);
    }
  }

  //broadcast the right answer to all procs
  broadcast(maxVelocity, srcProc);

  return m_cfl * m_dx / maxVelocity;
}

/*******/
Real
AMRLevelTracer::
computeInitialDt()
{

  // compute the max particle velocity locally
  RealVect velocity;
  m_maxVelocityLocal = 0.0;
  DataIterator dit(m_grids);
  for (dit.begin(); dit.ok(); ++dit)
  {
    for (ListIterator<Particle> li (m_PNew[dit].listItems()); li.ok(); ++li)
    {
      RealVect x=li().position();

      // grab velocity at Particle location
      Real velocityMag = 0.0;
      for (int dir = 0; dir < CH_SPACEDIM; ++dir)
      {
        velocity[dir] = m_advFunc(x, dir);
        velocityMag += velocity[dir]*velocity[dir];
      }

      velocityMag = sqrt(velocityMag);
      m_maxVelocityLocal = max(m_maxVelocityLocal, velocityMag);
    }
  }

  return computeDt();
}

/*******/
void
AMRLevelTracer::
levelSetup()
{
  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelTracer::levelSetup " << m_level << endl;
    }

  AMRLevelTracer* amrPartCoarserPtr = getCoarserLevel();
  AMRLevelTracer* amrPartFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrPartCoarserPtr != NULL);
  m_hasFiner   = (amrPartFinerPtr   != NULL);

  int nRefCrse = -1;

  if (m_hasCoarser)
    {
      amrPartCoarserPtr = dynamic_cast<AMRLevelTracer*>(m_coarser_level_ptr);

      nRefCrse = amrPartCoarserPtr->refRatio();

      m_PVR.define(m_grids, amrPartCoarserPtr->m_PVR.mask(), nRefCrse, m_buffer);

    }

  else
    {
      m_PVR.define(m_grids, NULL, 1, 0);
    }
}

/*******/
AMRLevelTracer*
AMRLevelTracer::
getCoarserLevel() const
{
  AMRLevelTracer* amrPartCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      amrPartCoarserPtr = dynamic_cast<AMRLevelTracer*>(m_coarser_level_ptr);

      if (amrPartCoarserPtr == NULL)
        {
          MayDay::Error("AMRLevelTracer::getCoarserLevel: dynamic cast failed");
        }
    }

  return amrPartCoarserPtr;
}

/*******/
AMRLevelTracer*
AMRLevelTracer::
getFinerLevel() const
{
  AMRLevelTracer* amrPartFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
    {
      amrPartFinerPtr = dynamic_cast<AMRLevelTracer*>(m_finer_level_ptr);

      if (amrPartFinerPtr == NULL)
        {
          MayDay::Error("AMRLevelTracer::getFinerLevel: dynamic cast failed");
        }
    }

  return amrPartFinerPtr;
}

int AMRLevelTracer::countItems() const
{

  return m_PNew.numValid();
}

int AMRLevelTracer::countOutcast() const
{

  return m_PNew.numOutcast();
}

#include "NamespaceFooter.H"
