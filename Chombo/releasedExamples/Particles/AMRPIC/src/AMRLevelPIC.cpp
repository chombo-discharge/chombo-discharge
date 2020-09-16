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
#include "AMRLevelPIC.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

// set the solver parameters here
#define SOLVER_RESET_PHI false
#define SOLVER_NORM_TYPE (0)
#define SOLVER_MIN_ITER  (2)
#define SOLVER_MAX_ITER  (15)
#if CHF_USE_FLOAT
#  define SOLVER_TOLERANCE (1.0e-5)
#else
#  define SOLVER_TOLERANCE (1.0e-10)
#endif

#define SOLVER_NUM_SMOOTH 4
#define SOLVER_NUM_MG     1
#define SOLVER_HANG       1.e-12
#define SOLVER_NORM_THRES 1.e-12

BiCGStabSolver<LevelData<FArrayBox> >* AMRLevelPIC::s_bottomSolver  = NULL;
RefCountedPtr<AMRPoissonOpFactory> AMRLevelPIC::s_gravOpFact(new AMRPoissonOpFactory());
RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > AMRLevelPIC::s_gravMGSolver(new AMRMultiGrid<LevelData<FArrayBox> >());

// This is a dummy function for periodic BCs
void NoOpBc(FArrayBox& a_state,
            const Box& a_box,
            const ProblemDomain& a_domain,
            Real a_dx,
            bool a_homogeneous)
{}

// Return the time interpolation parameter
Real getAlpha(const Real a_tnew, const Real a_time, const Real a_dt)
{

  const Real eps = (a_dt > 0.0) ? 1.e-2*a_dt : 1.e-10;

  const Real a_told = a_tnew - a_dt;

  Real alpha= - 1.0;
  // case 1: we are at the new time
  if (Abs(a_tnew - a_time) <= eps)
    {
      alpha = 1.0;
    }
  // case 2: we are at the old time
  else if (Abs(a_time - a_told) <= eps)
    {
      alpha = 0.0;
    }
  // case 3: we are at neither; compute alpha
  else
    {
      CH_assert( a_dt > 0.0);
      alpha= (a_time - a_told) / a_dt;
    }
  CH_assert(alpha >= 0.0);

  return alpha;
}

/*******/
AMRLevelPIC::
~AMRLevelPIC()
{

  if (s_bottomSolver != NULL)
    {
      delete s_bottomSolver;
      s_bottomSolver = NULL;
    }
  
  if (m_patchParticle != NULL)
    {
      delete m_patchParticle;
      m_patchParticle = NULL;
    }
  
  if (m_ibc != NULL)
    {
      delete m_ibc;
      m_ibc = NULL;
    }

  s_gravOpFact=RefCountedPtr<AMRPoissonOpFactory>();
  s_gravMGSolver=RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
}

/********/
void
AMRLevelPIC::
define(const ZeldovichIBC*       a_ibc,
       const Real&               a_cfl,
       const Real&               a_domainLength,
       const int&                a_maxBoxSize,
       const InterpType&         a_poissonInterpType,
       const Real&               a_refineThresh)
{
  m_isDefined = true;
  m_cfl = a_cfl;
  m_domainLength = a_domainLength;
  m_maxBoxSize = a_maxBoxSize;
  m_poissonInterpType = a_poissonInterpType;
  m_refineThresh = a_refineThresh;
  m_numForceGhost = 1; 

  // set up bottom solver
  if (s_bottomSolver == NULL)
    {
      s_bottomSolver = new BiCGStabSolver<LevelData<FArrayBox> >();
      s_bottomSolver->m_verbosity = 3;
      s_bottomSolver->m_eps = 1e-14;
    }

  if (m_patchParticle != NULL)
    {
      delete m_patchParticle;
      m_patchParticle = NULL;
    }
  
  m_patchParticle = new PatchParticle();
  
  if (m_ibc != NULL)
    {
      delete m_ibc;
      m_ibc = NULL;
    }

  m_ibc = a_ibc->new_zeldovichIBC();

  Real Hubble_0 = 0.5;
  Real Omega_m  = 1.0;
  m_cosmology.define(Hubble_0);
  m_cosmology.m_Omega_m = Omega_m;
  m_Cexp = 1.0e-2;
}

/********/
void AMRLevelPIC::define(AMRLevel*            a_coarserLevelPtr,
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
      AMRLevelPIC* amrPartPtr = dynamic_cast<AMRLevelPIC*>(a_coarserLevelPtr);

      if (amrPartPtr != NULL)
        {
          define(amrPartPtr->m_ibc,
                 amrPartPtr->m_cfl,
                 amrPartPtr->m_domainLength,
                 amrPartPtr->m_maxBoxSize,
                 amrPartPtr->m_poissonInterpType,
                 amrPartPtr->m_refineThresh);
        }
      else
        {
          MayDay::Error("AMRLevelPIC::define: a_coarserLevelPtr is not castable to AMRLevelPIC*");
        }
    }

  // Compute the grid spacing
  m_dx          = m_domainLength/a_problemDomain.domainBox().longside();
  m_meshSpacing = RealVect(D_DECL(m_dx, m_dx, m_dx));
  m_origin      = RealVect(D_DECL(0.0, 0.0, 0.0));

  m_patchParticle->define(m_problem_domain, m_dx, 
                          m_poissonInterpType, m_cosmology);
}

/*******/
Real AMRLevelPIC::advance()
{
  if (m_verbosity >= 2)
    {
      pout() << "AMRLevelPIC::advance " << m_level << endl;
    }

  // Absurdly large value, will be over-written later
  Real newDt = 1.0e9;

  // push particles using the current force field
  newDt = m_levelParticle.step(m_PNew, m_fieldNew, m_time, m_dt);

  // also advance Joint particles from finer levels
  if (m_hasFiner && !m_isThisFinestLev)
  {
    m_levelParticle.step(m_jointParticle, m_fieldNew, m_time, m_dt);
  }

  m_PNew.gatherOutcast();
  m_PNew.remapOutcast();

  m_jointParticle.gatherOutcast();
  m_jointParticle.remapOutcast();

  // Update the time
  m_time += m_dt;

  // save new Phi/Force -> old Phi/Force
  for (DataIterator di=m_phiNew.dataIterator(); di.ok(); ++di)
    {
      m_phiOld[di].copy(m_phiNew[di]);
      m_fieldOld[di].copy(m_fieldNew[di]);
    }

  // do gravity level solve unless we are on the finest level
  if (m_hasFiner && !m_isThisFinestLev)
  {
    ellipticSolver(m_level, true);
    computeForce(m_fieldNew, m_phiNew, m_time);
  }

  m_dtOld = m_dt;
  m_dtNew = m_cfl * newDt;
  return m_dtNew;
}

/*******/
void
AMRLevelPIC::
postTimeStep()
{
  // if this level just reached the same time as the next coarser level, transfer particles
  if (m_hasCoarser)
  {

    AMRLevelPIC* amrPartCoarserPtr = getCoarserLevel();

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

  // We also must now update Joint Particles, if there is a finer level.
  if (m_hasFiner)
    {
      AMRLevelPIC* amrFinePtr = getFinerLevel();
      amrFinePtr->aggregateParticles(m_jointParticle);
    } 

  // call amrSync for all levels that are now aligned in time
  const Real eps = 5.e-2;
  Real crseTime = -1.0;
  if (m_level > 0)
  {
    crseTime = m_coarser_level_ptr->time();
  }
  const bool stepsLeft = Abs(crseTime - m_time) > (eps*m_dt);

  if (m_level==0 || stepsLeft)
  {
    amrSync(m_level);
  }
}

/*******/
void
AMRLevelPIC::
tagCells(IntVectSet& a_tags)
{

  // We tag cells in which the mass exceeds some threshold value.
  m_rhs.exchange();
  IntVectSet localTags;
  const Real cellVolume = pow(m_dx, CH_SPACEDIM);
  DataIterator dit(m_grids);
  for (dit.reset(); dit.ok(); ++dit)
  {
    const Box& b = m_grids[dit()];
    BoxIterator bit(b);
    for (bit.reset(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      
      Real m = m_rhs[dit](iv);
      const Real H0 = m_cosmology.m_H0;

      // convert to a mass density
      m *= (2.0/3.0)*m_cosmology.a(m_time)/H0/H0;
      m += 1.0;
      m *= cellVolume;

      if (m > m_refineThresh)
      {
        localTags |= iv;
      }

    }
  }

  a_tags = localTags;

}

/*******/
void
AMRLevelPIC::
tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

/// Set up data on this level before regridding
void
AMRLevelPIC::
preRegrid(int a_base_level,
          const Vector<Vector<Box> >& a_newGrids)
{

  CH_TIME("AMRLevelPIC::preRegrid");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::preRegrid " << m_level << endl;
    }

  // get pointer to a_base_level
  AMRLevelPIC* baseLevelPtr = this;
  while (baseLevelPtr->m_level > a_base_level)
  {
    baseLevelPtr = baseLevelPtr->getCoarserLevel();
  }
  CH_assert(baseLevelPtr->m_level == a_base_level);

  // We are regridding on level a_base_level; all particles that
  // live on finer levels should be removed and added here. There is 
  // probably a better way to do this.
  AMRLevelPIC* amrFinerPtr = getFinerLevel();
  while (amrFinerPtr != NULL)
  {
    DataIterator dit(amrFinerPtr->m_grids);
    for (dit.reset(); dit.ok(); ++dit)
    {
      // Transfer all particles on this level to a_base_level
      for (ListIterator<Particle> li (amrFinerPtr->m_PNew[dit].listItems()); li.ok(); )
      {
        baseLevelPtr->m_PNew.outcast().transfer(li);
      }
    }

    amrFinerPtr = amrFinerPtr->getFinerLevel();

  }

  // remap outcast on a_base_level
  baseLevelPtr->m_PNew.remapOutcast();
}

/*******/
void
AMRLevelPIC::
regrid(const Vector<Box>& a_newGrids)
{
  
  // timer
  CH_TIME("AMRLevelPIC::regrid");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::regrid " << m_level << endl;
    }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);

  // Define old and new particle data structures
  m_PNew.define(m_grids, m_problem_domain, m_maxBoxSize, 
                m_meshSpacing, m_origin);
  m_POld.define(m_grids, m_problem_domain, m_maxBoxSize, 
                m_meshSpacing, m_origin);
  m_jointParticle.define(m_grids, m_problem_domain, m_maxBoxSize, 
                         m_meshSpacing, m_origin);

  // Define old and new mesh data structures
  m_rhs.define(m_grids, 1, m_numForceGhost*IntVect::Unit);

  m_phiNew.define(m_grids, 1, IntVect::Unit);
  m_phiOld.define(m_grids, 1, IntVect::Unit);
  m_phiInt.define(m_grids, 1, IntVect::Unit);

  m_fieldNew.define(m_grids, CH_SPACEDIM, m_numForceGhost*IntVect::Unit);
  m_fieldOld.define(m_grids, CH_SPACEDIM, m_numForceGhost*IntVect::Unit);

  // standard ghost exchange copier
  m_forwardCopier.define(m_grids, 
                         m_grids, 
                         m_problem_domain, 
                         m_numForceGhost*IntVect::Unit, 
                         true);

  // a reversed version of the above
  m_reverseCopier.define(m_grids, 
                         m_grids, 
                         m_problem_domain, 
                         m_numForceGhost*IntVect::Unit, 
                         true);

  m_reverseCopier.reverse();

  // initialize everything to zero (do I really want to do this?)
  DataIterator dit = m_grids.dataIterator();
  for(dit.reset(); dit.ok(); ++dit)
    {
      m_rhs[dit()].setVal(0.0);
      m_phiNew[dit()].setVal(0.0);
      m_phiOld[dit()].setVal(0.0);
      m_fieldNew[dit()].setVal(0.0);
      m_fieldOld[dit()].setVal(0.0);
    }

  // Set up data structures
  levelSetup();

  if (m_hasCoarser)
  {
    AMRLevelPIC* coarserLevelPtr = getCoarserLevel();
    DisjointBoxLayout crseGrids;
    coarsen(crseGrids,m_grids,coarserLevelPtr->refRatio());
    m_rhsCrseNew.define(crseGrids,1);
    m_rhsCrseOld.define(crseGrids,1);
  }
}

/// postRegridding ops
void 
AMRLevelPIC::
postRegrid(int a_baseLevel)
{

  // timer
  CH_TIME("AMRLevelPIC::postRegrid");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::postRegrid " << m_level << endl;
    }

  if (m_hasFiner)
  {
    m_isThisFinestLev = (m_grids.size()>0 && getFinerLevel()->m_grids.size()<=0);
  }

  if (m_isThisFinestLev)
  {
    AMRLevelPIC* coarserLevelPtr = getCoarserLevel();
    while (coarserLevelPtr != NULL)
    {
      coarserLevelPtr->m_isThisFinestLev = 0;
      coarserLevelPtr = coarserLevelPtr->getCoarserLevel();
    }

    AMRLevelPIC* finerLevelPtr = getFinerLevel();
    while (finerLevelPtr != NULL)
    {
      finerLevelPtr->m_isThisFinestLev = 0;
      finerLevelPtr = finerLevelPtr->getFinerLevel();
    }
  }

  // get pointer to a_base_level, and compute total refRatio between
  // this level and a_baseLevel.
  AMRLevelPIC* baseLevelPtr = this;
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

  // Here is where we should collect the Joint Particles from the finer level
  if (m_hasFiner)
  {
    AMRLevelPIC* amrFinePtr = getFinerLevel();
    amrFinePtr->aggregateParticles(m_jointParticle);
    m_jointParticle.remapOutcast();
  }

  // If we are on the finest level, define elliptic solvers
  if (m_isThisFinestLev)
    {
      defineEllipticSolvers(m_level);
    }

  // finally, do a full multilevel Poisson solve from the coarsest 
  // synchronized level. This involves walking through the level
  // hierarchy.
  int crsestSynchLev = m_level;
  AMRLevelPIC* amrSynchLevPtr = this;

  if (m_level > 0)
  {
    // look for the right level
    const Real eps = 5.e-2;
    const Real smallT = eps * m_dt;
    AMRLevelPIC* coarserLevelPtr = amrSynchLevPtr->getCoarserLevel();

    while (crsestSynchLev > 0)
    {
      const Real dTLevs = amrSynchLevPtr->time() - coarserLevelPtr->time();

      if (Abs(dTLevs) <= smallT)
      {
        --crsestSynchLev;
        amrSynchLevPtr= coarserLevelPtr;
        coarserLevelPtr = amrSynchLevPtr->getCoarserLevel();
      }
      else
      {
        break;
      }
    }
  }

  // need to call amrSync() from the synch lev because of
  // reference to base-level-data there
  const bool srceCorr= false;
  amrSynchLevPtr->amrSync(crsestSynchLev, srceCorr);
}

/*******/
void
AMRLevelPIC::
initialGrid(const Vector<Box>& a_newGrids)
{

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);

  // Define old and new particle data structures
  m_PNew.define(m_grids, m_problem_domain, m_maxBoxSize, 
                m_meshSpacing, m_origin);
  m_POld.define(m_grids, m_problem_domain, m_maxBoxSize, 
                m_meshSpacing, m_origin);
  m_jointParticle.define(m_grids, m_problem_domain, m_maxBoxSize, 
                         m_meshSpacing, m_origin);

  // Define old and new mesh data structures
  m_rhs.define(m_grids, 1, m_numForceGhost*IntVect::Unit);

  m_phiNew.define(m_grids, 1, IntVect::Unit);
  m_phiOld.define(m_grids, 1, IntVect::Unit);
  m_phiInt.define(m_grids, 1, IntVect::Unit);

  m_fieldNew.define(m_grids, CH_SPACEDIM, m_numForceGhost*IntVect::Unit);
  m_fieldOld.define(m_grids, CH_SPACEDIM, m_numForceGhost*IntVect::Unit);

  // standard ghost exchange copier
  m_forwardCopier.define(m_grids, 
                         m_grids, 
                         m_problem_domain, 
                         m_numForceGhost*IntVect::Unit, 
                         true);

  // a reversed version of the above
  m_reverseCopier.define(m_grids, 
                         m_grids, 
                         m_problem_domain, 
                         m_numForceGhost*IntVect::Unit, 
                         true);
  m_reverseCopier.reverse();

  // initialize everything to zero
  DataIterator dit = m_grids.dataIterator();
  for(dit.reset(); dit.ok(); ++dit)
    {
      m_rhs[dit()].setVal(0.0);
      m_phiNew[dit()].setVal(0.0);
      m_phiOld[dit()].setVal(0.0);
      m_fieldNew[dit()].setVal(0.0);
      m_fieldOld[dit()].setVal(0.0);
    }

  // Set up data structures
  levelSetup();

  if (m_hasCoarser)
  {
    AMRLevelPIC* coarserLevelPtr = getCoarserLevel();
    DisjointBoxLayout crseGrids;
    coarsen(crseGrids,m_grids,coarserLevelPtr->refRatio());
    m_rhsCrseNew.define(crseGrids,1);
    m_rhsCrseOld.define(crseGrids,1);
  }
}

/*******/
void
AMRLevelPIC::
initialData()
{

  // we generate particles on the coarsest level only
  if (m_level == 0)
    {
      m_ibc->initialize(m_PNew);

      // put the particles in the proper bins
      m_PNew.remapOutcast();
      CH_assert(m_PNew.isClosed());
    }

  else
  {

    // particles have already been created, grab them from the coarser level
    AMRLevelPIC* amrPartCoarserPtr = getCoarserLevel();

    collectValidParticles(m_PNew.outcast(), 
                          amrPartCoarserPtr->m_PNew, 
                          m_PVR.mask(),
                          m_meshSpacing, 
                          amrPartCoarserPtr->refRatio());

    // put the particles in the proper bins
    m_PNew.remapOutcast();
    CH_assert(m_PNew.isClosed());
  }
}

/*******/
void
AMRLevelPIC::
postInitialize()
{
  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::postInitialize " << m_level << endl;
    }

  if (m_hasFiner)
  {
    m_isThisFinestLev = (m_grids.size()>0 && getFinerLevel()->m_grids.size()<=0);
  }

  if (m_isThisFinestLev)
  {
    AMRLevelPIC* coarserLevelPtr = getCoarserLevel();
    while (coarserLevelPtr != NULL)
    {
      coarserLevelPtr->m_isThisFinestLev = 0;
      coarserLevelPtr = coarserLevelPtr->getCoarserLevel();
    }
  }

  // As in postRegrid, we also must do Joint particles here.
  if (m_hasFiner)
  {
    AMRLevelPIC* amrFinePtr = getFinerLevel();
    amrFinePtr->aggregateParticles(m_jointParticle);
    m_jointParticle.remapOutcast();
  } 

  // If we are on the finest level, define elliptic solvers
  if (m_isThisFinestLev)
    {
      defineEllipticSolvers(m_level);
    }

  if (m_level == 0)
    {
      const bool srceCorr = false;
      amrSync(m_level, srceCorr);
    }
}

#ifdef CH_USE_HDF5

/*******/
void
AMRLevelPIC::
writeCheckpointHeader(HDF5Handle& a_handle) const
{
  writePlotHeader(a_handle);
}

/*******/
void
AMRLevelPIC::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  writePlotLevel(a_handle);
}

/*******/
void
AMRLevelPIC::
readCheckpointHeader(HDF5Handle& a_handle)
{
  MayDay::Error("Restart not currently implemented.");
}

/*******/
void
AMRLevelPIC::
readCheckpointLevel(HDF5Handle& a_handle)
{
  MayDay::Error("Restart not currently implemented.");
}

/*******/
void
AMRLevelPIC::
writePlotHeader(HDF5Handle& a_handle) const
{

  // timer
  CH_TIME("AMRLevelPIC::writePlotHeader");

  if (m_verbosity >= 3)
  {
    pout() << "AMRLevelPIC::writePlotHeader" << endl;
  }

  a_handle.setGroup("/");

  // Setup the component names
  HDF5HeaderData header;

  Vector<string> vectNames;
  char field_name [50];
  char comp_name [50];
  char coords[3] = {'x', 'y', 'z'};

  // setup mesh field names
  int numMeshComps = 2 + CH_SPACEDIM;
  vectNames.push_back("rhs");
  vectNames.push_back("phi");
  for (int dir = 0; dir < CH_SPACEDIM; dir++)
  {
    sprintf(field_name, "gravitational_field_%c", coords[dir]);
    vectNames.push_back(field_name);
  }

  for (int i = 0; i < numMeshComps; ++i)
  {
    sprintf(comp_name, "component_%d", i);
    header.m_string[comp_name] = vectNames[i];
  }

  header.m_int["num_components"] = numMeshComps;

  // now do the same for the particles
  vectNames.clear();
  int numComps = 1 + 3*CH_SPACEDIM;

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

    AMRLevelPIC* nextFinestLevelPtr = getFinerLevel();

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
AMRLevelPIC::
writePlotLevel(HDF5Handle& a_handle) const
{
  // timer
  CH_TIME("AMRLevelPIC::writePlotLevel");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::writePlotLevel" << endl;
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

  // write out mesh fields
  LevelData<FArrayBox> outputData;
  int numComps = 2 + CH_SPACEDIM;
  outputData.define(m_grids, numComps, IntVect::Zero);

  // do copies
  m_rhs.copyTo(Interval(0, 0), outputData, Interval(0, 0));
  m_phiNew.copyTo(Interval(0, 0), outputData, Interval(1, 1));
  m_fieldNew.copyTo(Interval(0, CH_SPACEDIM - 1), outputData, Interval(2, CH_SPACEDIM + 1));

  write(a_handle, m_rhs.boxLayout());
  write(a_handle, outputData, "data");

}
#endif

/*******/
Real
AMRLevelPIC::
computeDt()
{
  const Real a = m_cosmology.a(m_time);
  const Real adot = m_cosmology.adot(m_time);
  return m_Cexp * m_dx * a / adot;
}

/*******/
Real
AMRLevelPIC::
computeInitialDt()
{
  const Real a = m_cosmology.a(m_time);
  const Real adot = m_cosmology.adot(m_time);
  return m_Cexp * m_dx * a / adot;
}

/*******/
void
AMRLevelPIC::
levelSetup()
{
  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::levelSetup " << m_level << endl;
    }

  AMRLevelPIC* amrPartCoarserPtr = getCoarserLevel();
  AMRLevelPIC* amrPartFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrPartCoarserPtr != NULL);
  m_hasFiner   = (amrPartFinerPtr   != NULL);
  m_isThisFinestLev = m_grids.size() > 0;

  int nRefCrse = -1;

  if (m_hasCoarser)
    {
      amrPartCoarserPtr = dynamic_cast<AMRLevelPIC*>(m_coarser_level_ptr);

      nRefCrse = amrPartCoarserPtr->refRatio();
      const DisjointBoxLayout& crseGrids = amrPartCoarserPtr->m_grids;

      m_PVR.define(m_grids, amrPartCoarserPtr->m_PVR.mask(), nRefCrse, m_buffer);

      const int numComp = 1;
      m_coarseAveragePhi.define(m_grids,
                                crseGrids,
                                numComp,
                                nRefCrse);

      m_quadCFInterpPhi.define(m_grids, &crseGrids, m_dx, 
                               nRefCrse, 1, m_problem_domain);

      m_fineInterpRhs.define(m_grids, 1, nRefCrse, m_problem_domain);
    }

  else
    {
      m_PVR.define(m_grids, NULL, 1, 0);
    }

  m_levelParticle.define(m_problem_domain,
                         m_dx,
                         m_patchParticle);
}

void 
AMRLevelPIC::
aggregateParticles(ParticleData<JointParticle>& a_particles)
{

  CH_TIME("AMRLevelPIC::aggregateParticles");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::aggregateParticles " << m_level << endl;
    }

  a_particles.clear();
  List<JointParticle>& plist = a_particles.outcast();

  // coarse dx
  const RealVect cdx = m_meshSpacing * getCoarserLevel()->m_ref_ratio;

  for (DataIterator di(m_grids); di.ok(); ++di)
    {
      // map part pos to bins
      std::map<IntVect,JointParticle,CompIntVect> mip;

      binmapParticles(mip,m_PNew[di].listItems(), cdx, m_origin);

      if (m_hasFiner)
        {
          binmapParticles(mip,m_jointParticle[di].listItems(), cdx, m_origin);
        }

      for (map<IntVect,JointParticle,CompIntVect>::iterator it=mip.begin(); it!=mip.end(); ++it)
        {
          plist.add(it->second);
        }
    }

  // now do outcasts
  {
    // map part pos to bins
    std::map<IntVect,JointParticle,CompIntVect> mip;

    binmapParticles(mip,m_PNew.outcast(), cdx, m_origin);

    for (map<IntVect,JointParticle,CompIntVect>::iterator it=mip.begin(); it!=mip.end(); ++it)
      {
        plist.add(it->second);
      }
  }

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::aggregateParticles done " << endl;
    }

}

void AMRLevelPIC::amrSync(const int a_baseLevel, const bool a_srceCorr)
{
  CH_TIME("AMRLevelPIC::amrSync");

#ifdef CH_MPI
  {
    CH_TIME("amrSync::start-barrier");
    MPI_Barrier( Chombo_MPI::comm );
  }
#endif

  if (m_verbosity >= 3)
    {
      pout() << " AMRLevelPIC::amrSync start: baseLevel = " << a_baseLevel << endl;
    }

  ellipticSolver(a_baseLevel, false);

  Vector<AMRLevelPIC*> amr_levels = getAMRHierarchy();

  // need to do the two calls in this order, to ensure that boundary
  // conditions can be properly applied
  if (a_srceCorr)
    {
      // applies 2nd order corrections and compute the new force
      for (int l=a_baseLevel; l<amr_levels.size(); ++l)
        {
          amr_levels[l]->secondOrderCorrection();
        }
    }
  else
    {
      for (int l=a_baseLevel; l<amr_levels.size(); ++l)
        {

          amr_levels[l]->computeForce(amr_levels[l]->m_fieldNew,
                                      amr_levels[l]->m_phiNew,
                                      m_time);
        }
    }

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::amrSync done " << endl;
    }
}

// Allow for second order source term corrections
void AMRLevelPIC::secondOrderCorrection()
{

  CH_TIME("AMRLevelPIC::secondOrderCorrection");

#ifdef CH_MPI
  {
    CH_TIME("secondOrderCorrection::start-barrier");
    MPI_Barrier( Chombo_MPI::comm );
  }
#endif

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::secondOrderCorrection: " << m_level << endl;
    }

  // new force ...
  computeForce(m_fieldNew, m_phiNew, m_time);

  // get the patchParticle object with coefficients already set
  PatchParticle* patchParticle = m_levelParticle.getPatchParticlePtr();

  for (DataIterator di(m_grids); di.ok(); ++di)
    {
      const Box& box = m_grids.get(di);
      patchParticle->correctState(m_PNew[di], m_fieldNew[di], m_dt, box);
    }

  // exchange+BC's and resetting of coeff's
  m_levelParticle.postCorrect(m_PNew);

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::secondOrderCorrection: done " << endl;
    }
}

void 
AMRLevelPIC::
defineEllipticSolvers(const int a_baseLevel)
{

  CH_TIME("AMRLevelPIC::defineEllipticSolvers");

  if (m_verbosity >= 3)
    {
      pout() << " AMRLevelPIC::defineEllipticSolver: " << a_baseLevel << endl;
    }

  //
  const int finestLevel = a_baseLevel;

  // setup level pointers and quantities up to finest level
  Vector<DisjointBoxLayout>  amrBoxes(finestLevel+1);
  Vector<int>                amrRefRat(finestLevel+1,1);

  Vector<AMRLevelPIC*> amrLevels = getAMRHierarchy();
  CH_assert(amrLevels.size() >= finestLevel+1);

  for (int lev=finestLevel; lev>=0; --lev)
    {
      amrBoxes[lev] = amrLevels[lev]->m_grids;
      amrRefRat[lev]= amrLevels[lev]->refRatio();
    }

  const Real dx0 = amrLevels[0]->m_dx;
  const ProblemDomain probDom0 = amrLevels[0]->problemDomain();

  s_gravOpFact->define(probDom0,amrBoxes,amrRefRat,dx0,&NoOpBc);

  s_gravMGSolver->define(probDom0,*s_gravOpFact,s_bottomSolver,finestLevel+1);
  s_gravMGSolver->setSolverParameters(SOLVER_NUM_SMOOTH,SOLVER_NUM_SMOOTH,
                                      SOLVER_NUM_SMOOTH,SOLVER_NUM_MG,
                                      SOLVER_MAX_ITER,SOLVER_TOLERANCE,
                                      SOLVER_HANG,SOLVER_NORM_THRES);
  s_gravMGSolver->m_verbosity = 3;
}

Vector<AMRLevelPIC*> 
AMRLevelPIC::
getAMRHierarchy()
{
  Vector<AMRLevelPIC*> retval;

  // First go to level 0
  AMRLevelPIC* levelPtr = this;

  while (levelPtr->hasCoarserLevel())
    {
      levelPtr = levelPtr->getCoarserLevel();
    }

  // Now can accumulate the pointers by chasing finer level
  retval.push_back(levelPtr);
  while (levelPtr->hasFinerLevel())
    {
      levelPtr = levelPtr->getFinerLevel();
      retval.push_back(levelPtr);
    }

  return retval;
}

LevelData<FArrayBox>* AMRLevelPIC::getPhi(const Real& a_time)
{

  if (m_verbosity >= 3)
    {
      pout() << " AMRLevelPIC::getPhi " << m_level << " a_time= " << a_time << " m_time " << m_time << endl;
    }

  LevelData<FArrayBox>* phiPtr = NULL;

  const Real eps = 0.01 * m_dt;
  if (Abs(a_time - m_time) <= eps) // case alpha=1; new synchronization point
    {
      phiPtr = &m_phiNew;
    }
  else if (Abs(a_time-(m_time-m_dtOld)) < eps) // case alpha=0; old synch point
    {
      phiPtr = &m_phiOld;
    }
  else
    {
      // need time interpolation
      interpolateInTime(m_phiInt, m_phiOld, m_phiNew,
                        a_time, m_time, m_time - m_dtOld, m_dt,
                        Interval(0,0),Interval(0,0));

      phiPtr = &m_phiInt;
    }

  if (phiPtr == NULL)
    {
      MayDay::Error("AMRLevelPIC::getPhi: something failed");
    }

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC:: getPhi done " << endl;
    }

  return phiPtr;
}

LevelData<FArrayBox>* AMRLevelPIC::getPoissonRhs(const Real& a_time)
{
  
  if (m_verbosity >= 3)
    {
      pout() << " AMRLevelPIC:: getPoissonRhs: level " << m_level << endl;
    }

  makePoissonRhs( m_rhs, m_PNew, m_jointParticle, a_time);

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::getPoissonRhs done " << endl;
    }

  return &m_rhs;
}

// function that computes F=-grad(phi)
void AMRLevelPIC::computeForce(LevelData<FArrayBox>&  a_force,
                                        LevelData<FArrayBox>&  a_phi,
                                        const Real&            a_time)
{

  CH_TIME("AMRLevelPIC::computeForce");

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::computeForce " << m_level << endl;
    }

  a_phi.exchange();

  // fill ghost cells abutting coarser level
  if (m_hasCoarser)
    {
      AMRLevelPIC* coarserLevelPtr = getCoarserLevel();

      // get potential from coarser level
      const LevelData<FArrayBox>* crsePhi = coarserLevelPtr->getPhi(a_time);

      // fill ghost cells on fine level
      m_quadCFInterpPhi.coarseFineInterp(a_phi, *crsePhi);
    }

  // compute negative gradient
  SingleLevelGradient::levelGradientCC(a_force, a_phi, -m_dx);

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::computeForce done "  << endl;
    }

  a_force.exchange();
}

//
void AMRLevelPIC::ellipticSolver(const int  a_baseLevel,
                                          const bool a_isLevelSolve)
{
  CH_TIME("AMRLevelPIC::ellipticSolver");
  CH_assert(a_baseLevel==m_level);

  if (m_verbosity >= 3)
    {
      pout() << " AMRLevelPIC::ellipticSolver start: baseLevel = " << a_baseLevel;
    }

  // get pointers to all the levels in the hierarchy
  Vector<AMRLevelPIC*> amrLevels = getAMRHierarchy();

  // If we're doing a level solve, the finest level and the 
  // base level are the same. If we're doing a full multilevel
  // solve, we want to go all the way to the top of the hierarchy.
  int finestLevel = m_level;
  if (!a_isLevelSolve) finestLevel=amrLevels.size() - 1;
  const int maxLevel = finestLevel;

  if (m_verbosity >= 3)
    {
      pout() << "... and maxLevel = " << maxLevel << endl;
    }

  // Pointers to the answer
  Vector<LevelData<FArrayBox>*> amrPhi(finestLevel+1,NULL);
  const int bndryLevel = Max(0, a_baseLevel - 1);
  for (int lev = bndryLevel; lev <= finestLevel; ++lev)
    {
      amrPhi[lev] = amrLevels[lev]->getPhi(m_time);
    }

  // get refinement ratios for all required levels
  Vector<int> amrRefRat(maxLevel+1, 1);
  for (int lev = maxLevel; lev >= 0; --lev)
  {
    amrRefRat[lev]  = amrLevels[lev]->refRatio();
  }

  // setup Poisson eq. right hand side
  Vector<LevelData<FArrayBox>*> amrRhs(maxLevel+1, NULL);
  for (int lev = a_baseLevel; lev <= maxLevel; ++lev)
    {
      amrRhs[lev]=amrLevels[lev]->getPoissonRhs(m_time);
    }

  // if we are using periodic boundary conditions, we have to ensure 
  // that the sum of the RHS is zero. We do this here. 
  bool isDomainCovered = (a_baseLevel == 0);
  if (m_problem_domain.isPeriodic())
    {
      // should be: sum(rhs)=0; if not, then offset residual
      // and synchronize
      if (!isDomainCovered)
        {
          //long numPtsDomain = amrProbDom[a_baseLevel].domainBox().numPts();
          long long numPtsDomain = m_problem_domain.domainBox().numPts();

          // count number of cells on this level
          long numPtsLevel = 0;
          //const DisjointBoxLayout& baseGrids = amrBoxes[a_baseLevel];
          for (LayoutIterator lit = m_grids.layoutIterator(); lit.ok(); ++lit)
            {
              numPtsLevel += m_grids[lit].numPts();
            }

          isDomainCovered = (numPtsDomain==numPtsLevel);
        }

      // compute sum(rhs)
      if (!a_isLevelSolve)
        {
          if (isDomainCovered)
            {
              m_rhsOffset = computeSum(amrRhs, amrRefRat,
                                       m_dx,
                                       Interval(0,0),
                                       a_baseLevel);

              if (m_verbosity >= 2)
                {
                   pout() << " gravity:: lev: " << a_baseLevel
                          << " rhs_resid(a): " << m_rhsOffset << endl;
                }

              for (int lev = maxLevel; lev > a_baseLevel; --lev)
                {
                  amrLevels[lev]->m_rhsOffset = m_rhsOffset;
                }
              }
          else
              {
                // we needed this if the levels were created after
                // rhsOffset was computed
                for (int lev = maxLevel; lev > bndryLevel; --lev)
                  {
                    amrLevels[lev]->m_rhsOffset = amrLevels[bndryLevel]->m_rhsOffset;
                  }
              }

          // enforce sum(rhs)=0
          for (int lev = a_baseLevel; lev <= maxLevel; ++lev)
            {
              Real domainVolume = pow(m_domainLength, CH_SPACEDIM);
              offset(*amrRhs[lev], m_rhsOffset / domainVolume);
            }

          if (m_verbosity >= 2)
            {
              if (isDomainCovered)
                {
                  // compute the new sum(rhs)
                  pout() << " gravity:: lev: " << a_baseLevel
                         << " rhs_resid(b): " << computeSum(amrRhs,
                                                            amrRefRat,
                                                            m_dx,
                                                            Interval(0,0),
                                                            a_baseLevel)
                         << endl;
                }
            }
          }

    else if (isDomainCovered) // if it *is* a level solve...
      {
        Real totalRHS = computeSum(*amrRhs[a_baseLevel], NULL, 1, m_dx);
        Real domainVolume = pow(m_domainLength, CH_SPACEDIM);
        offset(*amrRhs[a_baseLevel], totalRHS/domainVolume);
      }
    }

  bool reset = (a_baseLevel==0 ? true : SOLVER_RESET_PHI);

  // solve
  s_gravMGSolver->solve(amrPhi,amrRhs,maxLevel,a_baseLevel,reset);

  if (!a_isLevelSolve)
    {
      // do postelliptic kind of operations
      for (int lev = maxLevel; lev > a_baseLevel; --lev)
        {
          // Average Phi from finer level data
          amrLevels[lev]->m_coarseAveragePhi.averageToCoarse(*amrPhi[lev-1],
                                                             *amrPhi[lev]);

          amrLevels[lev]->m_coarseAveragePhi.averageToCoarse(*amrRhs[lev-1],
                                                             *amrRhs[lev]);
        }
    }

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::ellipticSolver done " << endl;
    }
}

//
void AMRLevelPIC::makePoissonRhs(LevelData<FArrayBox>&       a_rhs,
                                          const ParticleData<Particle>&  a_P,
                                          const ParticleData<JointParticle>& a_jointP,
                                          const Real&                 a_time)
{

  CH_TIME("AMRLevelPIC::makePoissonRhs");

  // I am assuming they are all defined on the same grids; assert()
  // this along the way

  if (m_verbosity >= 3)
    {
      pout() << "AMRLevelPIC::makePoissonRhs: level " << m_level << endl;
    }

  const DisjointBoxLayout& grids = a_rhs.getBoxes();
  DataIterator di = grids.dataIterator();

  // 0. set rhs = 0
  setToVal(a_rhs, 0.0);

  // 1. add contrib. from coarser lev particles whose clouds hang on this lev
  if (m_hasCoarser)
    {
      CH_TIME("makePoissonRhs::coarserhs");

      const AMRLevelPIC* coarserLevelPtr = getCoarserLevel();
      const Real tNew = coarserLevelPtr->m_time;
      const Real tOld = tNew - coarserLevelPtr->m_dtOld;

      // temporary storage for the coarse RHS interpolated to the current time
      const DisjointBoxLayout& crseFineBoxes = m_rhsCrseNew.getBoxes();
      LevelData<FArrayBox> rhsCrse(crseFineBoxes, 1);
      setToVal(rhsCrse, 0.0);

      // interpolate in time
      interpolateInTime(rhsCrse,m_rhsCrseOld,m_rhsCrseNew,
                        a_time,tNew,tOld,m_dt,
                        Interval(0,0),Interval(0,0), IntVect::Zero);     

      m_fineInterpRhs.pwcinterpToFine(a_rhs, rhsCrse);
    }

  // Deposit particles
  if (a_P.isDefined())
    {
      CH_assert(a_P.isClosed());
      CH_TIME("makePoissonRhs::particles");

      for (DataIterator di(m_grids); di.ok(); ++di)
        {
          m_patchParticle->deposit(a_rhs[di], a_P[di], m_grids[di]);
        }
    }

  // exchange ghost cells
  LDaddOp<FArrayBox> addOp;
  a_rhs.exchange(a_rhs.interval(), m_reverseCopier, addOp);

  // save this rhs into the coarseRhs for the next finest level
  if (m_hasFiner && !m_isThisFinestLev)
  {
      AMRLevelPIC* amrFinerPtr = getFinerLevel();
      LevelData<FArrayBox>& rhsFineNew = amrFinerPtr->m_rhsCrseNew;
      LevelData<FArrayBox>& rhsFineOld = amrFinerPtr->m_rhsCrseOld;

      // save data
      for (DataIterator di = rhsFineNew.dataIterator(); di.ok(); ++di)
        {
          rhsFineOld[di].copy(rhsFineNew[di]);
          rhsFineNew[di].setVal(0.0);
        }

      // new values: takes into account even contribution from ghost cells.
      // rhs_tmp.addTo(Interval(0,0),rhsFineNew,Interval(0,0),m_problem_domain);
      Copier fineRhsCopier;
      fineRhsCopier.define(getFinerLevel()->m_rhsCrseNew.getBoxes(),
                               m_grids,
                               m_problem_domain);
      fineRhsCopier.reverse();
      a_rhs.addTo(Interval(0,0),rhsFineNew,Interval(0,0),
                  m_problem_domain, fineRhsCopier);
  }

  // Handle particles on finer levels through Joint Particles
  if (m_hasFiner and a_jointP.isDefined())
  {
    CH_TIME("makePoissonRhs::jointParticles");

    for (DataIterator di(m_grids); di.ok(); ++di)
      {
        m_patchParticle->deposit(a_rhs[di], a_jointP[di], m_grids[di]);
      }
  }

  //a_rhs.exchange(a_rhs.interval(), m_reverseCopier, addOp);

  // cosmic expansion stuff here
  const Real H0 = m_cosmology.m_H0;
  const Real rhsScale = 1.5*H0*H0/m_cosmology.a(a_time);

  for (DataIterator dit(m_grids); dit.ok(); ++dit)
  {
    a_rhs[dit()]   -= 1.0;
    a_rhs[dit()] *= rhsScale;
  }

  a_rhs.exchange();

  if (m_verbosity >= 3)
    {
      pout() << " AMRLevelPIC::makePoissonRhs: done  " << endl;
    }
}

/*******/
AMRLevelPIC*
AMRLevelPIC::
getCoarserLevel() const
{
  AMRLevelPIC* amrPartCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      amrPartCoarserPtr = dynamic_cast<AMRLevelPIC*>(m_coarser_level_ptr);

      if (amrPartCoarserPtr == NULL)
        {
          MayDay::Error("AMRLevelPIC::getCoarserLevel: dynamic cast failed");
        }
    }

  return amrPartCoarserPtr;
}

/*******/
AMRLevelPIC*
AMRLevelPIC::
getFinerLevel() const
{
  AMRLevelPIC* amrPartFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
    {
      amrPartFinerPtr = dynamic_cast<AMRLevelPIC*>(m_finer_level_ptr);

      if (amrPartFinerPtr == NULL)
        {
          MayDay::Error("AMRLevelPIC::getFinerLevel: dynamic cast failed");
        }
    }

  return amrPartFinerPtr;
}

/*******/
int AMRLevelPIC::countItems() const
{
  return m_PNew.numValid();
}

/*******/
int AMRLevelPIC::countOutcast() const
{
  return m_PNew.numOutcast();
}

#include "NamespaceFooter.H"
