#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELPIC_H_
#define _AMRLEVELPIC_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "ProblemDomain.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "QuadCFInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "SingleLevelGradient.H"

// particle includes
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "Particle.H"
#include "ParticleData.H"
#include "ParticleValidRegion.H"
#include "MeshInterp.H"
#include "LevelParticle.H"
#include "JointParticle.H"
#include "Cosmology.H"
#include "ZeldovichIBC.H"

// solver stuff
#include "AMRMultiGrid.H"
#include "AMRPoissonOp.H"
#include "BiCGStabSolver.H"

#include "Convenience.H"

#include "NamespaceHeader.H"

/// AMRLevel for Particles
/**
 */
class AMRLevelPIC : public AMRLevel
{
public:

  /// Default constructor
  AMRLevelPIC()
  {
    m_patchParticle = NULL;
    m_ibc = NULL;
    m_isDefined = false;
  }

/// Full constructor. Arguments are same as in define()
AMRLevelPIC(const ZeldovichIBC* a_ibc,
	    const Real&         a_cfl,
	    const Real&         a_domainLength,
	    const int&          a_maxBoxSize,
	    const InterpType&   a_poissonInterpType,
	    const Real&         a_refineThresh)
  {
    
    m_patchParticle = NULL;
    m_ibc = NULL;
    m_isDefined = false;    
    
    define(a_ibc, 
           a_cfl, 
           a_domainLength, 
           a_maxBoxSize, 
           a_poissonInterpType,
           a_refineThresh);
  }
  
/// Defines this AMRLevelPIC
void define(const ZeldovichIBC* a_ibc,
	    const Real&         a_cfl,
	    const Real&         a_domainLength,
	    const int&          a_maxBoxSize,
	    const InterpType&   a_poissonInterpType,
	    const Real&         a_refineThresh);

/// Destructor
virtual ~AMRLevelPIC();

/// Never called: historical
virtual void define(AMRLevel*  a_coarserLevelPtr,
		    const Box& a_problemDomain,
		    int        a_level,
		    int        a_refRatio)
{
  MayDay::Error("never called--historical");
}

/// Define new AMRLevelPIC from coarser
virtual void define(AMRLevel*            a_coarserLevelPtr,
		    const ProblemDomain& a_problemDomain,
		    int                  a_level,
		    int                  a_refRatio);

/// Advance by one timestep
virtual Real advance();

/// Things to do after a timestep
virtual void postTimeStep();

/// Create tags for regridding
virtual void tagCells(IntVectSet& a_tags) ;

/// Create tags at initialization
virtual void tagCellsInit(IntVectSet& a_tags) ;

/// Set up data on this level before regridding
virtual void preRegrid(int a_base_level,
                       const Vector<Vector<Box> >& a_newGrids);

/// Set up data on this level after regridding
virtual void regrid(const Vector<Box>& a_newGrids);

/// postRegridding ops
virtual void postRegrid(int a_baseLevel);

/// Initialize grids
virtual void initialGrid(const Vector<Box>& a_newGrids);

/// Initialize data
virtual void initialData();

/// Things to do after initialization
virtual void postInitialize();

#ifdef CH_USE_HDF5
/// Write checkpoint header
virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

/// Write checkpoint data for this level
virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

/// Read checkpoint header
virtual void readCheckpointHeader(HDF5Handle& a_handle);

/// Read checkpoint data for this level
virtual void readCheckpointLevel(HDF5Handle& a_handle);

/// Write plotfile header
virtual void writePlotHeader(HDF5Handle& a_handle) const;

/// Write plotfile data for this level
virtual void writePlotLevel(HDF5Handle& a_handle) const;
#endif

/// Returns the dt computed earlier for this level
virtual Real computeDt();

/// Compute dt using initial data
virtual Real computeInitialDt();

protected:

/// Setup menagerie of data structures
void levelSetup();

/// collect Particles from finer levels into the Joint particle list
void aggregateParticles(ParticleData<JointParticle>& a_particles);

/// amrSync performs a full, multilevel Poisson solve starting on
/// a_baseLevel and including all finer levels. This both computes
/// the composite potential and updates the forces. If a_srceCorr is 
/// true, the resulting forces will also be used to apply 2nd-order
/// corrections to the particle positions and velocities on all affected
/// levels. 
void amrSync(const int a_baseLevel, const bool a_srceCorr=true);

/// Computes the force associated with the current potential and uses
/// it to apply second order corrections to the particle positions and 
/// velocities on this level.
void secondOrderCorrection();

/// Returns the total number of valid particles on this level, not
/// counting outcasts.
int countItems() const;

/// Returns the total number of outcast particles on this level.
int countOutcast() const;

/// define gravity solver
void defineEllipticSolvers(const int a_baseLevel);

/// return pointer to Phi, interpolated to a_time
LevelData<FArrayBox>* getPhi(const Real& a_time);

/// build the RHS for the Poisson equation and return
/// a pointer to it, interpolated to a_time
LevelData<FArrayBox>* getPoissonRhs(const Real& a_time);

/// compute a force from a potential
void computeForce(LevelData<FArrayBox>& a_force,
		  LevelData<FArrayBox>& a_phi,
		  const Real&           a_time);

/// actually solve Poisson's equation. If a_isLevelSolve is true, 
/// then this will only do this level. If it is false, it will do
/// this and all finder levels.
void ellipticSolver(const int  a_baseLevel,
		    const bool a_isLevelSolve);

/// This generates the RHS for Poisson's equations from the particles
/// on this level. It also includes a contribution for particles on 
/// coarser levels whose clouds hang over this level's grids. 
void makePoissonRhs(LevelData<FArrayBox>&          a_rhs,
		    const ParticleData<Particle>&  a_P,
		    const ParticleData<JointParticle>& a_jointP,
		    const Real&                    a_time);

// Returns a pointer to the AMR Level Hierarchy
Vector<AMRLevelPIC*> getAMRHierarchy();

// Get the next coarser level
AMRLevelPIC* getCoarserLevel() const;

// Get the next finer level
AMRLevelPIC* getFinerLevel() const;

bool                 m_isDefined;
Real                 m_cfl;
Real                 m_Cexp;
Real                 m_domainLength;
Real                 m_dx;
Real                 m_refineThresh;
RealVect             m_meshSpacing;
RealVect             m_origin;
Real                 m_rhsOffset;
Real                 m_dtNew, m_dtOld;
Real                 m_numForceGhost;
int                  m_maxBoxSize;
InterpType           m_poissonInterpType;
const static int     m_buffer = 2;
const static int     m_verbosity = 3;

// m_hasFiner tells you whether there are *allowed* 
// to be finer or levels in the AMR hierarchy. These
// levels do not neccessarily have to have any data yet,
// if, say, the refinement criteria have not yet required
// that a finer level be used. m_isThisFinestLev tracks 
// whether this the finest level that actually contains 
// data.
bool                 m_hasCoarser;
bool                 m_hasFiner;
bool                 m_isThisFinestLev;
DisjointBoxLayout    m_grids;

// Particles at the old and new times
ParticleData<Particle> m_POld, m_PNew;

// Joint Particles that represent the net effect of the particles
// that live on finer levels
ParticleData<JointParticle> m_jointParticle;

// Handles expansion of the universe
Cosmology m_cosmology;

// Initial conditions
ZeldovichIBC* m_ibc;

// Used to send the particles to the correct AMR level.
ParticleValidRegion  m_PVR;

// Right-hand side of the Poisson equation
LevelData<FArrayBox> m_rhs;      

// the contribution to the RHS from particles that live on the next
// coarser level, at the old and new times.
LevelData<FArrayBox> m_rhsCrseNew, m_rhsCrseOld;

// The gravitational potential at the Old, New, and an intermediate time
LevelData<FArrayBox> m_phiOld, m_phiNew, m_phiInt;

// The gravitational fields
LevelData<FArrayBox> m_fieldOld, m_fieldNew;

// Averaging the gravitational potential from fine to coarse level
CoarseAverage m_coarseAveragePhi;

// Filling in ghost cells for the potential using quadratic interpolation
QuadCFInterp m_quadCFInterpPhi;

// Used to add in the contributions from coarser particle's clouds in to
// the RHS for this level.
FineInterp m_fineInterpRhs;

// Level integrator for particles
LevelParticle m_levelParticle;

// Patch Integrator for particles
PatchParticle* m_patchParticle;
  
/// We have two different Copiers for ghost exchange operations. 
/// The first is used for the normal exchange, which copies from valid
/// to ghost cells. The second reverses this pattern of communication,
/// moving data FROM ghost TO valid cells. This is needed during the 
/// mass deposition step. 
Copier               m_forwardCopier;
Copier               m_reverseCopier;

// gravity solver stuff
static BiCGStabSolver<LevelData<FArrayBox> >* s_bottomSolver;
static RefCountedPtr<AMRPoissonOpFactory> s_gravOpFact;
static RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > s_gravMGSolver;

private:

// Disallowed for all the usual reasons
void operator=(const AMRLevelPIC&);
AMRLevelPIC(const AMRLevelPIC&);

};

#include "NamespaceFooter.H"

#endif
