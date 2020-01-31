#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELTRACER_H_
#define _AMRLEVELTRACER_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"

#include "ListBox.H"
#include "ListBoxFactory.H"
#include "Particle.H"
#include "AdvectionFunctions.H"
#include "ParticleData.H"
#include "ParticleValidRegion.H"

#include "ProblemDomain.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"

#include "NamespaceHeader.H"

/// AMRLevel for tracer particles
/**
 */
class AMRLevelTracer : public AMRLevel
{
public:

  /// Default constructor
  AMRLevelTracer()
  {
    m_isDefined = false;
  }

  /// Full constructor. Arguments are same as in define()
  AMRLevelTracer(AdvectionVelocityFunction a_advFunc,
		 const Real& a_cfl,
		 const Real& a_domainLength,
		 const int& a_maxBoxSize)
  {
    define(a_advFunc, a_cfl, a_domainLength, a_maxBoxSize);
  }

  /// Defines this AMRLevelTracer
  void define(AdvectionVelocityFunction a_advFunc,
              const Real& a_cfl,
              const Real& a_domainLength,
              const int& a_maxBoxSize);

  /// Destructor
  virtual ~AMRLevelTracer();

  /// Never called: historical
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio)
  {
    MayDay::Error("never called--historical");
  }

  /// Define new AMRLevelTracer from coarser
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

  // Setup menagerie of data structures
  void levelSetup();

  // the total number of particles on this level
  int countItems() const;

  // total number of outcast particles
  int countOutcast() const;

  // Get the next coarser level
  AMRLevelTracer* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelTracer* getFinerLevel() const;

  // Particles at the old and new times
  ParticleData<Particle> m_POld, m_PNew;

  bool                 m_isDefined;
  Real                 m_cfl;
  Real                 m_domainLength;
  Real                 m_dx;
  RealVect             m_meshSpacing;
  RealVect             m_origin;
  AdvectionVelocityFunction    m_advFunc;
  Real                 m_maxVelocityLocal;

  Real                 m_dtNew;
  int                  m_numGhost;
  int                  m_maxBoxSize;
  const static int     m_buffer = 2;
  const static int     m_verbosity = 1;

  ParticleValidRegion  m_PVR;

  bool                 m_hasCoarser;
  bool                 m_hasFiner;

  DisjointBoxLayout    m_grids;

private:

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelTracer&);
  AMRLevelTracer(const AMRLevelTracer&);
};

#include "NamespaceFooter.H"

#endif
