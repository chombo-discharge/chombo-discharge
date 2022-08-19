/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticleSolverImplem.H
  @brief  Implementation of CD_TracerParticleSolver.H
  @author Robert Marskar
*/

#ifndef CD_TracerParticleSolverImplem_H
#define CD_TracerParticleSolverImplem_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_TracerParticleSolver.H>
#include <CD_DischargeIO.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

template <typename P>
inline TracerParticleSolver<P>::TracerParticleSolver()
{
  CH_TIME("TracerParticleSolver::TracerParticleSolver()");

  // Default settings
  m_verbosity     = -1;
  m_name          = "TracerParticleSolver";
  m_className     = "TracerParticleSolver";
  m_realm         = Realm::primal;
  m_phase         = phase::gas;
  m_pvrBuffer     = 0;
  m_haloBuffer    = 1;
  m_plotVelocity  = true;
  m_plotMass      = true;
  m_volumeScale   = false;
  m_deposition    = DepositionType::CIC;
  m_velocityField = nullptr;
}

template <typename P>
inline TracerParticleSolver<P>::TracerParticleSolver(const RefCountedPtr<AmrMesh>&              a_amr,
                                                     const RefCountedPtr<ComputationalGeometry> a_compGeom)
  : TracerParticleSolver<P>()
{
  CH_TIME("TracerParticleSolver::TracerParticleSolver(RefCountedPtr<AmrMesh>, RefCountedPtr<ComputationalGeometry>)");

  m_amr                   = a_amr;
  m_computationalGeometry = a_compGeom;
}

template <typename P>
inline TracerParticleSolver<P>::~TracerParticleSolver()
{
  CH_TIME("TracerParticleSolver::~TracerParticleSolver()");
}

template <typename P>
inline void
TracerParticleSolver<P>::parseOptions()
{
  CH_TIME("TracerParticleSolver::parseOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions()" << endl;
  }

  this->parseDeposition();
  this->parsePlotVariables();
  this->parseVerbosity();
}

template <typename P>
inline void
TracerParticleSolver<P>::parseRuntimeOptions()
{
  CH_TIME("TracerParticleSolver::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }

  this->parseDeposition();
  this->parsePlotVariables();
  this->parseVerbosity();
}

template <typename P>
inline void
TracerParticleSolver<P>::parseDeposition()
{
  CH_TIME("TracerParticleSolver::parseDeposition()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDeposition()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  // Deposition for particle-mesh operations
  pp.get("deposition", str);
  if (str == "ngp") {
    m_haloBuffer = 0;

    m_deposition           = DepositionType::NGP;
    m_coarseFineDeposition = CoarseFineDeposition::Interp;
  }
  else if (str == "cic") {
    m_haloBuffer = 1;

    m_deposition           = DepositionType::CIC;
    m_coarseFineDeposition = CoarseFineDeposition::Halo;
  }
  else {
    MayDay::Error("TracerParticleSolver::parseDeposition - unknown deposition method requested.");
  }

  pp.get("volume_scale", m_volumeScale);
}

template <typename P>
inline void
TracerParticleSolver<P>::parsePlotVariables()
{
  CH_TIME("TracerParticleSolver::parsePlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("plot_mass", m_plotMass);
  pp.get("plot_velocity", m_plotVelocity);
}

template <typename P>
inline void
TracerParticleSolver<P>::parseVerbosity()
{
  CH_TIME("TracerParticleSolver::parseVerbosity()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseVerbosity()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("verbosity", m_verbosity);
}

template <typename P>
inline void
TracerParticleSolver<P>::registerOperators() const
{
  CH_TIME("TracerParticleSolver::registerOperators()");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators()" << endl;
  }

  if (m_amr.isNull()) {
    MayDay::Abort("TracerParticleSolver::registerOperators - need to set AmrMesh!");
  }
  else {
    m_amr->registerOperator(s_particle_mesh, m_realm, m_phase);
    m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);

    m_amr->registerMask(s_particle_halo, 1, m_realm); // For CIC halo deposition.
  }
}

template <typename P>
inline void
TracerParticleSolver<P>::allocate()
{
  CH_TIME("TracerParticleSolver::allocate()");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate()" << endl;
  }

  // Allocate data for storing the velocity field on the mesh, and storage for the particles.
  m_amr->allocate(m_particles, m_realm);
}

template <typename P>
inline void
TracerParticleSolver<P>::setRealm(const std::string& a_realm)
{
  CH_TIME("TracerParticleSolver::setRealm(std::string)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRealm(std::string)" << endl;
  }

  m_realm = a_realm;
}

template <typename P>
inline void
TracerParticleSolver<P>::setPhase(const phase::which_phase& a_phase)
{
  CH_TIME("TracerParticleSolver::setPhase(phase::which_phase)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase(phase::which_phase)" << endl;
  }

  m_phase = a_phase;
}

template <typename P>
void
TracerParticleSolver<P>::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("TracerParticleSolver::setTime(int, Real, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setTime(int, Real, Real)" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

template <typename P>
inline void
TracerParticleSolver<P>::setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh)
{
  CH_TIME("TracerParticleSolver::setAmr(RefCountedPtr<AmrMesh>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setAmr(RefCountedPtr<AmrMesh>)" << endl;
  }

  m_amr = a_amrMesh;
}

template <typename P>
inline void
TracerParticleSolver<P>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_compGeom)
{
  CH_TIME("TracerParticleSolver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  m_computationalGeometry = a_compGeom;
}

template <typename P>
inline void
TracerParticleSolver<P>::setVelocity(const EBAMRCellData& a_velocityField)
{
  CH_TIME("TracerParticleSolver::setVelocity(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVelocity(EBAMRCellData)" << endl;
  }

  m_velocityField = &a_velocityField;
}

template <typename P>
inline void
TracerParticleSolver<P>::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("TracerParticleSolver::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid(int, int)" << endl;
  }

  CH_assert(a_lbase >= 0);

  // Put particles in pre-regrid mode.
  m_particles.preRegrid(a_lbase);
}

template <typename P>
inline void
TracerParticleSolver<P>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("TracerParticleSolver::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }

  CH_assert(a_lmin >= 0);
  CH_assert(a_oldFinestLevel >= 0);
  CH_assert(a_newFinestLevel >= 0);

  // Regrid the particles.
  const Vector<DisjointBoxLayout>& grids   = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>&     domains = m_amr->getDomains();
  const Vector<Real>&              dx      = m_amr->getDx();
  const Vector<int>&               refRat  = m_amr->getRefinementRatios();

  m_particles.regrid(grids, domains, dx, refRat, a_lmin, a_newFinestLevel);
}

template <typename P>
inline void
TracerParticleSolver<P>::remap()
{
  CH_TIME("TracerParticleSolver::remap()");
  if (m_verbosity > 5) {
    pout() << m_name + "::remap()" << endl;
  }

  m_particles.remap();
}

template <typename P>
inline void
TracerParticleSolver<P>::deposit(EBAMRCellData& a_phi) const
{
  CH_TIME("TracerParticleSolver::deposit(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::deposit(EBAMRCellData)" << endl;
  }

  DataOps::setValue(a_phi, 0.0);
  this->depositParticles<P, &P::mass>(a_phi, m_particles, m_deposition, m_coarseFineDeposition);

  if (m_volumeScale) {
    DataOps::volumeScale(a_phi, m_amr->getDx());
  }
}

template <typename P>
inline void
TracerParticleSolver<P>::interpolateVelocities()
{
  CH_TIME("TracerParticleSolver::interpolateVelocities()");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateVelocities()" << endl;
  }

  CH_assert(!(m_velocityField == nullptr));

  m_amr->interpolateParticles<P, &P::velocity>(m_particles, m_realm, m_phase, *m_velocityField, m_deposition, false);
}

template <typename P>
inline int
TracerParticleSolver<P>::getNumberOfPlotVariables() const
{
  CH_TIME("TracerParticleSolver::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  int numPlotVars = 0;

  if (m_plotMass)
    numPlotVars += 1;
  if (m_plotVelocity)
    numPlotVars += SpaceDim;

  return numPlotVars;
}

template <typename P>
inline Vector<std::string>
TracerParticleSolver<P>::getPlotVariableNames() const
{
  CH_TIME("TracerParticleSolver::getPlotVariableNames()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames()" << endl;
  }

  // TLDR: Possible plot variables is the density (m_phi), diffusion coefficient, source term, velocity, and eb flux. This
  //       function returns the associated plot variable names, and will be used in the plot files.

  Vector<std::string> plotVarNames(0);

  if (m_plotMass) {
    plotVarNames.push_back(m_name + " density");
  }

  if (m_plotVelocity)
    plotVarNames.push_back("x-Velocity " + m_name);
  if (m_plotVelocity)
    plotVarNames.push_back("y-Velocity " + m_name);
#if CH_SPACEDIM == 3
  if (m_plotVelocity)
    plotVarNames.push_back("z-Velocity " + m_name);
#endif

  return plotVarNames;
}

template <typename P>
inline void
TracerParticleSolver<P>::writePlotFile()
{
  CH_TIME("TracerParticleSolver::writePlotFile()");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotFile()" << endl;
  }

  // Number of output components and their names
  const int                 numPlotVars  = this->getNumberOfPlotVariables();
  const Vector<std::string> plotVarNames = this->getPlotVariableNames();

  // Allocate storage
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, m_phase, numPlotVars);
  DataOps::setValue(output, 0.0);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0;
  this->writePlotData(output, icomp);

  // Filename
  char filename[100];
  sprintf(filename, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_timeStep, SpaceDim);
  std::string fname(filename);

  // Alias, because Chombo's EBAMRIO wants raw pointers (but we stick to smart pointers, like God intended).
  Vector<LevelData<EBCellFAB>*> outputPtr;
  m_amr->alias(outputPtr, output);

#ifdef CH_USE_HDF5
  constexpr int numPlotGhost = 0;

  DischargeIO::writeEBHDF5(fname,
                           plotVarNames,
                           m_amr->getGrids(m_realm),
                           outputPtr,
                           m_amr->getDomains(),
                           m_amr->getDx(),
                           m_amr->getRefinementRatios(),
                           m_dt,
                           m_time,
                           m_amr->getProbLo(),
                           m_amr->getFinestLevel() + 1,
                           numPlotGhost);
#endif
}

template <typename P>
inline void
TracerParticleSolver<P>::writePlotData(EBAMRCellData& a_output, int& a_icomp)
{
  CH_TIME("TracerParticleSolver::writePlotData(EBAMRCellData, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData(EBAMRCellData, int)" << endl;
  }

  // Write the particle scalars.
  if (m_plotMass) {
    EBAMRCellData mass;
    m_amr->allocate(mass, m_realm, m_phase, 1);

    this->deposit(mass);

    const Interval srcInterval(0, 0);
    const Interval dstInterval(a_icomp, a_icomp);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      mass[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
    }

    a_icomp++;
  }

  // Write the velocity field to the output data holder.
  if (m_plotVelocity) {
    CH_assert(!(m_velocityField == nullptr));

    const Interval srcInterval(0, SpaceDim - 1);
    const Interval dstInterval(a_icomp, a_icomp + SpaceDim - 1);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      (*m_velocityField)[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
    }

    a_icomp += SpaceDim;
  }
}

#ifdef CH_USE_HDF5
template <typename P>
inline void
TracerParticleSolver<P>::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("TracerParticleSolver::writeCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  writeParticlesToHDF(a_handle, m_particles[a_level], m_name + "_particles");
}
#endif

#ifdef CH_USE_HDF5
template <typename P>
inline void
TracerParticleSolver<P>::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("TracerParticleSolver::readCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevel(HDF5Handle, int)" << endl;
  }

  readParticlesFromHDF(a_handle, m_particles[a_level], m_name + "_particles");
}
#endif

template <typename P>
Real
TracerParticleSolver<P>::computeDt() const
{
  CH_TIME("TracerParticleSolver::computeDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDt()" << endl;
  }

  Real dt = std::numeric_limits<Real>::infinity();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real               dx  = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const List<P>& particles = m_particles[lvl][dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const RealVect& v = lit().velocity();

        for (int dir = 0; dir < SpaceDim; dir++) {
          dt = std::min(dt, dx / std::abs(v[dir]));
        }
      }
    }
  }

  // Get the minimum of dt over all MPI ranks.
  dt = ParallelOps::min(dt);

  return dt;
}

template <typename P>
inline ParticleContainer<P>&
TracerParticleSolver<P>::getParticles()
{
  CH_TIME("TracerParticleSolver::getParticles()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getParticles()" << endl;
  }

  return m_particles;
}

template <typename P>
inline const ParticleContainer<P>&
TracerParticleSolver<P>::getParticles() const
{
  CH_TIME("TracerParticleSolver::getParticles()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getParticles()" << endl;
  }

  return m_particles;
}

template <typename P>
inline const EBAMRCellData&
TracerParticleSolver<P>::getVelocityField() const
{
  CH_TIME("TracerParticleSolver::getVelocityField()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getVelocityField()" << endl;
  }

  CH_assert(!(m_velocityField == nullptr));

  return *m_velocityField;
}

template <typename P>
template <typename T, const Real& (T::*particleScalarFunction)() const>
void
TracerParticleSolver<P>::depositParticles(EBAMRCellData&             a_phi,
                                          const ParticleContainer<T> a_particles,
                                          const DepositionType       a_baseDeposition,
                                          const CoarseFineDeposition a_coarseFineDeposition) const
{
  CH_TIME("TracerParticleSolver::depositParticles()");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles()" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    m_amr->depositParticles<T, particleScalarFunction>(a_phi,
                                                       m_realm,
                                                       m_phase,
                                                       a_particles,
                                                       a_baseDeposition,
                                                       CoarseFineDeposition::Interp,
                                                       false);

    break;
  }
  case CoarseFineDeposition::Halo: {
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);
    a_particles.copyMaskParticles(mask);

    m_amr->depositParticles<T, particleScalarFunction>(a_phi,
                                                       m_realm,
                                                       m_phase,
                                                       a_particles,
                                                       a_baseDeposition,
                                                       CoarseFineDeposition::Halo,
                                                       false);

    a_particles.clearMaskParticles();

    break;
  }
  default: {
    MayDay::Error("TracerParticleSolver::depositParticles -- logic bust due to unsupported coarse-fine deposition");
  }
  }
}

#include <CD_NamespaceFooter.H>

#endif
