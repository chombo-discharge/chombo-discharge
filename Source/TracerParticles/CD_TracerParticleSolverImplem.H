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
  m_verbosity            = -1;
  m_name                 = "TracerParticleSolver";
  m_className            = "TracerParticleSolver";
  m_realm                = Realm::primal;
  m_phase                = phase::gas;
  m_haloBuffer           = 1;
  m_plotVelocity         = true;
  m_plotWeight           = true;
  m_volumeScale          = false;
  m_deposition           = DepositionType::CIC;
  m_interpolation        = DepositionType::CIC;
  m_coarseFineDeposition = CoarseFineDeposition::Halo;
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
    m_deposition = DepositionType::NGP;
  }
  else if (str == "cic") {
    m_deposition = DepositionType::CIC;
  }
  else {
    MayDay::Error("TracerParticleSolver::parseDeposition - unknown deposition method requested.");
  }

  pp.get("deposition_cf", str);
  if (str == "interp") {
    m_coarseFineDeposition = CoarseFineDeposition::Interp;
  }
  else if (str == "halo") {
    m_coarseFineDeposition = CoarseFineDeposition::Halo;
  }
  else if (str == "halo_ngp") {
    m_coarseFineDeposition = CoarseFineDeposition::HaloNGP;
  }
  else {
    MayDay::Error("TracerParticleSolver::parseDeposition - unknown coarse-fine deposition method requested.");
  }

  // Interpolation type particle-mesh operations
  pp.get("interpolation", str);
  if (str == "ngp") {
    m_interpolation = DepositionType::NGP;
  }
  else if (str == "cic") {
    m_interpolation = DepositionType::CIC;
  }
  else {
    MayDay::Error("TracerParticleSolver::parseDeposition - unknown interpolation method requested.");
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

  pp.get("plot_weight", m_plotWeight);
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

    // For CIC halo deposition.
    m_amr->registerMask(s_particle_halo, 1, m_realm);
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
  m_amr->allocate(m_velocityField, m_realm, m_phase, SpaceDim);
}

template <typename P>
inline void
TracerParticleSolver<P>::setName(const std::string& a_name) noexcept
{
  CH_TIME("TracerParticleSolver::setName(std::string)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setName(std::string)" << endl;
  }

  m_name = a_name;
}

template <typename P>
inline void
TracerParticleSolver<P>::setVolumeScale(const bool a_scale) noexcept
{
  CH_TIME("TracerParticleSolver::setVolumeScale(bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVolumeScale(bool)" << endl;
  }

  m_volumeScale = a_scale;
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

  DataOps::copy(m_velocityField, a_velocityField);
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

  m_amr->remapToNewGrids(m_particles, a_lmin, a_newFinestLevel);
  m_amr->allocate(m_velocityField, m_realm, m_phase, SpaceDim);
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
TracerParticleSolver<P>::setDeposition(const DepositionType a_deposition) noexcept
{
  CH_TIME("TracerParticleSolver::setDeposition");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDeposition" << endl;
  }

  m_deposition = a_deposition;
}

template <typename P>
inline void
TracerParticleSolver<P>::deposit(EBAMRCellData& a_phi) const noexcept
{
  CH_TIME("TracerParticleSolver::deposit(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::deposit(EBAMRCellData)" << endl;
  }

  DataOps::setValue(a_phi, 0.0);
  this->depositParticles<P, &P::weight>(a_phi, m_particles, m_deposition, m_coarseFineDeposition);

  if (m_volumeScale) {
    DataOps::volumeScale(a_phi, m_amr->getDx());
  }
}

template <typename P>
void
TracerParticleSolver<P>::interpolateWeight(const EBAMRCellData& a_scalar) noexcept
{
  CH_TIME("TracerParticleSolver::interpolateWeight(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateWeight(EBAMRCellData)" << endl;
  }

  m_amr->interpolateParticles<P, &P::weight>(m_particles, m_realm, m_phase, a_scalar, m_interpolation, true);
}

template <typename P>
inline void
TracerParticleSolver<P>::interpolateVelocities()
{
  CH_TIME("TracerParticleSolver::interpolateVelocities()");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateVelocities()" << endl;
  }

  m_amr->interpolateParticles<P, &P::velocity>(m_particles, m_realm, m_phase, m_velocityField, m_interpolation, true);
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
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->writePlotData(*output[lvl], icomp, m_realm, lvl);
  }

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
inline int
TracerParticleSolver<P>::getNumberOfPlotVariables() const
{
  CH_TIME("TracerParticleSolver::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  int numPlotVars = 0;

  if (m_plotWeight) {
    numPlotVars += 1;
  }
  if (m_plotVelocity) {
    numPlotVars += SpaceDim;
  }

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

  Vector<std::string> plotVarNames(0);

  if (m_plotWeight) {
    plotVarNames.push_back(m_name + " density");
  }
  if (m_plotVelocity) {
    plotVarNames.push_back("x-velocity " + m_name);
  }
  if (m_plotVelocity) {
    plotVarNames.push_back("y-velocity " + m_name);
  }
#if CH_SPACEDIM == 3
  if (m_plotVelocity) {
    plotVarNames.push_back("z-velocity " + m_name);
  }
#endif

  return plotVarNames;
}

template <typename P>
inline void
TracerParticleSolver<P>::writePlotData(LevelData<EBCellFAB>& a_output,
                                       int&                  a_comp,
                                       const std::string     a_outputRealm,
                                       const int             a_level) const noexcept
{
  CH_TIME("TracerParticleSolver::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  // Write the particle scalars.
  if (m_plotWeight) {
    EBAMRCellData weight;
    m_amr->allocate(weight, m_realm, m_phase, 1);

    this->deposit(weight);

    this->writeData(a_output, a_comp, weight, a_outputRealm, a_level, false, true);
  }

  // Write the velocity field to the output data holder.
  if (m_plotVelocity) {
    this->writeData(a_output, a_comp, m_velocityField, a_outputRealm, a_level, true, true);
  }
}

template <typename P>
void
TracerParticleSolver<P>::writeData(LevelData<EBCellFAB>& a_output,
                                   int&                  a_comp,
                                   const EBAMRCellData&  a_data,
                                   const std::string     a_outputRealm,
                                   const int             a_level,
                                   const bool            a_interpToCentroids,
                                   const bool            a_interpGhost) const noexcept
{
  CH_TIMERS("TracerParticleSolver::writeData");
  CH_TIMER("TracerParticleSolver::writeData::allocate", t1);
  CH_TIMER("TracerParticleSolver::writeData::local_copy", t2);
  CH_TIMER("TracerParticleSolver::writeData::interp_ghost", t3);
  CH_TIMER("TracerParticleSolver::writeData::interp_centroid", t4);
  CH_TIMER("TracerParticleSolver::writeData::final_copy", t5);
  if (m_verbosity > 5) {
    pout() << m_name + "::writeData" << endl;
  }

  // Number of components we are working with.
  const int numComp = a_data[a_level]->nComp();

  CH_START(t1);
  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, m_phase, a_level, numComp);
  CH_STOP(t1);

  CH_START(t2);
  m_amr->copyData(scratch, *a_data[a_level], a_level, m_realm, m_realm);
  CH_START(t2);

  // Interpolate ghost cells
  CH_START(t3);
  if (a_level > 0 && a_interpGhost) {
    m_amr->interpGhost(scratch, *a_data[a_level - 1], a_level, m_realm, m_phase);
  }
  CH_STOP(t3);

  CH_START(t4);
  if (a_interpToCentroids) {
    m_amr->interpToCentroids(scratch, m_realm, m_phase, a_level);
  }
  CH_STOP(t4);

  DataOps::setCoveredValue(scratch, 0.0);

  CH_START(t5);
  const Interval srcInterv(0, numComp - 1);
  const Interval dstInterv(a_comp, a_comp + numComp - 1);
  m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, dstInterv, srcInterv);
  CH_STOP(t5);

  a_comp += numComp;
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
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = m_amr->getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : dt)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& particles = m_particles[lvl][din].listItems();

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

  return m_velocityField;
}

template <typename P>
inline DepositionType
TracerParticleSolver<P>::getDepositionType() const
{
  CH_TIME("TracerParticleSolver::getDepositionType");
  if (m_verbosity > 5) {
    pout() << m_name + "::getDepositionType" << endl;
  }

  return m_deposition;
}

template <typename P>
inline CoarseFineDeposition
TracerParticleSolver<P>::getCoarseFineDepositionType() const
{
  CH_TIME("TracerParticleSolver::getCoarseFineDepositionType");
  if (m_verbosity > 5) {
    pout() << m_name + "::getCoarseFineDepositionType" << endl;
  }

  return m_coarseFineDeposition;
}

template <typename P>
inline DepositionType
TracerParticleSolver<P>::getInterpolationType() const
{
  CH_TIME("TracerParticleSolver::getInterpolationType");
  if (m_verbosity > 5) {
    pout() << m_name + "::getInterpolationType" << endl;
  }

  return m_interpolation;
}

template <typename P>
template <typename T, const Real& (T::*particleScalarFunction)() const>
void
TracerParticleSolver<P>::depositParticles(EBAMRCellData&              a_phi,
                                          const ParticleContainer<T>& a_particles,
                                          const DepositionType        a_baseDeposition,
                                          const CoarseFineDeposition  a_coarseFineDeposition) const noexcept
{
  CH_TIME("TracerParticleSolver::depositParticles()");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles()" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

  // I HATE it when other people do this. Yet I keep doing it.
  ParticleContainer<T>& mutableParticles = (ParticleContainer<T>&)a_particles;

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

    // Copy particles living on the mask.
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);

    mutableParticles.copyMaskParticles(mask);

    m_amr->depositParticles<T, particleScalarFunction>(a_phi,
                                                       m_realm,
                                                       m_phase,
                                                       a_particles,
                                                       a_baseDeposition,
                                                       CoarseFineDeposition::Halo,
                                                       false);

    // Clear out the mask particles.
    mutableParticles.clearMaskParticles();

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);

    // Transfer particles living on the mask.
    mutableParticles.transferMaskParticles(mask);

    m_amr->depositParticles<T, particleScalarFunction>(a_phi,
                                                       m_realm,
                                                       m_phase,
                                                       a_particles,
                                                       a_baseDeposition,
                                                       CoarseFineDeposition::HaloNGP,
                                                       false);

    // Transfer them back.
    mutableParticles.transferParticles(mutableParticles.getMaskParticles());

    break;
  }
  default: {
    MayDay::Error("TracerParticleSolver::depositParticles -- logic bust due to unsupported coarse-fine deposition");
  }
  }
}

#include <CD_NamespaceFooter.H>

#endif
