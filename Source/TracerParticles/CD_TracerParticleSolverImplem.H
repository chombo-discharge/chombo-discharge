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
#include <CD_NamespaceHeader.H>

template <size_t M, size_t N>
inline
TracerParticleSolver<M, N>::TracerParticleSolver() {
  CH_TIME("TracerParticleSolver::TracerParticleSolver()");
  
  // Default settings
  m_verbosity         = -1;
  m_name              = "TracerParticleSolver";
  m_className         = "TracerParticleSolver";
  m_realm             = Realm::primal;
  m_phase             = phase::gas;
  m_pvrBuffer         = 0;
  m_haloBuffer        = 1;
  m_plotVelocityField = true;
  m_deposition        = DepositionType::CIC;

  // Set default names for the particle scalars
  for (size_t i = 0; i < M; i++){
    m_scalarNames[i] = "particle_scalar #" + std::to_string(i);
  }

  // Set default names for the particle scalars
  for (size_t i = 0; i < N; i++){
    m_vectorNames[i] = "particle_vector #" + std::to_string(i);
  }
}

template <size_t M, size_t N>
inline
TracerParticleSolver<M, N>::TracerParticleSolver(const RefCountedPtr<AmrMesh>& a_amr,
						 const RefCountedPtr<ComputationalGeometry> a_compGeom) : TracerParticleSolver<M, N>() {
  CH_TIME("TracerParticleSolver::TracerParticleSolver(RefCountedPtr<AmrMesh>, RefCountedPtr<ComputationalGeometry>)");  

  m_amr                   = a_amr;
  m_computationalGeometry = a_compGeom;
}

template <size_t M, size_t N>
inline
TracerParticleSolver<M, N>::~TracerParticleSolver() {
  CH_TIME("TracerParticleSolver::~TracerParticleSolver()");
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::parseOptions() {
  CH_TIME("TracerParticleSolver::parseOptions()");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseOptions()" << endl;
  }

  this->parseDeposition();
  this->parsePlotVariables();
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::parseRuntimeOptions() {
  CH_TIME("TracerParticleSolver::parseRuntimeOptions()");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }

  this->parseDeposition();
  this->parsePlotVariables();  
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::parseDeposition() {
  CH_TIME("TracerParticleSolver::parseDeposition()");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseDeposition()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  // Deposition for particle-mesh operations
  pp.get("deposition", str);
  if(str == "ngp") {
    m_pvrBuffer  = 0;
    m_haloBuffer = 0;
    
    m_deposition           = DepositionType::NGP;
    m_coarseFineDeposition = CoarseFineDeposition::Halo;    
  }
  else if(str == "cic") {
    m_pvrBuffer  = 0;
    m_haloBuffer = 1;
    
    m_deposition           = DepositionType::CIC;
    m_coarseFineDeposition = CoarseFineDeposition::Halo;
  }
  else if(str == "tsc") {
    m_pvrBuffer  = 0;
    m_haloBuffer = 2;
    
    m_deposition           = DepositionType::TSC;
    m_coarseFineDeposition = CoarseFineDeposition::Halo;    
  }
  else{
    MayDay::Error("TracerParticleSolver::parseDeposition - unknown deposition method requested.");
  }  
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::parsePlotVariables() {
  CH_TIME("TracerParticleSolver::parsePlotVariables()");
  if(m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("plot_scalars", m_plotScalars);
  pp.get("plot_vectors", m_plotVectors);  
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::registerOperators() const {
  CH_TIME("TracerParticleSolver::registerOperators()");
  if(m_verbosity > 5) {
    pout() << m_name + "::registerOperators()" << endl;
  }

  if(m_amr.isNull()) {
    MayDay::Abort("TracerParticleSolver::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_particle_mesh,  m_realm, m_phase);

    // Register mask if using halo deposition
    if(m_haloBuffer > 0) {
      CH_assert(m_pvrBuffer == 0);
      
      m_amr->registerMask(s_particle_halo, m_haloBuffer, m_realm);
    }
  }
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::allocateInternals() {
  CH_TIME("TracerParticleSolver::allocateInternals()");
  if(m_verbosity > 5) {
    pout() << m_name + "::allocateInternals()" << endl;
  }

  // Allocate data for storing the velocity field on the mesh, and storage for the particles. 
  m_amr->allocate(m_particles, m_pvrBuffer, m_realm);
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::setRealm(const std::string& a_realm) {
  CH_TIME("TracerParticleSolver::setRealm(std::string)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setRealm(std::string)" << endl;
  }

  m_realm = a_realm;
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::setPhase(const phase::which_phase& a_phase) {
  CH_TIME("TracerParticleSolver::setPhase(phase::which_phase)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setPhase(phase::which_phase)" << endl;
  }

  m_phase = a_phase;
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh) {
  CH_TIME("TracerParticleSolver::setAmr(RefCountedPtr<AmrMesh>)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setAmr(RefCountedPtr<AmrMesh>)" << endl;
  }

  m_amr = a_amrMesh;
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_compGeom) {
  CH_TIME("TracerParticleSolver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  m_computationalGeometry = a_compGeom;
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::setVelocityField(const EBAMRCellData& a_velocityField) {
  CH_TIME("TracerParticleSolver::setVelocityField(EBAMRCellData)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setVelocityField(EBAMRCellData)" << endl;
  }

  m_velocityField = &a_velocityField;
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::preRegrid(const int a_lbase, const int a_oldFinestLevel) {
  CH_TIME("TracerParticleSolver::preRegrid(int, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::preRegrid(int, int)" << endl;
  }

  CH_assert(a_lbase >= 0);

  // Put particles in pre-regrid mode. 
  m_particles.preRegrid(a_lbase);
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("TracerParticleSolver::regrid(int, int, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }

  CH_assert(a_lmin           >= 0);
  CH_assert(a_oldFinestLevel >= 0);
  CH_assert(a_newFinestLevel >= 0);

  // Regrid the particles. 
  const Vector<DisjointBoxLayout>& grids   = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>&     domains = m_amr->getDomains();
  const Vector<Real>&              dx      = m_amr->getDx();
  const Vector<int>&               refRat  = m_amr->getRefinementRatios();

  m_particles.regrid(grids, domains, dx, refRat, a_lmin, a_newFinestLevel);
}

template <size_t M, size_t N>
inline
int TracerParticleSolver<M, N>::getNumberOfPlotVariables() const {
  CH_TIME("TracerParticleSolver::getNumberOfPlotVariables()");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  int numPlotVars = 0;

  if(m_plotVelocityField) numPlotVars += SpaceDim;
  if(m_plotScalars      ) numPlotVars += M;
  if(m_plotVectors      ) numPlotVars += N * SpaceDim;

  return numPlotVars;
}

template <size_t M, size_t N>
inline
Vector<std::string> TracerParticleSolver<M, N>::getPlotVariableNames() const {
  CH_TIME("TracerParticleSolver::getPlotVariableNames()");
  if(m_verbosity > 5){
    pout() << m_name + "::getPlotVariableNames()" << endl;
  }

  // TLDR: Possible plot variables is the density (m_phi), diffusion coefficient, source term, velocity, and eb flux. This
  //       function returns the associated plot variable names, and will be used in the plot files. 
  
  Vector<std::string> plotVarNames(0);
  
  if(m_plotVelocityField) plotVarNames.push_back("x-Velocity " + m_name);
  if(m_plotVelocityField) plotVarNames.push_back("y-Velocity " + m_name);
#if CH_SPACEDIM==3
  if(m_plotVelocityField) plotVarNames.push_back("z-Velocity " + m_name);
#endif

  if(m_plotScalars) {
    for (const auto& s : m_scalarNames) {
      plotVarNames.push_back(m_name + " " + s);
    }
  }

  if(m_plotVectors) {
    for (const auto& v : m_vectorNames) {
      plotVarNames.push_back("x-" + v + " " + m_name);
      plotVarNames.push_back("y-" + v + " " + m_name);      
#if CH_SPACEDIM==3
      plotVarNames.push_back("z-" + v + " " + m_name);
#endif
    }
  }  
  
  return plotVarNames;
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::setScalarNames(const std::array<std::string, M>& a_scalarNames) {
  CH_TIME("TracerParticleSolver::setScalarNames(std::array)");
  if(m_verbosity > 5){
    pout() << m_name + "::setScalarNames(std::array)" << endl;
  }

  m_scalarNames = a_scalarNames;
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::setVectorNames(const std::array<std::string, N>& a_vectorNames) {
  CH_TIME("TracerParticleSolver::setVectorNames(std::array)");
  if(m_verbosity > 5){
    pout() << m_name + "::setVectorNames(std::array)" << endl;
  }

  m_vectorNames = a_vectorNames;  
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::writePlotFile() {
  CH_TIME("TracerParticleSolver::writePlotFile()");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotFile()" << endl;
  }

  MayDay::Error("TracerParticleSolver::writePlotFile -- not implemented (yet)");  
}

template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::writePlotData(EBAMRCellData& a_output, int& a_icomp) {
  CH_TIME("TracerParticleSolver::writePlotData(EBAMRCellData, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotFile(EBAMRCellData, int)" << endl;
  }

  // Write the velocity field to the output data holder.
  if(m_plotVelocityField) {
    const Interval srcInterval(0, 0);
    const Interval dstInterval(a_icomp, a_icomp + SpaceDim - 1);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      (*m_velocityField)[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
    }
    
    a_icomp += SpaceDim;
  }

  // Write the particle scalars.
  if(m_plotScalars) {
    EBAMRCellData scratch;
    m_amr->allocate(scratch, m_realm, m_phase, 1);

    // Deposit everything onto scratch and copy to the output data holder.
    for (size_t i = 0; i < M; i++){
      DataOps::setValue(scratch, 0.0);
      this->depositScalar<TracerParticle<M, N>, &TracerParticle<M,N>::template getScalar<0> >(scratch, m_particles, m_deposition, m_coarseFineDeposition);
    }
  }

  MayDay::Error("TracerParticleSolver::writePlotData -- not implemented (yet)");
}

#ifdef CH_USE_HDF5
template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("TracerParticleSolver::writeCheckpointLevel(HDF5Handle, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  MayDay::Error("TracerParticleSolver::writeCheckpointLevel -- not implemented (yet)");  
}
#endif

#ifdef CH_USE_HDF5
template <size_t M, size_t N>
inline
void TracerParticleSolver<M, N>::readCheckpointLevel(HDF5Handle& a_handle, const int a_level) {
  CH_TIME("TracerParticleSolver::readCheckpointLevel(HDF5Handle, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::readCheckpointLevel(HDF5Handle, int)" << endl;
  }

  MayDay::Error("TracerParticleSolver::readCheckpointLevel -- not implemented (yet)");  
}
#endif

template <size_t M, size_t N>
inline
ParticleContainer<TracerParticle<M, N> >& TracerParticleSolver<M, N>::getParticles() {
  CH_TIME("TracerParticleSolver::getParticles()");
  if(m_verbosity > 5){
    pout() << m_name + "::getParticles()" << endl;
  }

  return m_particles;
}

template <size_t M, size_t N>
inline
const ParticleContainer<TracerParticle<M, N> >& TracerParticleSolver<M, N>::getParticles() const {
  CH_TIME("TracerParticleSolver::getParticles()");
  if(m_verbosity > 5){
    pout() << m_name + "::getParticles()" << endl;
  }

  return m_particles;
}

template <size_t M, size_t N>
inline
const EBAMRCellData& TracerParticleSolver<M, N>::getVelocityField() const {
  CH_TIME("TracerParticleSolver::getVelocityField()");
  if(m_verbosity > 5){
    pout() << m_name + "::getVelocityField()" << endl;
  }

  return *m_velocityField;
}

#include <CD_NamespaceFooter.H>

#endif
