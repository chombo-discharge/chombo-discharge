/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaStepper.cpp
  @brief  Implementation of CD_CdrPlasmaStepper.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBLevelDataOps.H>
#include <EBArith.H>
#include <PolyGeom.H>

// Our includes
#include <CD_CdrPlasmaStepper.H>
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include <CD_Units.H>
#include <CD_Timer.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

// C-style VLA methods can have occasional memory issues. Use them at your own peril. 
#define USE_FAST_REACTIONS  0
#define USE_FAST_VELOCITIES 0
#define USE_FAST_DIFFUSION  0

using namespace Physics::CdrPlasma;

CdrPlasmaStepper::CdrPlasmaStepper(){
  CH_TIME("CdrPlasmaStepper::CdrPlasmaStepper()");
  
  // Default settings
  m_className       = "CdrPlasmaStepper";
  m_verbosity       = -1;
  m_solverVerbosity = -1;
  m_phase           = phase::gas;
  m_realm           = Realm::Primal;
}

CdrPlasmaStepper::CdrPlasmaStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics) : CdrPlasmaStepper() {
  CH_TIME("CdrPlasmaStepper::CdrPlasmaStepper(RefCountedPtr<CdrPlasmaPhysics>)");
  
  m_physics = a_physics;
}

void CdrPlasmaStepper::postInitialize() {
  CH_TIME("CdrPlasmaStepper::postInitialize");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::postInitialize" << endl;
  }
}

void CdrPlasmaStepper::registerRealms(){
  CH_TIME("CdrPlasmaStepper::registerRealms");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::registerRealms" << endl;
  }
  
  m_amr->registerRealm(m_realm);
}

void CdrPlasmaStepper::postRegrid(){

}

void CdrPlasmaStepper::registerOperators(){
  CH_TIME("CdrPlasmaStepper::registerOperators()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::registerOperators()" << endl;
  }

  // Solvers know what they need so they can just register the operators themselves. 
  m_cdr        ->registerOperators();
  m_fieldSolver->registerOperators();
  m_rte        ->registerOperators();
  m_sigma      ->registerOperators();
}

CdrPlasmaStepper::~CdrPlasmaStepper(){
  CH_TIME("CdrPlasmaStepper::~CdrPlasmaStepper()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::~CdrPlasmaStepper()" << endl;
  }  
}

bool CdrPlasmaStepper::stationaryRTE(){
  CH_TIME("CdrPlasmaStepper::stationaryRTE()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::stationaryRTE()" << endl;
  }

  return m_rte->isStationary();
}

void CdrPlasmaStepper::computeSpaceChargeDensity(){
  CH_TIME("CdrPlasmaStepper::computeSpaceChargeDensity()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeSpaceChargeDensity()" << endl;
  }

  this->computeSpaceChargeDensity(m_fieldSolver->getRho(), m_cdr->getPhis());
}

void CdrPlasmaStepper::computeSpaceChargeDensity(MFAMRCellData&                 a_rho,
						 const Vector<EBAMRCellData*>&  a_densities) {
  CH_TIME("CdrPlasmaStepper::computeSpaceChargeDensity(MFAMRCellData, Vector(EBAMRCellData))");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeSpaceChargeDensity(MFAMRCellData, Vector(EBAMRCellData))" << endl;
  }

  // Set space charge to zero everywhere. We will then compute the space charge on
  // the gas phase (solid phase space charge is always zero in the CdrPlasma module). 
  DataOps::setValue(a_rho, 0.0);

  // Get the data holder for the gas phase (the MFAMRCellData holds data on both phases)
  EBAMRCellData rhoGas = m_amr->alias(phase::gas, a_rho);


  // Go through the CDR solvers and add their space charge density. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const EBAMRCellData&             density = *(a_densities[solverIt.index()]);
    const RefCountedPtr<CdrSpecies>& spec    = solverIt.getSpecies();
    const int                        Z       = spec->getChargeNumber();
      
    if(Z != 0){
      DataOps::incr(rhoGas, density, Real(Z));
    }
  }

  // Scale by electron charge. 
  DataOps::scale(a_rho, Units::Qe);

  // Above, we computed the cell-centered space charge. We must have centroid-centered. 
  m_amr->averageDown      (a_rho,  m_realm            );
  m_amr->interpGhostMG    (a_rho,  m_realm            );
  m_amr->interpToCentroids(rhoGas, m_realm, phase::gas);
}

void CdrPlasmaStepper::computeCellConductivity(EBAMRCellData& a_cellConductivity) const {
  CH_TIME("CdrPlasmaStepper::computeCellConductivity(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCellConductivity(EBAMRCellData)" << endl;
  }

  // Get a handle to the 
  const EBAMRCellData cellCenteredElectricField = m_amr->alias(m_phase, m_fieldSolver->getElectricField());

  // Allocate some storage for holding the electric field magnitude and one species conductivity.
  EBAMRCellData fieldMagnitude;      
  EBAMRCellData speciesConductivity; 

  m_amr->allocate(fieldMagnitude,      m_realm, m_phase, 1);
  m_amr->allocate(speciesConductivity, m_realm, m_phase, 1);

  // Compute the electric field magnitude
  DataOps::vectorLength(fieldMagnitude, cellCenteredElectricField);

  // Reset the conductivity.
  DataOps::setValue(a_cellConductivity, 0.0);

  // Increment the conductivity by the drift contribution from each CDR species. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver> & solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solverIt.getSpecies();

    // Compute the conductivity for this species. This is just Z*mu*n.
    if(solver->isMobile()){
      const EBAMRCellData& cellVel = solver->getCellCenteredVelocity();
      const EBAMRCellData& phi     = solver->getPhi();

      const int Z = species->getChargeNumber();

      // In the below comments, f = speciesConductivity. The loop just adds to the total conductivity. 
      if(Z != 0){
	DataOps::vectorLength  (speciesConductivity, cellVel       ); // Compute f = |v|
	DataOps::divideByScalar(speciesConductivity, fieldMagnitude); // Compute f = |v|/|E| = |mu|
	DataOps::multiply      (speciesConductivity, phi           ); // Compute f = mu*n

	// Add to total conductivity. 
	DataOps::incr(a_cellConductivity, speciesConductivity, 1.0*std::abs(Z));
      }
    }
  }

  // Need to scale by electron charge.
  DataOps::scale(a_cellConductivity, Units::Qe); 
}

void CdrPlasmaStepper::computeFaceConductivity(EBAMRFluxData&       a_conductivityFace,
					       EBAMRIVData&         a_conductivityEB,
					       const EBAMRCellData& a_conductivityCell) const {
  CH_TIME("CdrPlasmaStepper::computeFaceConductivity(EBAMRFluxData, EBAMRIVData, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeFaceConductivity(EBAMRFluxData, EBAMRIVData, EBAMRCellData)" << endl;
  }

  // Init to something stupid in order to catch errors. 
  DataOps::setValue(a_conductivityFace, std::numeric_limits<Real>::max());
  DataOps::setValue(a_conductivityEB,   std::numeric_limits<Real>::max());

#if 1
  // Average the cell-centered conductivity to faces. Note that this includes one "ghost face", which we need
  // because the multigrid solver will interpolate face-centered conductivities to face centroids. 
  DataOps::averageCellScalarToFaceScalar(a_conductivityFace, a_conductivityCell, m_amr->getDomains());
#else
  DataOps::averageCellToFace(a_conductivityFace, a_conductivityCell, m_amr->getDomains());
#endif


#if 1
  // Now compute the conductivity on the EB. 
  const auto& interpStencil = m_amr->getCentroidInterpolationStencils(m_realm, m_phase);
  interpStencil.apply(a_conductivityEB, a_conductivityCell);
#else
  DataOps::setValue(a_conductivityEB, 0.0);
  DataOps::incr(a_conductivityEB, a_conductivityCell, 1.0);
#endif
}

void CdrPlasmaStepper::setupSemiImplicitPoisson(const Real a_dt){
  CH_TIME("CdrPlasmaStepper::setupSemiImplicitPoisson(Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setupSemiImplicitPoisson(Real)" << endl;
  }

  // Storage that is needed. 
  EBAMRCellData conductivityCell;
  EBAMRFluxData conductivityFace;
  EBAMRIVData   conductivityEB;

  m_amr->allocate(conductivityCell, m_realm, phase::gas, 1);
  m_amr->allocate(conductivityFace, m_realm, phase::gas, 1);
  m_amr->allocate(conductivityEB,   m_realm, phase::gas, 1);

  // Compute the cell-centered conductivity first. Then compute the face- and EB-centered conductivities by averaging the
  // cell-centered conductivity. 
  this->computeCellConductivity(conductivityCell);
  this->computeFaceConductivity(conductivityFace, conductivityEB, conductivityCell);

  // Set up the semi implicit Poisson equation. Note that the Poisson equation will have modified face-centered permittivities.
  this->setupSemiImplicitPoisson(conductivityFace, conductivityEB, a_dt/Units::eps0);
}

void CdrPlasmaStepper::setupSemiImplicitPoisson(const EBAMRFluxData& a_conductivityFace,
						const EBAMRIVData&   a_conductivityEB,
						const Real           a_factor) {
  CH_TIME("CdrPlasmaStepper::setupSemiImplicitPoisson(EBAMRFluxData, EBAMRIVData, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setupSemiImplicitPoisson(EBAMRFluxData, EBAMRIVData, Real)" << endl;
  }

  // TLDR: This routine will set up the field solver from an equation
  //
  // div(B*grad(phi)) where B = epsr + a_factor*sigma and sigma = conductivity.

  // First, the field solver must set the permittivities as usual. The "permittivities" that we are after are
  // eps = epsr + dt*sigma/eps0 (but we only have plasma on the gas phase).
  m_fieldSolver->setPermittivities();

  // Get the permittivities on the faces.
  MFAMRFluxData& permFace = m_fieldSolver->getPermittivityFace();
  MFAMRIVData&   permEB   = m_fieldSolver->getPermittivityEB();

  // Get handles to the gas-phase permittivities. 
  EBAMRFluxData permFaceGas = m_amr->alias(phase::gas, permFace);
  EBAMRIVData   permEBGas   = m_amr->alias(phase::gas, permEB  );

  // Increment the field solver permittivities by a_factor*sigma. After this, the "permittivities" are
  // given by epsr + a_factor*sigma
  DataOps::incr(permFaceGas, a_conductivityFace, a_factor);
  DataOps::incr(permEBGas,   a_conductivityEB,   a_factor);

  // Set up the solver with the new "permittivities".
  m_fieldSolver->setupSolver();
}

bool CdrPlasmaStepper::solvePoisson(){
  CH_TIME("CdrPlasmaStepper::solvePoisson()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::solvePoisson()" << endl;
  }

  return this->solvePoisson(m_fieldSolver->getPotential(),
			    m_fieldSolver->getRho      (),
			    m_cdr        ->getPhis     (),
			    m_sigma      ->getPhi      ());
}

bool CdrPlasmaStepper::solvePoisson(MFAMRCellData&                a_potential,
				    MFAMRCellData&                a_rho,
				    const Vector<EBAMRCellData*>  a_cdrDensities,
				    const EBAMRIVData&            a_sigma){
  CH_TIME("CdrPlasmaStepper::solvePoisson(MFAMRCellData, MFAMRCellData, Vector<EBAMRCellData*>, EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::solvePoisson(MFAMRCellData, MFAMRCellData, Vector<EBAMRCellData*>, EBAMRIVData)" << endl;
  }

  // Compute the space charge density onto the input data holder.
  this->computeSpaceChargeDensity(a_rho, a_cdrDensities);

  // Field solver solves for the potential.
  const bool converged = m_fieldSolver->solve(a_potential, a_rho, a_sigma, false);

  // Return whether or not the field solve converged. 
  return converged;
}

void CdrPlasmaStepper::allocateInternals(){
  CH_TIME("CdrPlasmaStepper::allocateInternals"); 
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::allocateInternals" << endl;
  }

  /// Do nothing
}

void CdrPlasmaStepper::deallocateInternals(){
  CH_TIME("CdrPlasmaStepper::deallocateInternals"); 
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::deallocateInternals" << endl;
  }

  /// Do nothing
}

void CdrPlasmaStepper::advanceReactionNetwork(const Real a_time, const Real a_dt){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetwork(solvers)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetwork(solvers)" << endl;
  }

  // Compute the electric field.
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, m_cdr->getPhase(), SpaceDim);
  this->computeElectricField(E, m_cdr->getPhase(), m_fieldSolver->getPotential());

  // Get the source terms and densities from the solvers. 
  Vector<EBAMRCellData*> cdrSources   = m_cdr->getSources();
  Vector<EBAMRCellData*> rteSources   = m_rte->getSources();
  Vector<EBAMRCellData*> cdrDensities = m_cdr->getPhis();
  Vector<EBAMRCellData*> rteDensities = m_rte->getPhis();

  // Call the other version (without the gradient)
  this->advanceReactionNetwork(cdrSources,
			       rteSources,
			       cdrDensities,
			       rteDensities,
			       E,
			       a_time,
			       a_dt);
}

void CdrPlasmaStepper::advanceReactionNetwork(Vector<EBAMRCellData*>&       a_cdrSources,
					      Vector<EBAMRCellData*>&       a_rteSources,
					      const Vector<EBAMRCellData*>& a_cdrDensities,
					      const Vector<EBAMRCellData*>& a_rteDensities,
					      const EBAMRCellData&          a_E,
					      const Real&                   a_time,
					      const Real&                   a_dt){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetwork(Vector<EBAMRCellData>x4, EBAMRCellData, Real, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetwork(Vector<EBAMRCellData>x4, EBAMRCellData, Real, Real)" << endl;
  }

  // TLDR: This version of advanceReactionNetwork will compute the gradients of the plasma species densities and then
  //       call the general version. The gradients are needed for various purposes.

  // Allocate scratch data which we will use to compute the gradient. 
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_cdr->getPhase(), 1);

  // Number of CDR solvers. 
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  // Compute the cell-centered gradient for each plasma species. 
  Vector<EBAMRCellData*> cdrGradients(numCdrSpecies); 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();

    // Create some storage where we store the gradient. Note that we use the new operator so this
    // storage must be deleted later. 
    cdrGradients[idx] = new EBAMRCellData();
    m_amr->allocate(*cdrGradients[idx], m_realm, m_cdr->getPhase(), SpaceDim);  

    // Copy the cell-centered density to the scratch data holder. Need to do this because the CDR solvers may not
    // have updated their ghost cells. 
    DataOps::copy(scratch, *a_cdrDensities[idx]);
    m_amr->interpGhostMG(scratch, m_realm, m_cdr->getPhase());

    // Compute the gradient and coarsen the result. 
    m_amr->computeGradient(*cdrGradients[idx], scratch, m_realm, m_cdr->getPhase());
    m_amr->averageDown    (*cdrGradients[idx],          m_realm, m_cdr->getPhase());
    m_amr->interpGhost    (*cdrGradients[idx],          m_realm, m_cdr->getPhase());
  }

  // Call the other version.
  this->advanceReactionNetwork(a_cdrSources,
			       a_rteSources,
			       a_cdrDensities,
			       cdrGradients,
			       a_rteDensities,
			       a_E,
			       a_time,
			       a_dt);

  // Delete the extra storage since we didn't use smart pointers. 
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->deallocate(*cdrGradients[idx]);
    delete cdrGradients[idx];
  }
}

void CdrPlasmaStepper::advanceReactionNetwork(Vector<EBAMRCellData*>&       a_cdrSources,
					      Vector<EBAMRCellData*>&       a_rteSources,
					      const Vector<EBAMRCellData*>& a_cdrDensities,
					      const Vector<EBAMRCellData*>& a_cdrGradients,
					      const Vector<EBAMRCellData*>& a_rteDensities,
					      const EBAMRCellData&          a_E,
					      const Real&                   a_time,
					      const Real&                   a_dt){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetwork(Vector<EBAMRCellData>x5, EBAMRCellData, Real, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetwork(Vector<EBAMRCellData>x5, EBAMRCellData, Real, Real)" << endl;
  }

  // TLDR: This is the general version for computing the CDR and RTE source terms. This will call the level version. 

  const int numCdrSpecies = m_physics->getNumCdrSpecies();
  const int numRteSpecies = m_physics->getNumRtSpecies();


  // This is a special option in case we use upwinding. In that case we need to allocate
  // storage for holding the upwind-weighted data for each cdr solver.
  Vector<EBAMRCellData*> cdrStates(numCdrSpecies);
  
  if(m_whichSourceTermComputation == SourceTermComputation::Upwind){
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // Allocate some storage and compute an upwind approximation to the cell-centered density. This is the Villa et. al magic for stabilizing
      // the drift-reaction mechanism. Only use if you absolute know what you're doing. 
      cdrStates[idx] = new EBAMRCellData();
      m_amr->allocate(*cdrStates[idx], m_realm, m_cdr->getPhase(), 1);

      // Compute the upwind approximation. 
      solverIt()->weightedUpwind(*cdrStates[idx], m_upwindFactor);
    }
  }
  else{
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      cdrStates[idx] = a_cdrDensities[idx];
    }
  }


  // Level loop. We expose the entire reaction network method as a level version. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    // Things to be exposed by level. 
    Vector<LevelData<EBCellFAB>* > cdrSources  (numCdrSpecies);
    Vector<LevelData<EBCellFAB>* > rteSources  (numRteSpecies);    
    Vector<LevelData<EBCellFAB>* > cdrDensities(numCdrSpecies);
    Vector<LevelData<EBCellFAB>* > cdrGradients(numCdrSpecies);
    Vector<LevelData<EBCellFAB>* > rteDensities(numRteSpecies);

    // Get the CDR stuff on this level. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx     = solverIt.index();
      
      cdrSources  [idx] = (*a_cdrSources  [idx])[lvl];
      cdrDensities[idx] = (*cdrStates  [idx])[lvl];
      cdrGradients[idx] = (*a_cdrGradients[idx])[lvl];
    }

    // Get the RTE stuff on this level.     
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx     = solverIt.index();
      
      rteSources  [idx] = (*a_rteSources  [idx])[lvl];
      rteDensities[idx] = (*a_rteDensities[idx])[lvl];
    }

    // Now call the level version. 
    this->advanceReactionNetwork(cdrSources,
				 rteSources,
				 cdrDensities,
				 cdrGradients,
				 rteDensities,
				 *a_E[lvl],
				 a_time,
				 a_dt,
				 lvl);
  }

  // When we did the upwind approximation we explicitly allocate memory. Now release it. 
  if(m_whichSourceTermComputation == SourceTermComputation::Upwind){
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      delete cdrStates[idx];
    }
  }  
}

void CdrPlasmaStepper::advanceReactionNetwork(Vector<LevelData<EBCellFAB>* >&       a_cdrSources,
					      Vector<LevelData<EBCellFAB>* >&       a_rteSources,
					      const Vector<LevelData<EBCellFAB>* >& a_cdrDensities,
					      const Vector<LevelData<EBCellFAB>* >& a_cdrGradients,
					      const Vector<LevelData<EBCellFAB> *>& a_rteDensities,
					      const LevelData<EBCellFAB>&           a_E,
					      const Real&                           a_time,
					      const Real&                           a_dt,
					      const int                             a_lvl){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetwork(Vector<LD<EBCellFAB>x5, LD<EBCellFAB>, Real, Real, int)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetwork(Vector<LD<EBCellFAB>x5, LD<EBCellFAB>, Real, Real, int)" << endl;
  }

  const Real zero = 0.0;

  const int num_Photons  = m_physics->getNumRtSpecies();
  const int num_species  = m_physics->getNumCdrSpecies();


  // Stencils for extrapolating things to cell centroids
  const IrregAmrStencil<CentroidInterpolationStencil>& interp_stencils = m_amr->getCentroidInterpolationStencils(m_realm, m_cdr->getPhase());

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real dx                = m_amr->getDx()[a_lvl];
  const RealVect origin        = m_amr->getProbLo();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    Vector<EBCellFAB*> particle_sources(num_species);
    Vector<EBCellFAB*> cdrDensities(num_species);
    Vector<EBCellFAB*> cdrGradients(num_species);
    Vector<EBCellFAB*> particle_velocities(num_species, NULL);
    Vector<EBCellFAB*> rteSources(num_Photons);
    Vector<EBCellFAB*> rteDensities(num_Photons);

    
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<CdrSolver>& solver = solver_it();
      const int idx = solver_it.index();
      particle_sources[idx]   = &(*a_cdrSources[idx])[dit()];
      cdrDensities[idx] = &(*a_cdrDensities[idx])[dit()];
      cdrGradients[idx] = &(*a_cdrGradients[idx])[dit()];

      if(solver->isMobile()){
	const EBAMRCellData& velo = solver->getCellCenteredVelocity();
	particle_velocities[idx] = &((*velo[a_lvl])[dit()]);
      }
    }

    
    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      rteSources[idx]   = &(*a_rteSources[idx])[dit()];
      rteDensities[idx] = &(*a_rteDensities[idx])[dit()];
    }
    
    // This does all cells
#if USE_FAST_REACTIONS
    advanceReactionNetworkRegularCellsFast(particle_sources,
					   rteSources,
					   cdrDensities,
					   cdrGradients,
					   rteDensities,
					   a_E[dit()],
					   a_time,
					   a_dt,
					   dx,
					   dbl.get(dit()));

#else
    advanceReactionNetworkRegularCells(particle_sources,
				       rteSources,
				       cdrDensities,
				       cdrGradients,
				       rteDensities,
				       a_E[dit()],
				       a_time,
				       a_dt,
				       dx,
				       dbl.get(dit()));
#endif

    // This overwrites irregular celles
    advanceReactionNetworkIrreg(particle_sources,
				rteSources,
				cdrDensities,
				cdrGradients,
				particle_velocities,
				rteDensities,
				interp_stencils[a_lvl][dit()],
				a_E[dit()],
				a_time,
				a_dt,
				dx,
				dbl.get(dit()),
				a_lvl,
				dit());
  }
}

void CdrPlasmaStepper::advanceReactionNetworkRegularCells(Vector<EBCellFAB*>&       a_cdrSources,
							  Vector<EBCellFAB*>&       a_rteSources,
							  const Vector<EBCellFAB*>& a_cdrDensities,
							  const Vector<EBCellFAB*>& a_cdrGradients,
							  const Vector<EBCellFAB*>& a_rteDensities,
							  const EBCellFAB&          a_E,
							  const Real&               a_time,
							  const Real&               a_dt,
							  const Real&               a_dx,
							  const Box&                a_box){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkRegularCells(patch)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkRegularCells(patch)" << endl;
  }

  const Real zero = 0.0;
    

  const int num_species  = m_physics->getNumCdrSpecies();
  const int num_Photons  = m_physics->getNumRtSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     cdrDensities(num_species);
  Vector<RealVect> cdrGradients(num_species);
  Vector<Real>     rteSources(num_Photons);
  Vector<Real>     rteDensities(num_Photons);

  // Computed source terms onto here
  EBCellFAB part_src(ebisbox, a_E.getRegion(), num_species);
  EBCellFAB phot_src(ebisbox, a_E.getRegion(), num_Photons);

  part_src.setVal(0.0);
  phot_src.setVal(0.0);

  const BaseFab<Real>& EFab     = a_E.getSingleValuedFAB();

  
  // Grid loop
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    if(ebisbox.isRegular(iv)){
    
      // Position and E
      pos    = origin + RealVect(iv)*a_dx;
      E      = RealVect(D_DECL(EFab(iv, 0),     EFab(iv, 1),     EFab(iv, 2)));

      // Fill vectors with densities
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_cdrDensities[idx]).getSingleValuedFAB()(iv, 0);
	cdrDensities[idx] = Max(zero, phi);
	cdrGradients[idx] = RealVect(D_DECL((*a_cdrGradients[idx]).getSingleValuedFAB()(iv, 0),
						  (*a_cdrGradients[idx]).getSingleValuedFAB()(iv, 1),
						  (*a_cdrGradients[idx]).getSingleValuedFAB()(iv, 2)));
      }
      for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_rteDensities[idx]).getSingleValuedFAB()(iv, 0);
	rteDensities[idx] = Max(zero, phi);
      }

      // Compute source terms
      const Real kappa = 1.0; // Kappa for regular cells
      m_physics->advanceReactionNetwork(particle_sources,
					rteSources,
					cdrDensities,
					cdrGradients,
					rteDensities,
					E,
					pos,
					a_dx,
					a_dt,
					a_time,
					kappa);

      // Put vector into temporary holders
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	part_src.getSingleValuedFAB()(iv,idx) = particle_sources[idx];
      }
    
      // Put vector into temporary holders
      for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	phot_src.getSingleValuedFAB()(iv,idx) = rteSources[idx];
      }
    }
  }

  // Copy temporary storage back to solvers
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    (*a_cdrSources[idx]).setVal(0.0);
    (*a_cdrSources[idx]).plus(part_src, idx, 0, 1);

    // Covered cells are bogus. 
    (*a_cdrSources[idx]).setCoveredCellVal(0.0, 0);
  }

  // Copy temporary storage back to solvers
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    (*a_rteSources[idx]).setVal(0.0);
    (*a_rteSources[idx]).plus(phot_src, idx, 0, 1);

    // Covered cells are bogus. 
    (*a_rteSources[idx]).setCoveredCellVal(0.0, 0);
  }
}

void CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast(Vector<EBCellFAB*>&       a_cdrSources,
							      Vector<EBCellFAB*>&       a_rteSources,
							      const Vector<EBCellFAB*>& a_cdrDensities,
							      const Vector<EBCellFAB*>& a_cdrGradients,
							      const Vector<EBCellFAB*>& a_rteDensities,
							      const EBCellFAB&          a_E,
							      const Real&               a_time,
							      const Real&               a_dt,
							      const Real&               a_dx,
							      const Box&                a_box){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast(patch)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast(patch)" << endl;
  }

#if CH_SPACEDIM==2
  advanceReactionNetworkRegularCellsFast2D(a_cdrSources, a_rteSources, a_cdrDensities, a_cdrGradients,
					   a_rteDensities, a_E, a_time, a_dt, a_dx, a_box);
#elif CH_SPACEDIM==3
  advanceReactionNetworkRegularCellsFast3D(a_cdrSources, a_rteSources, a_cdrDensities, a_cdrGradients,
					   a_rteDensities, a_E, a_time, a_dt, a_dx, a_box);
#endif
  
}

void CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast2D(Vector<EBCellFAB*>&       a_cdrSources,
								Vector<EBCellFAB*>&       a_rteSources,
								const Vector<EBCellFAB*>& a_cdrDensities,
								const Vector<EBCellFAB*>& a_cdrGradients,
								const Vector<EBCellFAB*>& a_rteDensities,
								const EBCellFAB&          a_E,
								const Real&               a_time,
								const Real&               a_dt,
								const Real&               a_dx,
								const Box&                a_box){
#if CH_SPACEDIM==2
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast2D(patch)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast2D(patch)" << endl;
  }

  const Real zero = 0.0;
    

  const int num_species  = m_physics->getNumCdrSpecies();
  const int num_Photons  = m_physics->getNumRtSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics. 
  RealVect         pos = RealVect::Zero;
  RealVect         E   = RealVect::Zero;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     cdrDensities(num_species);
  Vector<RealVect> cdrGradients(num_species);
  Vector<Real>     rteSources(num_Photons);
  Vector<Real>     rteDensities(num_Photons);

  // I need contiguous memory for what is about to happen, so begin by copying stuff onto smaller data holders
  
  // Temps for source terms and particle densities
  FArrayBox cdr_src(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_gx(a_box, num_species);  // Gradient in x-direction
  FArrayBox cdr_gy(a_box, num_species);  // Gradient in y-direction
  cdr_phi.setVal(0.0);
  cdr_src.setVal(0.0);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_cdrDensities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
    cdr_gx.copy(a_cdrGradients[i]->getFArrayBox(),  a_box, 0, a_box, i, 1);
    cdr_gy.copy(a_cdrGradients[i]->getFArrayBox(),  a_box, 1, a_box, i, 1);
  }

  // Temps for Photon source terms and densities
  FArrayBox rte_phi(a_box, num_Photons);
  FArrayBox rte_src(a_box, num_Photons);
  rte_phi.setVal(0.0);
  rte_src.setVal(0.0);
  for (int i = 0; i < num_Photons; i++){
    rte_phi.copy(a_rteDensities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }

  // Temp for electric field
  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int offset     = lo[0] + n0*lo[1];

  // C style variable-length array conversion magic
  auto vla_cdr_src = (Real (*__restrict__)[n1][n0]) (cdr_src.dataPtr() - offset);
  auto vla_rte_src = (Real (*__restrict__)[n1][n0]) (rte_src.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_rte_phi = (Real (*__restrict__)[n1][n0]) (rte_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n1][n0]) (Efab.dataPtr()    - offset);
  auto vla_cdr_gx  = (Real (*__restrict__)[n1][n0]) (cdr_gx.dataPtr()  - offset);
  auto vla_cdr_gy  = (Real (*__restrict__)[n1][n0]) (cdr_gy.dataPtr()  - offset);

  for (int j = lo[1]; j <= hi[1]; ++j){
    for (int i = lo[0]; i <= hi[0]; ++i){
      
      // Particle densities
      for (int idx = 0; idx < num_species; ++idx){
	cdrDensities[idx] = Max(0.0, vla_cdr_phi[idx][j][i]);
	cdrGradients[idx][0] = vla_cdr_gx[idx][j][i];
	cdrGradients[idx][1] = vla_cdr_gy[idx][j][i];
      }

      // Photon densities
      for (int idx = 0; idx < num_Photons; ++idx){
	rteDensities[idx] = Max(0.0, vla_rte_phi[idx][j][i]);
      }

      E   = RealVect(vla_E[0][j][i], vla_E[1][j][i]);
      pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

      m_physics->advanceReactionNetwork(particle_sources,
					rteSources,
					cdrDensities,
					cdrGradients,
					rteDensities,
					E,
					pos,
					a_dx,
					a_dt,
					a_time,
					1.0);

      // Put result in correct palce
      for (int idx = 0; idx < num_species; ++idx){
	vla_cdr_src[idx][j][i] = particle_sources[idx];
      }

      // Put result in correct palce
      for (int idx = 0; idx < num_Photons; ++idx){
	vla_rte_src[idx][j][i] = rteSources[idx];
      }
    }
  }

  // Copy result back to solvers
  for (int i = 0; i < num_species; i++){
    FArrayBox& src = a_cdrSources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(cdr_src, a_box, i, a_box, 0, 1);
    a_cdrSources[i]->setCoveredCellVal(0.0, 0);
  }

  // Copy result back to solvers
  for (int i = 0; i < num_Photons; i++){
    FArrayBox& src = a_rteSources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(rte_src, a_box, i, a_box, 0, 1);
  }
#endif
}


void CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast3D(Vector<EBCellFAB*>&       a_cdrSources,
								Vector<EBCellFAB*>&       a_rteSources,
								const Vector<EBCellFAB*>& a_cdrDensities,
								const Vector<EBCellFAB*>& a_cdrGradients,
								const Vector<EBCellFAB*>& a_rteDensities,
								const EBCellFAB&          a_E,
								const Real&               a_time,
								const Real&               a_dt,
								const Real&               a_dx,
								const Box&                a_box){
#if CH_SPACEDIM==3
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast3D(patch)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkRegularCellsFast3D(patch)" << endl;
  }

  const Real zero = 0.0;

  const int num_species  = m_physics->getNumCdrSpecies();
  const int num_Photons  = m_physics->getNumRtSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics. 
  RealVect         pos = RealVect::Zero;
  RealVect         E   = RealVect::Zero;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     cdrDensities(num_species);
  Vector<RealVect> cdrGradients(num_species);
  Vector<Real>     rteSources(num_Photons);
  Vector<Real>     rteDensities(num_Photons);

  // I need contiguous memory for what is about to happen, so begin by copying stuff onto smaller data holders
  
  // Temps for source terms and particle densities
  FArrayBox cdr_src(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_gx(a_box, num_species);  // Gradient in x-direction
  FArrayBox cdr_gy(a_box, num_species);  // Gradient in y-direction
  FArrayBox cdr_gz(a_box, num_species);  // Gradient in z-direction
  cdr_phi.setVal(0.0);
  cdr_src.setVal(0.0);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_cdrDensities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
    cdr_gx.copy(a_cdrGradients[i]->getFArrayBox(),  a_box, 0, a_box, i, 1);
    cdr_gy.copy(a_cdrGradients[i]->getFArrayBox(),  a_box, 1, a_box, i, 1);
    cdr_gz.copy(a_cdrGradients[i]->getFArrayBox(),  a_box, 2, a_box, i, 1);
  }

  // Temps for Photon source terms and densities
  FArrayBox rte_phi(a_box, num_Photons);
  FArrayBox rte_src(a_box, num_Photons);
  rte_phi.setVal(0.0);
  rte_src.setVal(0.0);
  for (int i = 0; i < num_Photons; i++){
    rte_phi.copy(a_rteDensities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }

  // Temp for electric field
  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int n2         = dims[2];
  const int offset     = lo[0] + n0*(lo[1] + n1*lo[2]);

  // C style variable-length array conversion magic
  auto vla_cdr_src = (Real (*__restrict__)[n2][n1][n0]) (cdr_src.dataPtr() - offset);
  auto vla_rte_src = (Real (*__restrict__)[n2][n1][n0]) (rte_src.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n2][n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_rte_phi = (Real (*__restrict__)[n2][n1][n0]) (rte_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n2][n1][n0]) (Efab.dataPtr()    - offset);
  auto vla_cdr_gx  = (Real (*__restrict__)[n2][n1][n0]) (cdr_gx.dataPtr()  - offset);
  auto vla_cdr_gy  = (Real (*__restrict__)[n2][n1][n0]) (cdr_gy.dataPtr()  - offset);
  auto vla_cdr_gz  = (Real (*__restrict__)[n2][n1][n0]) (cdr_gz.dataPtr()  - offset);

  for (int k = lo[2]; k <= hi[2]; k++){
    for (int j = lo[1]; j <= hi[1]; ++j){
      for (int i = lo[0]; i <= hi[0]; ++i){
      
	// Particle densities
	for (int idx = 0; idx < num_species; ++idx){
	  cdrDensities[idx]    = Max(0.0, vla_cdr_phi[idx][k][j][i]);
	  cdrGradients[idx][0] = vla_cdr_gx[idx][k][j][i];
	  cdrGradients[idx][1] = vla_cdr_gy[idx][k][j][i];
	  cdrGradients[idx][2] = vla_cdr_gz[idx][k][j][i];
	}

	// Photon densities
	for (int idx = 0; idx < num_Photons; ++idx){
	  rteDensities[idx] = Max(0.0, vla_rte_phi[idx][k][j][i]);
	}

	E   = RealVect(vla_E[0][k][j][i], vla_E[1][k][j][i], vla_E[2][k][j][i]);
	pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

	m_physics->advanceReactionNetwork(particle_sources,
					  rteSources,
					  cdrDensities,
					  cdrGradients,
					  rteDensities,
					  E,
					  pos,
					  a_dx,
					  a_dt,
					  a_time,
					  1.0);

	// Put result in correct palce
	for (int idx = 0; idx < num_species; ++idx){
	  vla_cdr_src[idx][k][j][i] = particle_sources[idx];
	}

	// Put result in correct palce
	for (int idx = 0; idx < num_Photons; ++idx){
	  vla_rte_src[idx][k][j][i] = rteSources[idx];
	}
      }
    }
  }

  // Copy result back to solvers
  for (int i = 0; i < num_species; i++){
    FArrayBox& src = a_cdrSources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(cdr_src, a_box, i, a_box, 0, 1);
    a_cdrSources[i]->setCoveredCellVal(0.0, 0);
  }

  // Copy result back to solvers
  for (int i = 0; i < num_Photons; i++){
    FArrayBox& src = a_rteSources[i]->getFArrayBox();
    src.setVal(0.0);
    src.copy(rte_src, a_box, i, a_box, 0, 1);
  }

#endif
}



void CdrPlasmaStepper::advanceReactionNetworkIrreg(Vector<EBCellFAB*>&          a_cdrSources,
						   Vector<EBCellFAB*>&          a_rteSources,
						   const Vector<EBCellFAB*>&    a_cdrDensities,
						   const Vector<EBCellFAB*>&    a_cdrGradients,
						   const Vector<EBCellFAB*>&    a_particle_velocities,
						   const Vector<EBCellFAB*>&    a_rteDensities,
						   const BaseIVFAB<VoFStencil>& a_interp_stencils,
						   const EBCellFAB&             a_E,
						   const Real&                  a_time,
						   const Real&                  a_dt,
						   const Real&                  a_dx,
						   const Box&                   a_box,
						   const int                    a_lvl,
						   const DataIndex&             a_dit){
  switch(m_whichSourceTermComputation){
  case SourceTermComputation::Interpolated:
    {
      this->advanceReactionNetworkIrregInterp(a_cdrSources,
					      a_rteSources,
					      a_cdrDensities,
					      a_cdrGradients,
					      a_particle_velocities,
					      a_rteDensities,
					      a_interp_stencils,
					      a_E,
					      a_time,
					      a_dt,
					      a_dx,
					      a_box,
					      a_lvl,
					      a_dit);
      break;
    }
  case SourceTermComputation::InterpolatedStable:
    {
      this->advanceReactionNetworkIrregUpwind(a_cdrSources,
					      a_rteSources,
					      a_cdrDensities,
					      a_cdrGradients,
					      a_particle_velocities,					      
					      a_rteDensities,
					      a_interp_stencils,
					      a_E,
					      a_time,
					      a_dt,
					      a_dx,
					      a_box,
					      a_lvl,
					      a_dit);      
      break;	
    }
  case SourceTermComputation::CellAverage:
    {
      this->advanceReactionNetworkIrregKappa(a_cdrSources,
					     a_rteSources,
					     a_cdrDensities,
					     a_cdrGradients,
					     a_rteDensities,
					     a_interp_stencils,
					     a_E,
					     a_time,
					     a_dt,
					     a_dx,
					     a_box,
					     a_lvl,
					     a_dit);
      break;
    }
  case SourceTermComputation::Upwind:
    {
      this->advanceReactionNetworkIrregUpwind(a_cdrSources,
					      a_rteSources,
					      a_cdrDensities,
					      a_cdrGradients,
					      a_particle_velocities,					      
					      a_rteDensities,
					      a_interp_stencils,
					      a_E,
					      a_time,
					      a_dt,
					      a_dx,
					      a_box,
					      a_lvl,
					      a_dit);      
      break;	
    }
  default:
    MayDay::Error("CdrPlasmaStepper::advanceReactionNetworkIrreg - logic bust");
    break;
  }
}

void CdrPlasmaStepper::advanceReactionNetworkIrregInterp(Vector<EBCellFAB*>&          a_cdrSources,
							 Vector<EBCellFAB*>&          a_rteSources,
							 const Vector<EBCellFAB*>&    a_cdrDensities,
							 const Vector<EBCellFAB*>&    a_cdrGradients,
							 const Vector<EBCellFAB*>&    a_cdrVelocities,
							 const Vector<EBCellFAB*>&    a_rteIsotropicTerms,
							 const BaseIVFAB<VoFStencil>& a_interpStencils,
							 const EBCellFAB&             a_E,
							 const Real&                  a_time,
							 const Real&                  a_dt,
							 const Real&                  a_dx,
							 const Box&                   a_box,
							 const int                    a_lvl,
							 const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkIrregInterp(patch)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkIrregInterp(patch)" << endl;
  }

  constexpr int  comp = 0;
  constexpr Real zero = 0.0;
  
  const int numRteSolvers  = m_physics->getNumRtSpecies();
  const int numCdrSolvers  = m_physics->getNumCdrSpecies();

  // Things that are passed into CdrPlasmaPhysics
  RealVect         E;
  Vector<Real>     cdrSources      (numCdrSolvers, 0.            );
  Vector<Real>     cdrDensities    (numCdrSolvers, 0.            );
  Vector<Real>     rteSources      (numRteSolvers, 0.            );  
  Vector<Real>     rteIsotropicTerm(numRteSolvers, 0.            );
  Vector<RealVect> cdrGradients    (numCdrSolvers, RealVect::Zero);  

  VoFIterator&   vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];
  
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex&   vof     = vofit();
    const Real        kappa   = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interpStencils(vof, comp);
    const RealVect    pos     = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, a_dx);
    const RealVect    normal  = ebisbox.normal(vof);

    // I assume that the electric field is on the cell center (which is what FieldSolverMultigrid currently does), but we do want the field on the centroid.
    E = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += stencil.weight(i) * a_E(stencil.vof(i), dir);
      }
    }

    // Compute RTE densities on the centroids. Again, I sort of assume that the solver being used
    // is a cell-centered solver so we interpolate the isotropic term to the centroid. 
    for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // Apply the stencil but always leave a positive number of photons (interpolators might not always have non-negative weights). 
      Real tmp = 0.0;
      for (int i = 0; i < stencil.size(); i++){
	tmp += stencil.weight(i) * (*a_rteIsotropicTerms[idx])(stencil.vof(i), comp);
      }
      rteIsotropicTerm[idx] = std::max(tmp, zero);
    }

    // Compute cdr_densities and their gradients on centroids.
    for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      cdrDensities[idx] = zero;
      for (int i = 0; i < stencil.size(); i++){
	cdrDensities[idx] += stencil.weight(i) * (*a_cdrDensities[idx])(stencil.vof(i), comp);
      }
      cdrDensities[idx] = std::max(cdrDensities[idx], zero);

      // Interpolate gradients to centroids.
      cdrGradients[idx] = RealVect::Zero;
      for (int i = 0; i < stencil.size(); i++){
	for (int dir = 0; dir < SpaceDim; dir++){
	  cdrGradients[idx][dir] += stencil.weight(i) * (*a_cdrGradients[idx])(stencil.vof(i), dir);
	}
      }
    }


    // Call the plasma-kinetics framework and have it fill the source terms over a time step a_dt. 
    m_physics->advanceReactionNetwork(cdrSources,
				      rteSources,
				      cdrDensities,
				      cdrGradients,
				      rteIsotropicTerm,
				      E,
				      pos,
				      a_dx,
				      a_dt,
				      a_time,
				      kappa);

    // Iterate through the CDR and RTE solvers and set the source terms accordingly. 
    for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      (*a_cdrSources[idx])(vof, comp) = cdrSources[idx];
    }
    
    for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      (*a_rteSources[idx])(vof, comp) = rteSources[idx];
    }
  }
}

void CdrPlasmaStepper::advanceReactionNetworkIrregKappa(Vector<EBCellFAB*>&          a_cdrSources,
							Vector<EBCellFAB*>&          a_rteSources,
							const Vector<EBCellFAB*>&    a_cdrDensities,
							const Vector<EBCellFAB*>&    a_cdrGradients,
							const Vector<EBCellFAB*>&    a_rteDensities,
							const BaseIVFAB<VoFStencil>& a_interp_stencils,
							const EBCellFAB&             a_E,
							const Real&                  a_time,
							const Real&                  a_dt,
							const Real&                  a_dx,
							const Box&                   a_box,
							const int                    a_lvl,
							const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkIrregKappa(patch)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkIrregKappa(patch)" << endl;
  }

  const Real zero = 0.0;
  
  const int num_Photons  = m_physics->getNumRtSpecies();
  const int num_species  = m_physics->getNumCdrSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  RealVect         pos, E;
  Vector<Real>     particle_sources(num_species);
  Vector<Real>     rteSources(num_Photons);
  Vector<Real>     cdrDensities(num_species);
  Vector<RealVect> cdrGradients(num_species);
  Vector<Real>     rteDensities(num_Photons);

  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof       = vofit();
    const Real kappa          = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interp_stencils(vof, 0);

    pos = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);

    // Compute electric field on centroids
    E = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){
      const VolIndex& ivof = stencil.vof(i);
      const Real& iweight  = stencil.weight(i);
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += a_E(ivof, dir)*iweight;
      }
    }

    // Compute cdr_densities and their gradients on centroids
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      // Gradient on centroids for some reason
      RealVect grad = RealVect::Zero;
      for (int i = 0; i < stencil.size(); i++){
	const VolIndex& ivof = stencil.vof(i);
	const Real& iweight  = stencil.weight(i);
	for (int dir = 0; dir < SpaceDim; dir++){
	  grad[dir] += (*a_cdrGradients[idx])(ivof, dir)*iweight;
	}
      }
      
      cdrDensities[idx] = Max(zero, (*a_cdrDensities[idx])(vof,0));
      cdrGradients[idx] = grad;
    }

    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      rteDensities[idx] = Max(zero, (*a_rteDensities[idx])(vof, 0));
    }

    // Compute source terms
    m_physics->advanceReactionNetwork(particle_sources,
				      rteSources,
				      cdrDensities,
				      cdrGradients,
				      rteDensities,
				      E,
				      pos,
				      a_dx,
				      a_dt,
				      a_time,
				      kappa);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_cdrSources[idx])(vof, 0) = particle_sources[idx];
    }
    
    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_rteSources[idx])(vof, 0) = rteSources[idx];
    }
  }
}

void CdrPlasmaStepper::advanceReactionNetworkIrregUpwind(Vector<EBCellFAB*>&          a_cdrSources,
							 Vector<EBCellFAB*>&          a_rteSources,
							 const Vector<EBCellFAB*>&    a_cdrDensities,
							 const Vector<EBCellFAB*>&    a_cdrGradients,
							 const Vector<EBCellFAB*>&    a_cdrVelocities,
							 const Vector<EBCellFAB*>&    a_rteIsotropicTerms,
							 const BaseIVFAB<VoFStencil>& a_interpStencils,
							 const EBCellFAB&             a_E,
							 const Real&                  a_time,
							 const Real&                  a_dt,
							 const Real&                  a_dx,
							 const Box&                   a_box,
							 const int                    a_lvl,
							 const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkIrregUpwind(patch)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkIrregUpwind(patch)" << endl;
  }

  constexpr int  comp = 0;
  constexpr Real zero = 0.0;
  
  const int numRteSolvers  = m_physics->getNumRtSpecies();
  const int numCdrSolvers  = m_physics->getNumCdrSpecies();

  // Things that are passed into CdrPlasmaPhysics
  RealVect         E;
  Vector<Real>     cdrSources      (numCdrSolvers, 0.            );
  Vector<Real>     cdrDensities    (numCdrSolvers, 0.            );
  Vector<Real>     rteSources      (numRteSolvers, 0.            );  
  Vector<Real>     rteIsotropicTerm(numRteSolvers, 0.            );
  Vector<RealVect> cdrGradients    (numCdrSolvers, RealVect::Zero);  

  VoFIterator&   vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];
  
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex&   vof     = vofit();
    const Real        kappa   = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interpStencils(vof, comp);
    const RealVect    pos     = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, a_dx);

    // I assume that the electric field is on the cell center (which is what FieldSolverMultigrid currently does), but we do want the field on the centroid.
    E = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += stencil.weight(i) * a_E(stencil.vof(i), dir);
      }
    }

    // Compute RTE densities on the centroids. Again, I sort of assume that the solver being used
    // is a cell-centered solver so we interpolate the isotropic term to the centroid. 
    for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // Apply the stencil but always leave a positive number of photons (interpolators might not always have non-negative weights). 
      Real tmp = 0.0;
      for (int i = 0; i < stencil.size(); i++){
	tmp += stencil.weight(i) * (*a_rteIsotropicTerms[idx])(stencil.vof(i), comp);
      }
      rteIsotropicTerm[idx] = std::max(tmp, zero);
    }

    // Compute cdr_densities and their gradients on centroids. Note that since we upwind, we check if the
    // flow is into or away from the boundary. If the flow is away from the boundary we don't really have an upwind
    // side so we set the density to zero in that case. This is not really captured by CdrSolver::weightedUpwind because
    // the fallback option is to use the cell-centered value (that is the correct design, when we don't have reactive plasmas). We fix
    // that here.
    const Real EdotN = PolyGeom::dot(E, ebisbox.normal(vof));
    for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      const RefCountedPtr<CdrSpecies> species = solverIt.getSpecies();
      const int Z = species->getChargeNumber();
      
      bool inflow = false;
      if(solverIt()->isMobile()){
	if(Real(Z) * EdotN >= 0.0){
	  inflow = true;
	}
      }

      cdrDensities[idx] = zero;
      if(!inflow){
	for (int i = 0; i < stencil.size(); i++){
	  cdrDensities[idx] += stencil.weight(i) * (*a_cdrDensities[idx])(stencil.vof(i), comp);
	}
	cdrDensities[idx] = (*a_cdrDensities[idx])(vof, comp);
	cdrDensities[idx] = std::max(cdrDensities[idx], zero);
      }

      // Interpolate gradients to centroids.
      cdrGradients[idx] = RealVect::Zero;      
      for (int i = 0; i < stencil.size(); i++){
	for (int dir = 0; dir < SpaceDim; dir++){
	  cdrGradients[idx][dir] += stencil.weight(i) * (*a_cdrGradients[idx])(stencil.vof(i), dir);
	}
      }
    }


    // Call the plasma-kinetics framework and have it fill the source terms over a time step a_dt. 
    m_physics->advanceReactionNetwork(cdrSources,
				      rteSources,
				      cdrDensities,
				      cdrGradients,
				      rteIsotropicTerm,
				      E,
				      pos,
				      a_dx,
				      a_dt,
				      a_time,
				      kappa);

    // Iterate through the CDR and RTE solvers and set the source terms accordingly. 
    for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      (*a_cdrSources[idx])(vof, comp) = cdrSources[idx];
    }
    
    for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      (*a_rteSources[idx])(vof, comp) = rteSources[idx];
    }
  }
}

void CdrPlasmaStepper::computeCdrDiffusionFace(Vector<EBAMRFluxData*>&       a_diffusionCoefficient_face,
					       const Vector<EBAMRCellData*>& a_cdr_densities,
					       const EBAMRCellData&          a_E,
					       const Real&                   a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionFace(full)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionFace(full)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->getNumCdrSpecies();
  const int finest_level = m_amr->getFinestLevel();

  // Allocate data for cell-centered diffusion coefficients
  Vector<EBAMRCellData> diffco(num_species);
  Vector<Real> cdr_densities(num_species);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();

    m_amr->allocate(diffco[idx], m_realm, m_cdr->getPhase(), ncomp);
  }

  // Call the cell version
  computeCdrDiffusionCell(diffco, a_cdr_densities, a_E, a_time);

  // Now compute face-centered things by taking the average of cell-centered things
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();

    if(solver->isDiffusive()){ // Only need to do this for diffusive things
      DataOps::setValue(*a_diffusionCoefficient_face[idx], std::numeric_limits<Real>::max());
      m_amr->averageDown  (diffco[idx], m_realm, m_cdr->getPhase());
      m_amr->interpGhostMG(diffco[idx], m_realm, m_cdr->getPhase()); 

#if 0 // Doesn't fill ghost faces
      DataOps::averageCellToFace(*a_diffusionCoefficient_face[idx], diffco[idx], m_amr->getDomains());
#else // Fills ghost faces
      DataOps::averageCellScalarToFaceScalar(*a_diffusionCoefficient_face[idx], diffco[idx], m_amr->getDomains());
#endif
    }
  }
}

void CdrPlasmaStepper::computeCdrDiffusionCell(Vector<EBAMRCellData>&       a_diffusionCoefficient_cell,
					       const Vector<EBAMRCellData*>& a_cdr_densities,
					       const EBAMRCellData&          a_E,
					       const Real&                   a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCell(amr)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCell(amr)" << endl;
  }
  
  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->getNumCdrSpecies();
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<LevelData<EBCellFAB>* > diffco(num_species);
    Vector<LevelData<EBCellFAB>* > cdr_densities(num_species);
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx] = a_diffusionCoefficient_cell[idx][lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    // Cell the level version
    computeCdrDiffusionCell(diffco, cdr_densities, *a_E[lvl], lvl, a_time);
  }
}

void CdrPlasmaStepper::computeCdrDiffusionCell(Vector<LevelData<EBCellFAB>* >&       a_diffusionCoefficient_cell,
					       const Vector<LevelData<EBCellFAB>* >& a_cdr_densities,
					       const LevelData<EBCellFAB>&           a_E,
					       const int                             a_lvl,
					       const Real&                           a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCell(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCell(level)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->getNumCdrSpecies();

  const IrregAmrStencil<CentroidInterpolationStencil>& interp_stencils = m_amr->getCentroidInterpolationStencils(m_realm, m_cdr->getPhase());

  // Call the level version
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const RealVect origin         = m_amr->getProbLo();
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = IntVectSet(box);
    const EBCellFAB& E     = a_E[dit()];

    Vector<EBCellFAB*> diffco(num_species);
    Vector<EBCellFAB*> cdr_densities(num_species);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx]        = &(*a_diffusionCoefficient_cell[idx])[dit()];
      cdr_densities[idx] = &(*a_cdr_densities[idx])[dit()];
    }

    // Regular cells
#if USE_FAST_DIFFUSION
    computeCdrDiffusionCellRegularFast(diffco, cdr_densities, E, dbl.get(dit()), m_amr->getDx()[a_lvl], a_time);
#else
    computeCdrDiffusionCellRegular(diffco, cdr_densities, E, dbl.get(dit()), m_amr->getDx()[a_lvl], a_time);
#endif

    // Irregular cells
    computeCdrDiffusionCellIrregular(diffco, cdr_densities, E, dbl.get(dit()), m_amr->getDx()[a_lvl],
				     interp_stencils[a_lvl][dit()], a_time, a_lvl, dit());


  }
}

void CdrPlasmaStepper::computeCdrDiffusionCellRegular(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
						      const Vector<EBCellFAB*>& a_cdr_densities,
						      const EBCellFAB&          a_E,
						      const Box                 a_box,
						      const Real                a_dx,
						      const Real                a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCellRegular");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCellRegular" << endl;
  }

  const int comp        = 0;
  const int num_species = m_physics->getNumCdrSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  RealVect         pos, E;
  Vector<Real>     cdr_densities(num_species);

  // Computed coefficients terms onto here
  EBCellFAB tmp(ebisbox, a_E.getRegion(), num_species);

  const BaseFab<Real>& EFab  = a_E.getSingleValuedFAB();
  
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    
    pos = origin + iv*a_dx;
    E   = RealVect(D_DECL(EFab(iv, 0),     EFab(iv, 1),     EFab(iv, 2)));

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      const Real phi = (*a_cdr_densities[idx]).getSingleValuedFAB()(iv, comp);
      cdr_densities[idx] = Max(0.0, phi);
    }

    const Vector<Real> coeffs = m_physics->computeCdrDiffusionCoefficients(a_time,
									   pos,
									   E,
									   cdr_densities);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->isDiffusive()){
	tmp.getSingleValuedFAB()(iv,idx) = coeffs[idx];
      }
    }
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    if(solver_it()->isDiffusive()){
      (*a_diffusionCoefficient_cell[idx]).setVal(0.0);
      (*a_diffusionCoefficient_cell[idx]).plus(tmp, idx, 0, 1);

      // Covered cells are bogus. 
      (*a_diffusionCoefficient_cell[idx]).setCoveredCellVal(0.0, 0);
    }
  }
}

void CdrPlasmaStepper::computeCdrDiffusionCellRegularFast(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
							  const Vector<EBCellFAB*>& a_cdr_densities,
							  const EBCellFAB&          a_E,
							  const Box                 a_box,
							  const Real                a_dx,
							  const Real                a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCellRegularFast");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCellRegularFast" << endl;
  }

#if CH_SPACEDIM==2
  computeCdrDiffusionCellRegularFast2D(a_diffusionCoefficient_cell, a_cdr_densities, a_E, a_box, a_dx, a_time);
#elif CH_SPACEDIM==3
  computeCdrDiffusionCellRegularFast3D(a_diffusionCoefficient_cell, a_cdr_densities, a_E, a_box, a_dx, a_time);
#endif
}

void CdrPlasmaStepper::computeCdrDiffusionCellRegularFast2D(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
							    const Vector<EBCellFAB*>& a_cdr_densities,
							    const EBCellFAB&          a_E,
							    const Box                 a_box,
							    const Real                a_dx,
							    const Real                a_time){
#if CH_SPACEDIM==2
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCellRegularFast2D");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCellRegularFast2D" << endl;
  }

  const int comp        = 0;
  const int num_species = m_physics->getNumCdrSpecies();
  const RealVect origin = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  RealVect     pos, E;
  Vector<Real> cdr_densities(num_species);

  // I need contiguous memory for the nasty that is about to happen. So begin by copying things onto smaller data holders
  FArrayBox cdr_dco(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox Efab(a_box, SpaceDim);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_cdr_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims = a_box.size();
  const IntVect lo   = a_box.smallEnd();
  const IntVect hi   = a_box.bigEnd();
  const int n0       = dims[0];
  const int n1       = dims[1];
  const int offset   = lo[0] + n0*lo[1];

  // C style variable-length array conversion magic
  auto vla_cdr_dco = (Real (*__restrict__)[n1][n0]) (cdr_dco.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n1][n0]) (Efab.dataPtr()    - offset);

  for (int j = lo[1]; j <= hi[1]; ++j){
    for (int i = lo[0]; i <= hi[0]; ++i){
      
      // Particle densities
      for (int idx = 0; idx < num_species; ++idx){
	cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][j][i]);
      }

      E   = RealVect(vla_E[0][j][i], vla_E[1][j][i]);
      pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

      const Vector<Real> coeffs = m_physics->computeCdrDiffusionCoefficients(a_time,
									     pos,
									     E,
									     cdr_densities);

      // Put result in correct palce
      for (int idx = 0; idx < num_species; ++idx){
	vla_cdr_dco[idx][j][i] = coeffs[idx];
      }
    }
  }
  
  // Copy result back to solvers
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    if(solver->isDiffusive()){
      const int idx = solver_it.index();
      FArrayBox& dco = a_diffusionCoefficient_cell[idx]->getFArrayBox();
      dco.setVal(0.0);
      dco.copy(cdr_dco, a_box, idx, a_box, 0, 1);
      a_diffusionCoefficient_cell[idx]->setCoveredCellVal(0.0, 0);
    }
  }
#endif
}

void CdrPlasmaStepper::computeCdrDiffusionCellRegularFast3D(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
							    const Vector<EBCellFAB*>& a_cdr_densities,
							    const EBCellFAB&          a_E,
							    const Box                 a_box,
							    const Real                a_dx,
							    const Real                a_time){
#if CH_SPACEDIM==3
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCellRegularFast3D");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCellRegularFast3D" << endl;
  }

  const int comp        = 0;
  const int num_species = m_physics->getNumCdrSpecies();
  const RealVect origin = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  RealVect     pos, E;
  Vector<Real> cdr_densities(num_species);

  // I need contiguous memory for the nasty that is about to happen. So begin by copying things onto smaller data holders
  FArrayBox cdr_dco(a_box, num_species);
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox Efab(a_box, SpaceDim);
  for (int i = 0; i < num_species; i++){
    cdr_phi.copy(a_cdr_densities[i]->getFArrayBox(), a_box, 0, a_box, i, 1);
  }
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);

  // Pointer offsets
  const IntVect dims = a_box.size();
  const IntVect lo   = a_box.smallEnd();
  const IntVect hi   = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int n2         = dims[2];
  const int offset     = lo[0] + n0*(lo[1] + n1*lo[2]);


  // C style variable-length array conversion magic
  auto vla_cdr_dco = (Real (*__restrict__)[n2][n1][n0]) (cdr_dco.dataPtr() - offset);
  auto vla_cdr_phi = (Real (*__restrict__)[n2][n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_E       = (Real (*__restrict__)[n2][n1][n0]) (Efab.dataPtr()    - offset);

  for (int k = lo[2]; k <= hi[2]; k++){
    for (int j = lo[1]; j <= hi[1]; ++j){
      for (int i = lo[0]; i <= hi[0]; ++i){
      
	// Particle densities
	for (int idx = 0; idx < num_species; ++idx){
	  cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][k][j][i]);
	}
	
	E   = RealVect(vla_E[0][k][j][i], vla_E[1][k][j][i], vla_E[2][k][j][i]);
	pos = origin + RealVect(D_DECL(i,j,k))*a_dx;

	const Vector<Real> coeffs = m_physics->computeCdrDiffusionCoefficients(a_time,
									       pos,
									       E,
									       cdr_densities);
	
	// Put result in correct palce
	for (int idx = 0; idx < num_species; ++idx){
	  vla_cdr_dco[idx][k][j][i] = coeffs[idx];
	}
      }
    }
  }
    
  // Copy result back to solvers
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    if(solver->isDiffusive()){
      const int idx = solver_it.index();
      FArrayBox& dco = a_diffusionCoefficient_cell[idx]->getFArrayBox();
      dco.setVal(0.0);
      dco.copy(cdr_dco, a_box, idx, a_box, 0, 1);

      a_diffusionCoefficient_cell[idx]->setCoveredCellVal(0.0, 0);
    }
  }
#endif
}

void CdrPlasmaStepper::computeCdrDiffusionCellIrregular(Vector<EBCellFAB*>&          a_diffusionCoefficient_cell,
							const Vector<EBCellFAB*>&    a_cdr_densities,
							const EBCellFAB&             a_E,
							const Box                    a_box,
							const Real                   a_dx,
							const BaseIVFAB<VoFStencil>& a_interp_stencils,
							const Real&                  a_time,
							const int                    a_lvl,
							const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCellIrregular");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCellIrregular" << endl;
  }

  const int comp = 0;
  const int num_species  = m_physics->getNumCdrSpecies();

  // EBISBox and graph
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  RealVect         pos, E;
  Vector<Real>     cdr_densities(num_species);

  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
      
    const RealVect pos = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);
    const RealVect E   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));

    Vector<Real> cdr_densities(num_species);
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities[idx] = (*a_cdr_densities[idx])(vof, comp);
    }
      
    const Vector<Real> coeffs = m_physics->computeCdrDiffusionCoefficients(a_time,
									   pos,
									   E,
									   cdr_densities);
      
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->isDiffusive()){
	(*a_diffusionCoefficient_cell[idx])(vof, comp) = coeffs[idx];
      }
    }
  }
}

void CdrPlasmaStepper::computeCdrDiffusionEb(Vector<EBAMRIVData*>&       a_ebDiffusionCoefficient,
					     const Vector<EBAMRIVData*>& a_cdr_densities,
					     const EBAMRIVData&          a_E,
					     const Real&                 a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionEb(full)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionEb(full)" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int num_species  = m_physics->getNumCdrSpecies();
  const int finest_level = m_amr->getFinestLevel();
  const Real zero        = 0.0;

  Vector<LevelData<BaseIVFAB<Real> >* > diffco(num_species);
  Vector<LevelData<BaseIVFAB<Real> >* > cdr_densities(num_species);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<BaseIVFAB<Real> >* > diffco(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > cdr_densities(num_species);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      diffco[idx]        = (*a_ebDiffusionCoefficient[idx])[lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    // Call level versions
    computeCdrDiffusionEb(diffco, cdr_densities, *a_E[lvl], a_time, lvl);
  }
}


void CdrPlasmaStepper::computeCdrDiffusionEb(Vector<LevelData<BaseIVFAB<Real> >* >&       a_ebDiffusionCoefficient,
					     const Vector<LevelData<BaseIVFAB<Real> >* >& a_cdr_densities,
					     const LevelData<BaseIVFAB<Real> >&           a_E,
					     const Real&                                  a_time,
					     const int                                    a_lvl){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionEb(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionEb(level)" << endl;
  }

  const int comp = 0;
  const int num_species = m_physics->getNumCdrSpecies();

  
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const RealVect origin         = m_amr->getProbLo();

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box            = dbl.get(dit());
    const EBISBox& ebisbox   = ebisl[dit()];
    const EBGraph& ebgraph   = ebisbox.getEBGraph();
    const IntVectSet ivs     = ebisbox.getIrregIVS(box);
    const BaseIVFAB<Real>& E = a_E[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const RealVect cntr = ebisbox.bndryCentroid(vof);
      const RealVect e    = RealVect(D_DECL(E(vof,0), E(vof,1), E(vof,2)));
      const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin) + cntr*dx;
      
      Vector<Real> cdr_densities(num_species);
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx  = solver_it.index();
	const Real phi = (*a_cdr_densities[idx])[dit()](vof, 0);
	cdr_densities[idx] = Max(0.0, phi);
      }


      const Vector<Real> diffco = m_physics->computeCdrDiffusionCoefficients(a_time,
									     pos,
									     e,
									     cdr_densities);
										  
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	if(solver_it()->isDiffusive()){
#if 0 // Original code
	  (*a_ebDiffusionCoefficient[idx])[dit()](vof, comp) = diffco[idx];
#else // Debug code
	  (*a_ebDiffusionCoefficient[idx])[dit()](vof, comp) = 0.0;//1.E-5;
#endif
	}
      }
    }
  }
}

void CdrPlasmaStepper::computeCdrFluxes(Vector<LevelData<BaseIVFAB<Real> >*>&       a_fluxes,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_fluxes,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_densities,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_velocities,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_gradients,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_rte_fluxes,
					const LevelData<BaseIVFAB<Real> >&          a_E,
					const Real&                                 a_time,
					const int                                   a_lvl){
  CH_TIME("CdrPlasmaStepper::computeCdrFluxes(full, level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrFluxes(full, level)" << endl;
  }

  const int num_species  = m_physics->getNumCdrSpecies();
  const int num_Photons  = m_physics->getNumRtSpecies();
  const int comp         = 0;
  const int ncomp        = 1;

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_Photons);

  // Grid stuff
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const ProblemDomain& domain   = m_amr->getDomains()[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const MFLevelGrid& mflg       = *(m_amr->getMFLevelGrid(m_realm)[a_lvl]);
  const RealVect origin         = m_amr->getProbLo();

  // Patch loop
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box              = dbl.get(dit());
    const EBISBox& ebisbox      = ebisl[dit()];
    const EBGraph& ebgraph      = ebisbox.getEBGraph();
    const IntVectSet& diel_ivs  = mflg.interfaceRegion(box, dit());
    const IntVectSet& elec_ivs  = ebisbox.getIrregIVS(box) - diel_ivs;

    // Loop over conductor cells
    for (VoFIterator vofit(elec_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();


      // Define the electric field
      const BaseIVFAB<Real>& E = a_E[dit()];
      const RealVect centroid  = ebisbox.bndryCentroid(vof);
      const RealVect normal    = ebisbox.normal(vof);
      const RealVect e         = RealVect(D_DECL(E(vof,0), E(vof,1), E(vof,2)));
      const RealVect pos       = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin) + centroid*dx;

      // Build ion densities and velocities
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](vof,comp);
	extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](vof,comp);
	extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](vof,comp);
	extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](vof,comp);
      }

      // Build Photon intensities
      for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[dit()](vof,comp);
      }

      const Vector<Real> fluxes = m_physics->computeCdrElectrodeFluxes(a_time,
								       pos,
								       normal,
								       e,
								       extrap_cdr_densities,
								       extrap_cdr_velocities,
								       extrap_cdr_gradients,
								       extrap_rte_fluxes,
								       extrap_cdr_fluxes);
	
      // Put the fluxes in their respective place
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	(*a_fluxes[idx])[dit()](vof, comp) = fluxes[idx];
      }
    }

    // Loop over dielectric cells
    for (VoFIterator vofit(diel_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      // Define the electric field
      const BaseIVFAB<Real>& E = a_E[dit()];
      const RealVect centroid  = ebisbox.bndryCentroid(vof);
      const RealVect normal    = ebisbox.normal(vof);
      const RealVect e         = RealVect(D_DECL(E(vof,0), E(vof,1), E(vof,2)));
      const RealVect pos       = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin) + centroid*dx;
      const Real     time      = 0.0;

      // Build ion densities and velocities
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](vof,comp);
	extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](vof,comp);
	extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](vof,comp);
	extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](vof,comp);
      }

      // Build Photon intensities
      for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[dit()](vof,comp);
      }

      const Vector<Real> fluxes = m_physics->computeCdrDielectricFluxes(a_time,
									pos,
									normal,
									e,
									extrap_cdr_densities,
									extrap_cdr_velocities,
									extrap_cdr_gradients,
									extrap_rte_fluxes,
									extrap_cdr_fluxes);
	
      // Put the fluxes in their respective place
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	(*a_fluxes[idx])[dit()](vof, comp) = fluxes[idx];
      }
    }
  }
}

void CdrPlasmaStepper::computeCdrFluxes(Vector<EBAMRIVData*>&       a_fluxes,
					const Vector<EBAMRIVData*>& a_extrap_cdr_fluxes,
					const Vector<EBAMRIVData*>& a_extrap_cdr_densities,
					const Vector<EBAMRIVData*>& a_extrap_cdr_velocities,
					const Vector<EBAMRIVData*>& a_extrap_cdr_gradients,
					const Vector<EBAMRIVData*>& a_extrap_rte_fluxes,
					const EBAMRIVData&          a_E,
					const Real&                 a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrFluxes(full)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrFluxes(full)" << endl;
  }

  const int num_species  = m_physics->getNumCdrSpecies();
  const int num_Photons  = m_physics->getNumRtSpecies();
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){


    Vector<LevelData<BaseIVFAB<Real> >* > fluxes(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_fluxes(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_densities(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_velocities(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_cdr_gradients(num_species);
    Vector<LevelData<BaseIVFAB<Real> >* > extrap_rte_fluxes(num_Photons);
    
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      fluxes[idx]                = (*a_fluxes[idx])[lvl];
      extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[lvl];
      extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[lvl];
      extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[lvl];
      extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[lvl];
    }


    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[lvl];
    }

    // Call the level version
    computeCdrFluxes(fluxes, extrap_cdr_fluxes, extrap_cdr_densities, extrap_cdr_velocities,
		     extrap_cdr_gradients, extrap_rte_fluxes, *a_E[lvl], a_time, lvl);
  }
}

void CdrPlasmaStepper::computeCdrDomainFluxes(Vector<EBAMRIFData*>&       a_fluxes,
					      const Vector<EBAMRIFData*>& a_extrap_cdr_fluxes,
					      const Vector<EBAMRIFData*>& a_extrap_cdr_densities,
					      const Vector<EBAMRIFData*>& a_extrap_cdr_velocities,
					      const Vector<EBAMRIFData*>& a_extrap_cdr_gradients,
					      const Vector<EBAMRIFData*>& a_extrap_rte_fluxes,
					      const EBAMRIFData&          a_E,
					      const Real&                 a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDomainFluxes(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDomainFluxes(level)" << endl;
  }

  const int num_species  = m_physics->getNumCdrSpecies();
  const int num_Photons  = m_physics->getNumRtSpecies();
  const int finest_level = m_amr->getFinestLevel();

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_Photons);

  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<DomainFluxIFFAB>* > fluxes(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_fluxes(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_densities(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_velocities(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_cdr_gradients(num_species);
    Vector<LevelData<DomainFluxIFFAB>* > extrap_rte_fluxes(num_Photons);
    
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      fluxes[idx]                = (*a_fluxes[idx])[lvl];
      extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[lvl];
      extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[lvl];
      extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[lvl];
      extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[lvl];
    }

    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[lvl];
    }

    // Call the level version
    computeCdrDomainFluxes(fluxes, extrap_cdr_fluxes, extrap_cdr_densities, extrap_cdr_velocities,
			   extrap_cdr_gradients, extrap_rte_fluxes, *a_E[lvl], a_time, lvl);
    
  }
}

void CdrPlasmaStepper::computeCdrDomainFluxes(Vector<LevelData<DomainFluxIFFAB>*>        a_fluxes,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_fluxes,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_densities,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_velocities,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_gradients,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_rte_fluxes,
					      const LevelData<DomainFluxIFFAB>&          a_E,
					      const Real&                                a_time,
					      const int                                  a_lvl){
  CH_TIME("CdrPlasmaStepper::computeCdrDomainFluxes(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDomainFluxes(level)" << endl;
  }

  const int num_species  = m_physics->getNumCdrSpecies();
  const int num_Photons  = m_physics->getNumRtSpecies();
  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();

  // Things that will be passed into physics
  Vector<Real> extrap_cdr_fluxes(num_species);
  Vector<Real> extrap_cdr_densities(num_species);
  Vector<Real> extrap_cdr_velocities(num_species);
  Vector<Real> extrap_cdr_gradients(num_species);
  Vector<Real> extrap_rte_fluxes(num_Photons);

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const ProblemDomain& domain   = m_amr->getDomains()[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const RealVect origin         = m_amr->getProbLo();

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box         = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    
    const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
    
    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	const IntVectSet& ivs  = (*a_fluxes[0])[dit()](dir, sit()).getIVS();

	for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  const RealVect pos    = EBArith::getFaceLocation(face, dx*RealVect::Unit, origin);
	    
	  // Define the electric field
	  const RealVect E = RealVect(D_DECL((a_E)[dit()](dir, sit())(face,0),
					     (a_E)[dit()](dir, sit())(face,1),
					     (a_E)[dit()](dir, sit())(face,2)));

	  // Ion densities. 
	  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    extrap_cdr_fluxes[idx]     = (*a_extrap_cdr_fluxes[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_densities[idx]  = (*a_extrap_cdr_densities[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_velocities[idx] = (*a_extrap_cdr_velocities[idx])[dit()](dir, sit())(face,comp);
	    extrap_cdr_gradients[idx]  = (*a_extrap_cdr_gradients[idx])[dit()](dir, sit())(face,comp);
	  }

	  // Photon fluxes
	  for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    extrap_rte_fluxes[idx] = (*a_extrap_rte_fluxes[idx])[dit()](dir, sit())(face,comp);
	  }

	  // Call CdrPlasmaPhysics
	  const Vector<Real> fluxes = m_physics->computeCdrDomainFluxes(a_time,
									pos,
									dir,
									sit(),
									E,
									extrap_cdr_densities,
									extrap_cdr_velocities,
									extrap_cdr_gradients,
									extrap_rte_fluxes,
									extrap_cdr_fluxes);

	  // Put fluxes where they belong
	  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.index();
	    (*a_fluxes[idx])[dit()](dir, sit())(face, comp) = fluxes[idx];
	  }
	}
      }
    }
  }
}

void CdrPlasmaStepper::computeGradientsAtEb(Vector<EBAMRIVData*>&         a_grad,
					    const phase::which_phase&     a_phase,
					    const Vector<EBAMRCellData*>& a_phi){
  CH_TIME("CdrPlasmaStepper::computeGradientsAtEb");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeGradientsAtEb" << endl;
  }

  CH_assert(a_grad.size() == a_phi.size());

  
  EBAMRIVData   eb_gradient;
  EBAMRCellData gradient;
  EBAMRCellData scratch;  
  m_amr->allocate(eb_gradient, m_realm, a_phase, SpaceDim);
  m_amr->allocate(gradient,    m_realm, a_phase, SpaceDim);
  m_amr->allocate(scratch,     m_realm, a_phase, 1       );  

  for (int i = 0; i < a_phi.size(); i++){
    EBAMRIVData& grad_density    = *a_grad[i];
    const EBAMRCellData& density = *a_phi[i];

    CH_assert(grad_density[0]->nComp() == 1);
    CH_assert(density[0]->nComp()      == 1);

    DataOps::copy(scratch, density);
    m_amr->interpGhostMG(scratch, m_realm, a_phase);
    
    m_amr->computeGradient(gradient, scratch, m_realm, a_phase);       // Compute cell-centered gradient
    m_amr->averageDown(gradient, m_realm, a_phase);                    // Average down - shouldn't be necesasry
    m_amr->interpGhost(gradient, m_realm, a_phase);                    // Interpolate ghost cells (have to do this before interp)
    this->extrapolateToEb(eb_gradient, a_phase, gradient);            // Extrapolate to EB
    this->projectFlux(grad_density, eb_gradient);                      // Project onto EB
  }
}

void CdrPlasmaStepper::computeGradientsAtDomainFaces(Vector<EBAMRIFData*>&         a_grad,
						     const phase::which_phase&     a_phase,
						     const Vector<EBAMRCellData*>& a_phi){
  CH_TIME("CdrPlasmaStepper::computeGradientsAtEb");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeGradientsAtEb" << endl;
  }

  CH_assert(a_grad.size() == a_phi.size());

  EBAMRIFData   domain_gradient;
  EBAMRCellData gradient;
  EBAMRCellData scratch;
  m_amr->allocate(domain_gradient, m_realm, a_phase, SpaceDim);
  m_amr->allocate(gradient,        m_realm, a_phase, SpaceDim);
  m_amr->allocate(scratch,         m_realm, a_phase, 1       );  

  for (int i = 0; i < a_phi.size(); i++){
    EBAMRIFData& grad_density    = *a_grad[i];
    const EBAMRCellData& density = *a_phi[i];

    CH_assert(grad_density[0]->nComp() == 1);
    CH_assert(density[0]->nComp()      == 1);

    DataOps::copy(scratch, density);
    m_amr->interpGhostMG(scratch, m_realm, a_phase);
    m_amr->computeGradient(gradient, scratch, m_realm, a_phase);                         
    m_amr->averageDown(gradient, m_realm, a_phase);                             
    m_amr->interpGhost(gradient, m_realm, a_phase);                             
    
    this->extrapolateToDomainFaces(domain_gradient, a_phase, gradient);  // Extrapolate to EB
    this->projectDomain(grad_density, domain_gradient);                    // Project normal compoent
  }
}

void CdrPlasmaStepper::extrapolateToVectorDomainFaces(EBAMRIFData&             a_extrap,
						      const phase::which_phase a_phase,
						      const EBAMRCellData&     a_data){
  CH_TIME("CdrPlasmaStepper::extrapolateToDomainFaces");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToDomainFaces" << endl;
  }

  CH_assert(a_extrap[0]->nComp() == 1);
  CH_assert(a_data[0]->nComp()   == SpaceDim);

  EBAMRIFData domain_vector;
  m_amr->allocate(domain_vector, m_realm, a_phase, SpaceDim);
  this->extrapolateToDomainFaces(domain_vector, a_phase, a_data);
  this->projectDomain(a_extrap, domain_vector);
}

void CdrPlasmaStepper::extrapolateToVectorDomainFaces(Vector<EBAMRIFData*>&         a_extrap,
						      const phase::which_phase      a_phase,
						      const Vector<EBAMRCellData*>& a_data){
  CH_TIME("CdrPlasmaStepper::extrapolateToVectorDomainFaces");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToVectorDomainFaces" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    extrapolateToVectorDomainFaces(*a_extrap[idx], a_phase, *a_data[idx]);
  }
}

void CdrPlasmaStepper::extrapolateVelocitiesVectorDomainFaces(Vector<EBAMRIFData*>&         a_extrap,
							      const phase::which_phase      a_phase,
							      const Vector<EBAMRCellData*>& a_velocities){
  CH_TIME("CdrPlasmaStepper::extrapolateVelocitiesVectorDomainFaces");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateVelocitiesVectorDomainFaces)" << endl;
  }

  //  for(int i = 0; i < a_extrap.size(); i++){
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      extrapolateToVectorDomainFaces(*a_extrap[idx], a_phase, *a_velocities[idx]);
    }
    else{
      DataOps::setValue(*a_extrap[idx], 0.0);
    }
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocities(Vector<EBAMRCellData*>&       a_velocities,
						 const Vector<EBAMRCellData*>& a_cdr_densities,
						 const EBAMRCellData&          a_E,
						 const Real&                   a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocities(full)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocities(full)" << endl;
  }
  
  const int num_species  = m_physics->getNumCdrSpecies();
  const int finest_level = m_amr->getFinestLevel();

  // Interpolate E to centroids
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, phase::gas, SpaceDim);
  DataOps::copy(E, a_E);
  //  m_amr->interpToCentroids(E, m_realm, phase::gas);

  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<EBCellFAB>* > velocities(num_species);
    Vector<LevelData<EBCellFAB>* > cdr_densities(num_species);

    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      velocities[idx]    = (*a_velocities[idx])[lvl];
      cdr_densities[idx] = (*a_cdr_densities[idx])[lvl];
    }

    computeCdrDriftVelocities(velocities, cdr_densities, *E[lvl], lvl, a_time);
  }

  // Average down and interpolate ghost cells
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    if(solver_it()->isMobile()){
      m_amr->averageDown(*a_velocities[idx], m_realm, m_cdr->getPhase()); 
      m_amr->interpGhostMG(*a_velocities[idx], m_realm, m_cdr->getPhase());
    }
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocities(Vector<LevelData<EBCellFAB> *>&       a_velocities,
						 const Vector<LevelData<EBCellFAB> *>& a_cdr_densities,
						 const LevelData<EBCellFAB> &          a_E,
						 const int                             a_lvl,
						 const Real&                           a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocities(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocities(level)" << endl;
  }

  const phase::which_phase cdr_phase = m_cdr->getPhase();
  const int finest_level             = m_amr->getFinestLevel();

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, cdr_phase)[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];

  const int num_species = m_physics->getNumCdrSpecies();
    
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    Vector<EBCellFAB*> vel(num_species);
    Vector<EBCellFAB*> phi(num_species);;
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      if(solver_it()->isMobile()){
	vel[idx] = &(*a_velocities[idx])[dit()];
      }
      phi[idx] = &(*a_cdr_densities[idx])[dit()];
    }

    // Separate calls for regular and irregular
#if USE_FAST_VELOCITIES
    computeCdrDriftVelocitiesRegularFast(vel,   phi, a_E[dit()], dbl.get(dit()), a_time, dx);
#else
    computeCdrDriftVelocitiesRegular(vel,   phi, a_E[dit()], dbl.get(dit()), a_time, dx);
#endif
    computeCdrDriftVelocitiesIrregular(vel, phi, a_E[dit()], dbl.get(dit()), a_time, dx, a_lvl, dit());
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocitiesRegular(Vector<EBCellFAB*>&       a_velocities,
							const Vector<EBCellFAB*>& a_cdr_densities,
							const EBCellFAB&          a_E,
							const Box&                a_box,
							const Real&               a_time,
							const Real&               a_dx){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocitiesRegular");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocitiesRegular" << endl;
  }

  const int comp         = 0;
  const RealVect origin  = m_amr->getProbLo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();
  const EBISBox& ebisbox = a_E.getEBISBox();


  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv    = bit();
    const RealVect pos  = origin + a_dx*iv;
    const RealVect e    = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));

    Vector<RealVect> velocities(a_cdr_densities.size(), RealVect::Zero);

    if(ebisbox.isRegular(iv)){

      // Get densities
      Vector<Real> cdr_densities;
      for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	cdr_densities.push_back((*a_cdr_densities[idx]).getSingleValuedFAB()(iv, comp));
      }

      // Compute velocities
      velocities = m_physics->computeCdrDriftVelocities(a_time, pos, e, cdr_densities);
    }

    // Put velocities in the appropriate place. 
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver = solver_it();
      const int idx = solver_it.index();
      if(solver->isMobile()){
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velocities[idx]).getSingleValuedFAB()(iv, dir) = velocities[idx][dir];
	}
      }
    }
  }


  // Covered is always bogus.
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->isMobile()){
      const int idx = solver_it.index();
      for (int dir = 0; dir < SpaceDim; dir++){
	a_velocities[idx]->setCoveredCellVal(0.0, dir);
      }
    }
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocitiesRegularFast(Vector<EBCellFAB*>&       a_velocities,
							    const Vector<EBCellFAB*>& a_cdr_densities,
							    const EBCellFAB&          a_E,
							    const Box&                a_box,
							    const Real&               a_time,
							    const Real&               a_dx){
#if CH_SPACEDIM==2
  computeCdrDriftVelocitiesRegularFast2D(a_velocities, a_cdr_densities, a_E, a_box, a_time, a_dx);
#elif CH_SPACEDIM==3
  computeCdrDriftVelocitiesRegularFast3D(a_velocities, a_cdr_densities, a_E, a_box, a_time, a_dx);
#endif
}

void CdrPlasmaStepper::computeCdrDriftVelocitiesRegularFast2D(Vector<EBCellFAB*>&       a_velocities,
							      const Vector<EBCellFAB*>& a_cdr_densities,
							      const EBCellFAB&          a_E,
							      const Box&                a_box,
							      const Real&               a_time,
							      const Real&               a_dx){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocitiesRegularFast2D");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocitiesRegularFast2D" << endl;
  }
  
#if CH_SPACEDIM==2
  const int num_species  = m_physics->getNumCdrSpecies();
  const int comp         = 0;
  const RealVect origin  = m_amr->getProbLo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // I need contiguous memory for the nasty stuff that is about to happen, so begin by copying things onto smaller
  // data holders
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_vx(a_box,  num_species);
  FArrayBox cdr_vy(a_box,  num_species);
  for (int idx = 0; idx < num_species; idx++){
    cdr_phi.copy(a_cdr_densities[idx]->getFArrayBox(), a_box, 0, a_box, idx, 1);
  }

  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);


  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int offset     = lo[0] + n0*lo[1];

  // C style variable-length array conversion magic
  auto vla_cdr_phi = (Real (*__restrict__)[n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_cdr_vx  = (Real (*__restrict__)[n1][n0]) (cdr_vx.dataPtr()  - offset);
  auto vla_cdr_vy  = (Real (*__restrict__)[n1][n0]) (cdr_vy.dataPtr()  - offset);
  auto vla_E       = (Real (*__restrict__)[n1][n0]) (Efab.dataPtr()    - offset);
  
  for (int j = lo[1]; j <= hi[1]; ++j){
    for (int i = lo[0]; i <= hi[0]; ++i){

      const RealVect pos = origin + RealVect(D_DECL(i,j,k))*a_dx;
      const RealVect E   = RealVect(vla_E[0][j][i], vla_E[1][j][i]);

      Vector<Real> cdr_densities(num_species);
      for (int idx = 0; idx < num_species; idx++){
	cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][j][i]);
      }

      const Vector<RealVect> velocities = m_physics->computeCdrDriftVelocities(a_time, pos, E, cdr_densities);

      // Put result in correct palce
      for (int idx = 0; idx < num_species; ++idx){
	vla_cdr_vx[idx][j][i] = velocities[idx][0];
	vla_cdr_vy[idx][j][i] = velocities[idx][1];
      }
    }
  }

  // Put the results back into solvers
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      FArrayBox& v = a_velocities[idx]->getFArrayBox();

      v.copy(cdr_vx, a_box, idx, a_box, 0, 1);
      v.copy(cdr_vy, a_box, idx, a_box, 1, 1);

      for (int dir = 0; dir < SpaceDim; dir++){
	a_velocities[idx]->setCoveredCellVal(0.0, dir);
      }
    }
  }
#endif
}



void CdrPlasmaStepper::computeCdrDriftVelocitiesRegularFast3D(Vector<EBCellFAB*>&       a_velocities,
							      const Vector<EBCellFAB*>& a_cdr_densities,
							      const EBCellFAB&          a_E,
							      const Box&                a_box,
							      const Real&               a_time,
							      const Real&               a_dx){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocitiesRegularFast3D");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocitiesRegularFast3D" << endl;
  }

#if CH_SPACEDIM==3
  const int num_species  = m_physics->getNumCdrSpecies();
  const int comp         = 0;
  const RealVect origin  = m_amr->getProbLo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // I need contiguous memory for the nasty stuff that is about to happen, so begin by copying things onto smaller
  // data holders
  FArrayBox cdr_phi(a_box, num_species);
  FArrayBox cdr_vx(a_box,  num_species);
  FArrayBox cdr_vy(a_box,  num_species);
  FArrayBox cdr_vz(a_box,  num_species);
  for (int idx = 0; idx < num_species; idx++){
    cdr_phi.copy(a_cdr_densities[idx]->getFArrayBox(), a_box, 0, a_box, idx, 1);
  }

  FArrayBox Efab(a_box, SpaceDim);
  Efab.copy(a_E.getFArrayBox(), a_box, 0, a_box, 0, SpaceDim);


  // Pointer offsets
  const IntVect dims   = a_box.size();
  const IntVect lo     = a_box.smallEnd();
  const IntVect hi     = a_box.bigEnd();
  const int n0         = dims[0];
  const int n1         = dims[1];
  const int n2         = dims[2];
  const int offset     = lo[0] + n0*(lo[1] + n1*lo[2]);

  auto vla_cdr_phi = (Real (*__restrict__)[n2][n1][n0]) (cdr_phi.dataPtr() - offset);
  auto vla_cdr_vx  = (Real (*__restrict__)[n2][n1][n0]) (cdr_vx.dataPtr()  - offset);
  auto vla_cdr_vy  = (Real (*__restrict__)[n2][n1][n0]) (cdr_vy.dataPtr()  - offset);
  auto vla_cdr_vz  = (Real (*__restrict__)[n2][n1][n0]) (cdr_vz.dataPtr()  - offset);
  auto vla_E       = (Real (*__restrict__)[n2][n1][n0]) (Efab.dataPtr()    - offset);
  
  for (int k = lo[2]; k <= hi[2]; k++){
    for (int j = lo[1]; j <= hi[1]; ++j){
      for (int i = lo[0]; i <= hi[0]; ++i){

	const RealVect pos = origin + RealVect(D_DECL(i,j,k))*a_dx;
	const RealVect E   = RealVect(vla_E[0][k][j][i], vla_E[1][k][j][i], vla_E[2][k][j][i]);

	Vector<Real> cdr_densities(num_species);
	for (int idx = 0; idx < num_species; idx++){
	  cdr_densities[idx] = Max(0.0, vla_cdr_phi[idx][k][j][i]);
	}

	const Vector<RealVect> velocities = m_physics->computeCdrDriftVelocities(a_time, pos, E, cdr_densities);

	// Put result in correct palce
	for (int idx = 0; idx < num_species; ++idx){
	  vla_cdr_vx[idx][k][j][i] = velocities[idx][0];
	  vla_cdr_vy[idx][k][j][i] = velocities[idx][1];
	  vla_cdr_vz[idx][k][j][i] = velocities[idx][2];
	}
      }
    }
  }
  
  // Put the results back into solvers
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      FArrayBox& v = a_velocities[idx]->getFArrayBox();
      
      v.copy(cdr_vx, a_box, 0, a_box, 0, 1);
      v.copy(cdr_vy, a_box, 0, a_box, 1, 1);
      v.copy(cdr_vz, a_box, 0, a_box, 2, 1);

      for (int dir = 0; dir < SpaceDim; dir++){
	a_velocities[idx]->setCoveredCellVal(0.0, dir);
      }      
    }
  }
#endif
}


void CdrPlasmaStepper::computeCdrDriftVelocitiesIrregular(Vector<EBCellFAB*>&       a_velocities,
							  const Vector<EBCellFAB*>& a_cdr_densities,
							  const EBCellFAB&          a_E,
							  const Box&                a_box,
							  const Real&               a_time,
							  const Real&               a_dx,
							  const int                 a_lvl,
							  const DataIndex&          a_dit){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocitiesIrregular");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocitiesIrregular" << endl;
  }

  const int comp         = 0;
  const EBISBox& ebisbox = a_E.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const RealVect origin  = m_amr->getProbLo();

  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, origin);

    // Get densities
    Vector<Real> cdr_densities;
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      cdr_densities.push_back((*a_cdr_densities[idx])(vof, comp));
    }
    
    // Compute velocities
    const Vector<RealVect> velocities = m_physics->computeCdrDriftVelocities(a_time, pos, e, cdr_densities);

    // Put velocities in the appropriate place. 
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      if(solver_it()->isMobile()){
	const int idx = solver_it.index();
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velocities[idx])(vof, dir) = velocities[idx][dir];
	}
      }
    }
  }
}

void CdrPlasmaStepper::preRegrid(const int a_lmin, const int a_finestLevel){
  CH_TIME("CdrPlasmaStepper::preRegrid");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::preRegrid" << endl;
  }

  // Solvers do pre-regridding shit. 
  m_cdr->preRegrid(a_lmin, a_finestLevel);
  m_fieldSolver->preRegrid(a_lmin, a_finestLevel);
  m_rte->preRegrid(a_lmin, a_finestLevel);
  m_sigma->preRegrid(a_lmin, a_finestLevel);
}

void CdrPlasmaStepper::preRegridInternals(const int a_lbase, const int a_finestLevel){
  CH_TIME("CdrPlasmaStepper::preRegridInternals");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::preRegridInternals" << endl;
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocities(){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocities()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocities()" << endl;
  }

  // Compute the electric field (again)
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, m_cdr->getPhase(), SpaceDim);
  this->computeElectricField(E, m_cdr->getPhase(), m_fieldSolver->getPotential());

  Vector<EBAMRCellData*> states     = m_cdr->getPhis();
  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();

  this->computeCdrDriftVelocities(velocities, states, E, m_time);
}

void CdrPlasmaStepper::computeCdrDiffusion(){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusion");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusion" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_physics->getNumCdrSpecies();

  EBAMRCellData E_cell;
  EBAMRIVData   E_eb;
  m_amr->allocate(E_cell, m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(E_eb,   m_realm, m_cdr->getPhase(), SpaceDim);
  
  this->computeElectricField(E_cell, m_cdr->getPhase(), m_fieldSolver->getPotential());
  this->computeElectricField(E_eb,   m_cdr->getPhase(), E_cell);

  Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();

  CdrPlasmaStepper::computeCdrDiffusion(E_cell, E_eb);
#if 0

  // Extrapolate states to the EB
  Vector<EBAMRIVData*> cdr_extrap(num_species);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    cdr_extrap[idx] = new EBAMRIVData();  // This must be deleted
    m_amr->allocate(*cdr_extrap[idx], m_realm, m_cdr->getPhase(), ncomp);

    const IrregAmrStencil<EbCentroidInterpolationStencil>& stencil = m_amr->getEbCentroidInterpolationStencils(m_realm, m_cdr->getPhase());
    stencil.apply(*cdr_extrap[idx], *cdr_states[idx]);
  }
  
  Vector<EBAMRFluxData*> diffco_face = m_cdr->getFaceCenteredDiffusionCoefficient();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->getEbCenteredDiffusionCoefficient();

  this->computeCdrDiffusionFace(diffco_face, cdr_states, E_cell, m_time);
  this->computeCdrDiffusionEb(diffco_eb,     cdr_extrap, E_eb,   m_time);

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    delete cdr_extrap[idx];
  }
#endif
}

void CdrPlasmaStepper::computeCdrDiffusion(const EBAMRCellData& a_E_cell, const EBAMRIVData& a_E_eb){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusion");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusion" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_physics->getNumCdrSpecies();

  Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();

  // Extrapolate states to the EB
  Vector<EBAMRIVData*> cdr_extrap(num_species);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    cdr_extrap[idx] = new EBAMRIVData();  // This must be deleted
    m_amr->allocate(*cdr_extrap[idx], m_realm, m_cdr->getPhase(), ncomp);

    const IrregAmrStencil<EbCentroidInterpolationStencil>& stencil = m_amr->getEbCentroidInterpolationStencils(m_realm, m_cdr->getPhase());
    stencil.apply(*cdr_extrap[idx], *cdr_states[idx]);
  }
  
  Vector<EBAMRFluxData*> diffco_face = m_cdr->getFaceCenteredDiffusionCoefficient();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->getEbCenteredDiffusionCoefficient();

  this->computeCdrDiffusionFace(diffco_face, cdr_states, a_E_cell, m_time);
  this->computeCdrDiffusionEb(diffco_eb,     cdr_extrap, a_E_eb,   m_time);

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_amr->deallocate(*cdr_extrap[idx]);
    delete cdr_extrap[idx];
  }
}

void CdrPlasmaStepper::computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(mfamrcell, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(mfamrcell, mfamrcell)" << endl;
  }

  m_fieldSolver->computeElectricField(a_E, a_potential);

  m_amr->averageDown(a_E, m_realm);
  m_amr->interpGhostMG(a_E, m_realm);
}

void CdrPlasmaStepper::computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(ebamrcell, phase)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(ebamrcell, phase)" << endl;
  }

  this->computeElectricField(a_E, a_phase, m_fieldSolver->getPotential());
}

void CdrPlasmaStepper::computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(ebamrcell, phase, mfamrcell)" << endl;
  }

  m_fieldSolver->computeElectricField(a_E, a_phase, a_potential);

  m_amr->averageDown(a_E, m_realm, a_phase);
  m_amr->interpGhostMG(a_E, m_realm, a_phase);
}

void CdrPlasmaStepper::computeElectricField(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(ebamrflux, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(ebamrflux, phase, ebamrcell)" << endl;
  }

  CH_assert(a_E_face[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, a_phase)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& E_cell = (*a_E_cell[lvl])[dit()];
      const EBISBox& ebisbox  = ebisl[dit()];
      const EBGraph& ebgraph  = ebisbox.getEBGraph();
      const Box& box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& E_face = (*a_E_face[lvl])[dit()][dir];
	E_face.setVal(0.0);

	EBLevelDataOps::averageCellToFace(E_face,
					  E_cell,
					  ebgraph,
					  box,
					  0,
					  dir,
					  domain,
					  dir,
					  dir);
      }

    }
    //    DataOps::averageCellToFace(*a_E_face[lvl], *a_E_cell[lvl], m_amr->getDomains()[lvl]);
    a_E_face[lvl]->exchange();
  }
}

void CdrPlasmaStepper::computeElectricField(EBAMRIVData& a_E_eb, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(ebamriv, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(ebamriv, phase, ebamrcell)" << endl;
  }

  CH_assert(a_E_eb[0]->nComp()   == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const IrregAmrStencil<EbCentroidInterpolationStencil>& interp_stencil = m_amr->getEbCentroidInterpolationStencils(m_realm, a_phase);
  interp_stencil.apply(a_E_eb, a_E_cell);
}

void CdrPlasmaStepper::computeMaxElectricField(Real& a_Emax, const phase::which_phase a_phase){
  CH_TIME("CdrPlasmaStepper::computeMaxElectricField(Real, phase)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeMaxElectricField(Real, phase)" << endl;
  }

  EBAMRCellData E;
  m_amr->allocate(E, m_realm, a_phase, SpaceDim);

  this->computeElectricField(E, a_phase, m_fieldSolver->getPotential());
  m_amr->interpToCentroids(E, m_realm, a_phase);

  Real max, min;
  DataOps::getMaxMinNorm(max, min, E);

  a_Emax = max;
}

void CdrPlasmaStepper::computeChargeFlux(EBAMRIVData& a_flux, Vector<EBAMRIVData*>& a_cdr_fluxes){
  CH_TIME("CdrPlasmaStepper::computeChargeFlux");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeChargeFlux" << endl;
  }

  MayDay::Abort("CdrPlasmaStepper::computeChargeFlux - I'm suspecting that this is deprecated code which is no longer used");

  DataOps::setValue(a_flux, 0.0);

  for (CdrIterator<CdrSolver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = *a_cdr_fluxes[solver_it.index()];

    DataOps::incr(a_flux, solver_flux, spec->getChargeNumber()*Units::Qe);
  }

  m_sigma->resetCells(a_flux);
}

void CdrPlasmaStepper::computeExtrapolatedFluxes(Vector<EBAMRIVData*>&        a_fluxes,
						 const Vector<EBAMRCellData*> a_densities,
						 const Vector<EBAMRCellData*> a_velocities,
						 const phase::which_phase     a_phase){
  CH_TIME("CdrPlasmaStepper::computeExtrapolatedFluxes");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeExtrapolatedFluxes" << endl;
  }

#if 0 // This code computes the cell-centered flux which is then extrapolated to the EB. Please don't remove this code. 
  EBAMRCellData cell_flux;
  EBAMRIVData   eb_flux;

  m_amr->allocate(cell_flux, m_realm, a_phase, SpaceDim);
  m_amr->allocate(eb_flux,   m_realm, a_phase, SpaceDim);

  for (int i = 0; i < a_fluxes.size(); i++){
    this->computeFlux(cell_flux, *a_densities[i], *a_velocities[i]);
    
    m_amr->averageDown(cell_flux, m_realm, a_phase);
    m_amr->interpGhostMG(cell_flux, m_realm, a_phase);

    this->extrapolateToEb(eb_flux, a_phase, cell_flux);
    
    this->projectFlux(*a_fluxes[i], eb_flux);
  }
#else // New way of doing this. Extrapolate everything to the BC first. Then compute the flux and project it. We do this because
      // extrapolation stencils may have negative weights, and v*n may therefore nonphysically change sign. Better to compute
      // F = v_extrap*Max(0.0, phi_extrap) since we expect v to be "smooth" and phi_extrap to be a noisy bastard
  EBAMRIVData eb_flx; 
  EBAMRIVData eb_vel;
  EBAMRIVData eb_phi;

  m_amr->allocate(eb_flx, m_realm, a_phase, SpaceDim);
  m_amr->allocate(eb_vel, m_realm, a_phase, SpaceDim);
  m_amr->allocate(eb_phi, m_realm, a_phase, 1);

  const IrregAmrStencil<EbCentroidInterpolationStencil>& interp_stencils = m_amr->getEbCentroidInterpolationStencils(m_realm, a_phase);

  //  for (int i = 0; i < a_fluxes.size(); i++){
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      interp_stencils.apply(eb_vel, *a_velocities[idx]);
      interp_stencils.apply(eb_phi, *a_densities[idx]);

      DataOps::floor(eb_phi, 0.0);

      DataOps::setValue(eb_flx, 0.0);
      DataOps::incr(eb_flx, eb_vel, 1.0);
      DataOps::multiplyScalar(eb_flx, eb_phi);

      this->projectFlux(*a_fluxes[idx], eb_flx);

      m_amr->averageDown(*a_fluxes[idx], m_realm, a_phase);
    }
    else{
      DataOps::setValue(*a_fluxes[idx], 0.0);
    }
  }
#endif
}

void CdrPlasmaStepper::computeExtrapolatedVelocities(Vector<EBAMRIVData*>&        a_ebvelo,
						     const Vector<EBAMRCellData*> a_cell_vel,
						     const phase::which_phase     a_phase){
  CH_TIME("CdrPlasmaStepper::computeExtrapolatedVelocities");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeExtrapolatedVelocities" << endl;
  }

  EBAMRIVData scratch;
  m_amr->allocate(scratch, m_realm, a_phase, SpaceDim);

  //  for (int i = 0; i < a_ebvelo.size(); i++){
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      this->extrapolateToEb(scratch, a_phase, *a_cell_vel[idx]);
      this->projectFlux(*a_ebvelo[idx], scratch);
    }
  }
}

void CdrPlasmaStepper::computeExtrapolatedDomainFluxes(Vector<EBAMRIFData*>&        a_fluxes,
						       const Vector<EBAMRCellData*> a_densities,
						       const Vector<EBAMRCellData*> a_velocities,
						       const phase::which_phase     a_phase){
  CH_TIME("CdrPlasmaStepper::computeExtrapolatedDomainFluxes");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeExtrapolatedDomainFluxes" << endl;
  }

  EBAMRCellData cell_flux;
  EBAMRIFData domain_flux;

  m_amr->allocate(cell_flux,   m_realm, a_phase, SpaceDim);
  m_amr->allocate(domain_flux, m_realm, a_phase, SpaceDim);

  //  for (int i = 0; i < a_fluxes.size(); i++){
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      DataOps::copy(cell_flux, *a_velocities[idx]);
      DataOps::multiplyScalar(cell_flux, *a_densities[idx]); // cell_flux = n*v

      // Extrapolate cell-centered to domain faces
      this->extrapolateToDomainFaces(domain_flux, a_phase, cell_flux);

      // Project normal component onto domain face
      this->projectDomain(*a_fluxes[idx], domain_flux);
    }
    else{
      DataOps::setValue(*a_fluxes[idx], 0.0);
    }
  }
}

void CdrPlasmaStepper::computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_density, const EBAMRCellData& a_velocity){
  CH_TIME("CdrPlasmaStepper::computeFlux(amr)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeFlux(amr)" << endl;
  }
  
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    CH_assert(a_flux[lvl]->nComp()     == SpaceDim);
    CH_assert(a_density[lvl]->nComp()  == 1);
    CH_assert(a_velocity[lvl]->nComp() == SpaceDim);

    // Call level versions
    computeFlux(*a_flux[lvl], *a_density[lvl], *a_velocity[lvl], lvl);
  }
}

void CdrPlasmaStepper::computeFlux(LevelData<EBCellFAB>&       a_flux,
				   const LevelData<EBCellFAB>& a_density,
				   const LevelData<EBCellFAB>& a_velocity,
				   const int                   a_lvl){
  CH_TIME("CdrPlasmaStepper::computeFlux(amr)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeFlux(amr)" << endl;
  }

  DataOps::setValue(a_flux, 0.0);
  DataOps::incr(a_flux, a_velocity, 1.0);
  DataOps::multiplyScalar(a_flux, a_density);
}

void CdrPlasmaStepper::computeJ(EBAMRCellData& a_J) const{
  CH_TIME("CdrPlasmaStepper::computeJ(amr)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeJ(amr)" << endl;
  }
  
  const int finest_level = m_amr->getFinestLevel();

  // Call level versions
  for (int lvl = 0; lvl <= finest_level; lvl++){
    computeJ(*a_J[lvl], lvl);
  }

  m_amr->averageDown(a_J, m_realm, m_cdr->getPhase());
  m_amr->interpGhostMG(a_J, m_realm, m_cdr->getPhase());
}

void CdrPlasmaStepper::computeJ(LevelData<EBCellFAB>& a_J, const int a_lvl) const{
  CH_TIME("CdrPlasmaStepper::computeJ(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeJ(level)" << endl;
  }

  DataOps::setValue(a_J, 0.0);
  
  const int density_comp = 0;

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();

    if(solver->isMobile()){
      const int q                       = spec->getChargeNumber();
      const EBAMRCellData& density      = solver->getPhi();
      const EBAMRCellData& velo         = solver->getCellCenteredVelocity();

      if(q != 0){
	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  const Box& box         = dbl.get(dit());
	  const EBISBox& ebisbox = ebisl[dit()];
	  const EBGraph& ebgraph = ebisbox.getEBGraph();
	  const IntVectSet ivs(box);
      
	  EBCellFAB& J       = a_J[dit()];
	  const EBCellFAB& n = (*density[a_lvl])[dit()];
	  const EBCellFAB& v = (*velo[a_lvl])[dit()];
      
	  EBCellFAB cdr_j(ebisbox, box, SpaceDim);
	  cdr_j.setVal(0.0);
	  for (int comp = 0; comp < SpaceDim; comp++){
	    cdr_j.plus(n, 0, comp, 1);
	  }
	  cdr_j *= v;
	  cdr_j *= q;
      
	  J += cdr_j;

	  // Should we monkey with irregular cells??? 
	}
      }
    }
  }
  
  DataOps::scale(a_J, Units::Qe);
}



void CdrPlasmaStepper::deallocateSolverInternals(){
  CH_TIME("CdrPlasmaStepper::deallocateSolverInternals");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::deallocateSolverInternals" << endl;
  }

  m_cdr->deallocateInternals();
  m_rte->deallocateInternals();
  m_fieldSolver->deallocateInternals();
  m_sigma->deallocateInternals();
}

void CdrPlasmaStepper::extrapolateToEb(Vector<EBAMRIVData*>&         a_extrap,
				       const phase::which_phase      a_phase,
				       const Vector<EBAMRCellData*>& a_data){
  CH_TIME("CdrPlasmaStepper::extrapolateToEb(vec)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToEb(vec)" << endl;
  }

  CH_assert(a_extrap.size() == a_data.size());

  for (int i = 0; i < a_extrap.size(); i++){
    this->extrapolateToEb(*a_extrap[i], a_phase, *a_data[i]);
  }
}

void CdrPlasmaStepper::extrapolateToEb(EBAMRIVData& a_extrap, const phase::which_phase a_phase, const EBAMRCellData& a_data){
  CH_TIME("CdrPlasmaStepper::extrapolateToEb");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToEb" << endl;
  }

  const IrregAmrStencil<EbCentroidInterpolationStencil>& stencils = m_amr->getEbCentroidInterpolationStencils(m_realm, a_phase);
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    extrapolateToEb(*a_extrap[lvl], a_phase, *a_data[lvl], lvl);
  }
}

void CdrPlasmaStepper::extrapolateToEb(LevelData<BaseIVFAB<Real> >& a_extrap,
				       const phase::which_phase     a_phase,
				       const LevelData<EBCellFAB>&  a_data,
				       const int                    a_lvl){
  CH_TIME("CdrPlasmaStepper::extrapolateToEb(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToEb(level)" << endl;
  }

  const IrregAmrStencil<EbCentroidInterpolationStencil>& stencils = m_amr->getEbCentroidInterpolationStencils(m_realm, a_phase);
  stencils.apply(a_extrap, a_data, a_lvl);
}

void CdrPlasmaStepper::extrapolateToDomainFaces(EBAMRIFData&             a_extrap,
						const phase::which_phase a_phase,
						const EBAMRCellData&     a_data){
  CH_TIME("CdrPlasmaStepper::extrapolateToDomainFaces(amr)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToDomainFaces(amr)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    extrapolateToDomainFaces(*a_extrap[lvl], a_phase, *a_data[lvl], lvl);
  }
}

void CdrPlasmaStepper::extrapolateToDomainFaces(LevelData<DomainFluxIFFAB>& a_extrap,
						const phase::which_phase    a_phase,
						const LevelData<EBCellFAB>& a_data,
						const int                   a_lvl){
  CH_TIME("CdrPlasmaStepper::extrapolateToDomainFaces(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToDomainFaces(level)" << endl;
  }

  const int ncomp = a_data.nComp();
      
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, a_phase)[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBCellFAB& data         = a_data[dit()];
    const EBISBox& ebisbox        = ebisl[dit()];
    const BaseFab<Real>& data_fab = data.getSingleValuedFAB();
      
    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	BaseIFFAB<Real>& extrap = a_extrap[dit()](dir, sit());

	const IntVectSet& ivs  = extrap.getIVS();
	const EBGraph& ebgraph = extrap.getEBGraph();

	// Extrapolate to the boundary. Use face-centered stuff for all faces (also multivalued ones)
	const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
	for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  const int sgn = sign(sit()); // Lo = -1, Hi = 1
	    
	  const VolIndex& vof = face.getVoF(flip(sit()));
	  const IntVect iv0   = vof.gridIndex();
	  const IntVect iv1   = iv0 - sgn*BASISV(dir);

	  if(ebisbox.isCovered(iv0)){ // Just provide some bogus data because the face 
	    for (int comp = 0; comp < ncomp; comp++){
	      extrap(face, comp) = 0.0;
	    }
	  }
	  else{
	    if(!ebisbox.isCovered(iv1)){ // linear extrapolation
	      for (int comp = 0; comp < ncomp; comp++){
		extrap(face, comp) = 1.5*data_fab(iv0, comp) - 0.5*data_fab(iv1, comp); // Should be ok
	      }
	    }
	    else{ // Not enough cells available, use cell-centered only
	      for (int comp = 0; comp < ncomp; comp++){
		extrap(face, comp) = data_fab(iv0, comp);
	      }
	    }
	  }
	}
      }
    }
  }
}

void CdrPlasmaStepper::extrapolateToDomainFaces(Vector<EBAMRIFData*>&         a_extrap,
						const phase::which_phase      a_phase,
						const Vector<EBAMRCellData*>& a_data){
  CH_TIME("CdrPlasmaStepper::extrapolateToDomainFaces(vec)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::extrapolateToDomainFaces(vec)" << endl;
  }

  CH_assert(a_extrap.size() == a_data.size());

  for (int i = 0; i < a_extrap.size(); i++){
    this->extrapolateToDomainFaces(*a_extrap[i], a_phase, *a_data[i]);
  }
}

void CdrPlasmaStepper::getCdrMax(Real& a_cdr_max, std::string& a_solver_name){
  CH_TIME("CdrPlasmaStepper::getCdrMax");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::getCdrMax" << endl;
  }

  const int comp = 0;
  
  a_cdr_max = -1.E99;
  a_solver_name = "invalid solver";
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    Real max, min;
    DataOps::getMaxMin(max, min, solver->getPhi(), comp);

    if(max > a_cdr_max){
      a_cdr_max     = max;
      a_solver_name = solver->getName();
    }
  }
}

void CdrPlasmaStepper::setCdrSolvers(RefCountedPtr<CdrLayout<CdrSolver>>& a_cdr){
  CH_TIME("CdrPlasmaStepper::setCdrSolvers");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setCdrSolvers" << endl;
  }
  m_cdr = a_cdr;
}

void CdrPlasmaStepper::setFieldSolver(RefCountedPtr<FieldSolver>& a_poisson){
  CH_TIME("CdrPlasmaStepper::setFieldSolver");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setFieldSolver" << endl;
  }
  m_fieldSolver = a_poisson;
}

void CdrPlasmaStepper::setRadiativeTransferSolvers(RefCountedPtr<RtLayout<RtSolver>>& a_rte){
  CH_TIME("CdrPlasmaStepper::setRadiativeTransferSolvers");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setRadiativeTransferSolvers" << endl;
  }
  m_rte = a_rte;
}

void CdrPlasmaStepper::setupSolvers(){
  CH_TIME("CdrPlasmaStepper::setupSolvers");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setupSolvers" << endl;
  }
  
  parseOptions();

  // Make solvers
  this->setupCdr();
  this->setupRadiativeTransfer(); 
  this->setupPoisson();
  this->setupSigma();

  this->setSolverVerbosity();

  // Allocate internal memory
  this->allocateInternals();
}

void CdrPlasmaStepper::initialData(){
  CH_TIME("CdrPlasmaStepper::initialData");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::initialData" << endl;
  }


  m_cdr->initialData();               // Initial data comes in through CdrSpecies, in this case supplied by physics
  m_fieldSolver->setPermittivities(); // Set permittivities for Poisson operator
  if(!m_rte->isStationary()){
    m_rte->initialData();
  }
  this->initialSigma();

  // Solve Poisson equation
  this->solvePoisson();

  // Fill solvers with velocity and diffusion
  this->computeCdrDriftVelocities();
  this->computeCdrDiffusion();

  // Do stationary RTE solve if we must
  if(this->stationaryRTE()){                  // Solve RTE equations by using initial data and electric field
    const Real dummy_dt = 1.0;

    this->solveRadiativeTransfer(dummy_dt);                 // Argument does not matter, it's a stationary solver.
  }
}

void CdrPlasmaStepper::initialSigma(){
  CH_TIME("CdrPlasmaStepper::initialSigma");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::initialSigma" << endl;
  }

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  EBAMRIVData& sigma = m_sigma->getPhi();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::gas)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& state = (*sigma[lvl])[dit()];

      const EBISBox& ebisbox = ebisl[dit()];
      const IntVectSet& ivs  = state.getIVS();
      const EBGraph& ebgraph = state.getEBGraph();
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = origin + vof.gridIndex()*dx + ebisbox.bndryCentroid(vof)*dx;
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_physics->initialSigma(m_time, pos);
	}
      }
    }
  }

  m_amr->averageDown(sigma, m_realm, phase::gas);
  m_sigma->resetCells(sigma);
}

void CdrPlasmaStepper::projectFlux(LevelData<BaseIVFAB<Real> >&       a_projected_flux,
				   const LevelData<BaseIVFAB<Real> >& a_flux,
				   const int                          a_lvl){
  CH_TIME("CdrPlasmaStepper::projectFlux(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::projectFlux(level)" << endl;
  }

  CH_assert(a_projected_flux.nComp() == 1);
  CH_assert(a_flux.nComp()           == SpaceDim);

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box         = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    
    BaseIVFAB<Real>& proj_flux  = a_projected_flux[dit()];
    const BaseIVFAB<Real>& flux = a_flux[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof     = vofit();
      const RealVect& normal  = ebisbox.normal(vof);
      const RealVect vec_flux = RealVect(D_DECL(flux(vof,0), flux(vof,1), flux(vof,2)));
      
      // For EB's, the geometrical normal vector is opposite of the finite volume method normal
      proj_flux(vof,0) = PolyGeom::dot(vec_flux, -normal);
    }
  }
}

void CdrPlasmaStepper::projectFlux(EBAMRIVData& a_projected_flux, const EBAMRIVData& a_flux){
  CH_TIME("CdrPlasmaStepper::projectFlux");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::projectFlux" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    CdrPlasmaStepper::projectFlux(*a_projected_flux[lvl], *a_flux[lvl], lvl);
  }
}

void CdrPlasmaStepper::projectDomain(EBAMRIFData& a_projected_flux, const EBAMRIFData& a_flux){
  CH_TIME("CdrPlasmaStepper::projectDomain");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::projectDomain" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    CH_assert(a_projected_flux[lvl]->nComp() == 1);
    CH_assert(a_flux[lvl]->nComp()           == SpaceDim);

    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
      for (int dir = 0; dir < SpaceDim; dir++){
	for (SideIterator sit; sit.ok(); ++sit){
	  BaseIFFAB<Real>& normal_comp    = (*a_projected_flux[lvl])[dit()](dir, sit());
	  const BaseIFFAB<Real>& gradient = (*a_flux[lvl])[dit()](dir, sit());

	  const int sgn          = sign(sit());
	  const EBGraph& ebgraph = normal_comp.getEBGraph();
	  const IntVectSet& ivs  = normal_comp.getIVS();

	  for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	    const FaceIndex& face = faceit();
	    normal_comp(face, comp) = gradient(face, dir);
	  }
	}
      }
    }
  }
}

void CdrPlasmaStepper::regrid(const int a_lmin, const int a_old_finest, const int a_new_finest){
  CH_TIME("CdrPlasmaStepper::regrid");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::regrid" << endl;
  }

  this->allocateInternals(); // Allocate memory for time stepper
  this->regridSolvers(a_lmin, a_old_finest, a_new_finest);
  this->regridInternals(a_lmin, a_old_finest, a_new_finest);

  // Solvers have been regridded. Now resolve the Poisson equation with the new data
  bool converged = this->solvePoisson();

  // If we don't converge, try new Poisson solver settings
  if(!converged){ 
    if(m_verbosity > 0){
      pout() << "CdrPlasmaStepper::regrid - Poisson solver failed to converge." << endl;
    }
  }

  // Compute stuff that is important for the CDR solvers. 
  this->computeCdrDriftVelocities();
  this->computeCdrDiffusion();

  // If we're doing a stationary RTE solve, recompute source terms
  if(this->stationaryRTE()){     // Solve RTE equations by using data that exists inside solvers
    const Real dummy_dt = 1.0;

    // Need new source terms for RTE equations
    this->advanceReactionNetwork(m_time, dummy_dt);
    this->solveRadiativeTransfer(dummy_dt);    // Argument does not matter, it's a stationary solver.
  }
}

void CdrPlasmaStepper::regridSolvers(const int a_lmin, const int a_old_finest, const int a_new_finest){
  CH_TIME("CdrPlasmaStepper::regridSolvers");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::regridSolvers" << endl;
  }

  m_cdr->regrid(a_lmin,     a_old_finest, a_new_finest);
  m_fieldSolver->regrid(a_lmin, a_old_finest, a_new_finest);
  m_rte->regrid(a_lmin,     a_old_finest, a_new_finest);
  m_sigma->regrid(a_lmin,   a_old_finest, a_new_finest);
}

void CdrPlasmaStepper::resetDielectricCells(EBAMRIVData& a_data){
  CH_TIME("CdrPlasmaStepper::resetDielectricCells");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::resetDielectricCells" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const MFLevelGrid& mflg      = *m_amr->getMFLevelGrid(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      BaseIVFAB<Real>& data  = (*a_data[lvl])[dit()];
      const IntVectSet ivs   = data.getIVS() & mflg.interfaceRegion(box, dit());
      const EBGraph& ebgraph = data.getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	for (int comp = 0; comp < data.nComp(); comp++){
	  data(vof, comp) = 0.0;
	}
      }
    }
  }
}

void CdrPlasmaStepper::sanityCheck(){
  CH_TIME("CdrPlasmaStepper::sanityCheck");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::sanityCheck" << endl;
  }

  CH_assert(!m_computationalGeometry.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_physics.isNull());
}

void CdrPlasmaStepper::setCdrPlasmaPhysics(const RefCountedPtr<CdrPlasmaPhysics>& a_physics){
  CH_TIME("CdrPlasmaStepper::setCdrPlasmaPhysics");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setCdrPlasmaPhysics" << endl;
  }

  m_physics = a_physics;
}

void CdrPlasmaStepper::setVoltage(std::function<Real(const Real a_time)> a_potential){
  CH_TIME("CdrPlasmaStepper::setVoltage");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setVoltage" << endl;
  }

  m_potential = a_potential;
}

void CdrPlasmaStepper::parseVerbosity(){
  CH_TIME("CdrPlasmaStepper::parseVerbosity");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::parseVerbosity" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  pp.get("verbosity", m_verbosity);
}

void CdrPlasmaStepper::parseSolverVerbosity(){
  CH_TIME("CdrPlasmaStepper::parseSolverVerbosity");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::parseSolverVerbosity" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  pp.get("solver_verbosity", m_solverVerbosity);
}

void CdrPlasmaStepper::parseCFL(){

  ParmParse pp(m_className.c_str());
  pp.get("cfl", m_cfl);
  if(m_cfl < 0.0){
    MayDay::Abort("CdrPlasmaStepper::parseCFL - CFL cannot be negative!");
  }
}

void CdrPlasmaStepper::parseRelaxationTime(){
  ParmParse pp(m_className.c_str());
  pp.get("relax_time", m_relax_time);
  if(m_relax_time < 0.0){
    MayDay::Abort("CdrPlasmaStepper::parseRelaxationTime - relaxation time cannot be negative");
  }
}

void CdrPlasmaStepper::parseMinDt(){
  ParmParse pp(m_className.c_str());
  pp.get("min_dt", m_minDt);
  if(m_minDt < 0.0){
    MayDay::Abort("CdrPlasmaStepper::parseMinDt - value cannot be negative");
  }
}

void CdrPlasmaStepper::parseMaxDt(){
  ParmParse pp(m_className.c_str());
  pp.get("max_dt", m_maxDt);
  if(m_maxDt < 0.0){
    MayDay::Abort("CdrPlasmaStepper::parseMaxDt - value cannot be negative");
  }
}

void CdrPlasmaStepper::parseFastRadiativeTransfer(){
  ParmParse pp(m_className.c_str());
  pp.get("fast_rte", m_fast_rte);
  if(m_fast_rte <= 0){
    MayDay::Abort("CdrPlasmaStepper::parseFastRadiativeTransfer - value cannot be negative");
  }
}

void CdrPlasmaStepper::parseFastPoisson(){
  ParmParse pp(m_className.c_str());
  pp.get("fast_poisson", m_fast_poisson);
  if(m_fast_poisson <= 0){
    MayDay::Abort("CdrPlasmaStepper::parseFastPoisson - value cannot be negative");
  }
}

void CdrPlasmaStepper::parseSourceComputation(){
  CH_TIME("CdrPlasmaStepper::parseSourceComputation");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::parseSourceComputation" << endl;
  }

  // Parse the source term computation. We support:
  //    'Interpolated', which interpolates cell-centered data to centroids before reaction terms.
  //    'CellAverage', which uses centroid-centered data
  //    'Upwind', which uses the Villa et. al upwind methodology.

  ParmParse pp(m_className.c_str());

  std::string str;

  const int numArgs = pp.countval("source_comp");
  if(numArgs == 1){
    pp.get("source_comp", str);

    if(str == "interp"){
      m_whichSourceTermComputation = SourceTermComputation::Interpolated;
    }
    else if (str == "interp2"){
      m_whichSourceTermComputation = SourceTermComputation::InterpolatedStable;
    }
    else if(str == "cell_ave"){
      m_whichSourceTermComputation = SourceTermComputation::CellAverage;
    }
    else{
      MayDay::Abort("CdrPlasmaStepper::parseSourceComputation -- logic bust. Expected 'interp' or 'cell_ave'");
    }
  }
  else if(numArgs == 2){
    pp.get("source_comp", str,            0);
    pp.get("source_comp", m_upwindFactor, 1);

    if(str == "upwind"){
      m_whichSourceTermComputation = SourceTermComputation::Upwind;
    }
    else{
      MayDay::Error("CdrPlasmaStepper::parseSourceComputation -- logic bust. Expected 'upwind <integer>'");
    }
  }
  else{
    MayDay::Error("CdrPlasmaStepper::parseSourceComputation -- logic bust. Expected at least one parameter");
  }
}

void CdrPlasmaStepper::setSolverVerbosity(){
  CH_TIME("CdrPlasmaStepper::setSolverVerbosity");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setSolverVerbosity" << endl;
  }

  if(!m_cdr.isNull()){
    m_cdr->setVerbosity(m_solverVerbosity);
  }
  if(!m_fieldSolver.isNull()){
    m_fieldSolver->setVerbosity(m_solverVerbosity);
  }
  if(!m_rte.isNull()){
    m_rte->setVerbosity(m_solverVerbosity);
  }
  if(!m_sigma.isNull()){
    m_sigma->setVerbosity(m_solverVerbosity);
  }
}

void CdrPlasmaStepper::setFastRadiativeTransfer(const int a_fast_rte){
  m_fast_rte = a_fast_rte;

  ParmParse pp("CdrPlasmaStepper");
  pp.query("fast_rte", m_fast_rte);
  if(m_fast_rte <= 0){
    m_fast_rte = a_fast_rte;
  }
}

void CdrPlasmaStepper::setFastPoisson(const int a_fast_poisson){
  m_fast_poisson = a_fast_poisson;

  ParmParse pp("CdrPlasmaStepper");
  pp.query("fast_poisson", m_fast_poisson);
  if(m_fast_poisson <= 0){
    m_fast_poisson = a_fast_poisson;
  }
}

void CdrPlasmaStepper::setMinDt(const Real a_min_dt){
  const Real zero = 0.0;
  m_minDt = Max(a_min_dt, zero);

  ParmParse pp("CdrPlasmaStepper");
  pp.query("min_dt", m_minDt);
  if(m_minDt < 0.0){
    m_minDt = a_min_dt;
  }
}

void CdrPlasmaStepper::setMaxDt(const Real a_max_dt){
  const Real zero = 0.0;
  m_maxDt = Max(a_max_dt, zero);

  ParmParse pp("CdrPlasmaStepper");
  pp.query("max_dt", m_maxDt);
  if(m_maxDt < 0.0){
    m_maxDt = a_max_dt;
  }
}

void CdrPlasmaStepper::setCFL(const Real a_cfl){
  m_cfl = a_cfl;

  ParmParse pp("CdrPlasmaStepper");
  pp.query("cfl", m_cfl);
  if(m_cfl < 0.0){
    m_cfl = a_cfl;
  }
}

void CdrPlasmaStepper::setRelaxTime(const Real a_relax_time){
  m_relax_time = a_relax_time;

  ParmParse pp("CdrPlasmaStepper");
  pp.query("relax_time", m_relax_time);
  if(m_relax_time < 0.0){
    m_relax_time = a_relax_time;
  }
}

void CdrPlasmaStepper::setupCdr(){
  CH_TIME("CdrPlasmaStepper::setupCdr");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setupCdr" << endl;
  }

  m_cdr->setVerbosity(m_solverVerbosity);
  m_cdr->parseOptions();
  m_cdr->setAmr(m_amr);
  m_cdr->setComputationalGeometry(m_computationalGeometry);
  m_cdr->setPhase(phase::gas);
  m_cdr->setRealm(m_realm);
}

void CdrPlasmaStepper::allocate() {
  m_cdr->allocateInternals();
  m_fieldSolver->allocateInternals();
  m_rte->allocateInternals();
  m_sigma->allocateInternals();
}

void CdrPlasmaStepper::setupPoisson(){
  CH_TIME("CdrPlasmaStepper::setupPoisson");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setupPoisson" << endl;
  }

  m_fieldSolver->setVerbosity(m_solverVerbosity);
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setRealm(m_realm);
  m_fieldSolver->setVoltage(m_potential); // Needs to happen AFTER setFieldSolver_wall_func
}

void CdrPlasmaStepper::setupRadiativeTransfer(){
  CH_TIME("CdrPlasmaStepper::setupRadiativeTransfer");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setupRadiativeTransfer" << endl;
  }

  m_rte->setVerbosity(m_solverVerbosity);
  m_rte->parseOptions();
  m_rte->setPhase(phase::gas);
  m_rte->setAmr(m_amr);
  m_rte->setComputationalGeometry(m_computationalGeometry);
  m_rte->sanityCheck();
  m_rte->setRealm(m_realm);
  //  m_rte->allocateInternals();
}

void CdrPlasmaStepper::setupSigma(){
  CH_TIME("CdrPlasmaStepper::setupSigma");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::setupSigma" << endl;
  }

  m_sigma = RefCountedPtr<SigmaSolver> (new SigmaSolver());
  m_sigma->setAmr(m_amr);
  m_sigma->setVerbosity(m_solverVerbosity);
  m_sigma->setComputationalGeometry(m_computationalGeometry);
  m_sigma->setRealm(m_realm);
  //  m_sigma->allocateInternals();
}

void CdrPlasmaStepper::solverDump(){
  CH_TIME("CdrPlasmaStepper::solverDump");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::solverDump" << endl;
  }

  m_cdr->writePlotFile();
  m_fieldSolver->writePlotFile();
  m_rte->writePlotFile();
}

void CdrPlasmaStepper::solveRadiativeTransfer(const Real a_dt){
  CH_TIME("CdrPlasmaStepper::solveRadiativeTransfer()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::solveRadiativeTransfer()" << endl;
  }

  const phase::which_phase rte_phase = m_rte->getPhase();

  EBAMRCellData E;
  m_amr->allocate(E, m_realm, rte_phase, SpaceDim);
  this->computeElectricField(E, rte_phase, m_fieldSolver->getPotential());

  Vector<EBAMRCellData*> states     = m_rte->getPhis();
  Vector<EBAMRCellData*> rhs        = m_rte->getSources();
  Vector<EBAMRCellData*> cdr_states = m_cdr->getPhis();

  this->solveRadiativeTransfer(states, rhs, cdr_states, E, m_time, a_dt, Centering::CellCenter);
}

void CdrPlasmaStepper::solveRadiativeTransfer(Vector<EBAMRCellData*>&       a_rte_states,
					      Vector<EBAMRCellData*>&       a_rte_sources,
					      const Vector<EBAMRCellData*>& a_cdr_states,
					      const EBAMRCellData&          a_E,
					      const Real                    a_time,
					      const Real                    a_dt,
					      const Centering               a_centering){
  CH_TIME("CdrPlasmaStepper::solveRadiativeTransfer(full)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::solveRadiativeTransfer(full)" << endl;
  }

  //  this->compute_rte_sources(a_rte_sources, a_cdr_states, a_E, a_time, a_centering);
  //  advanceReactionNetwork(a_time, a_dt);

  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    
    RefCountedPtr<RtSolver>& solver = solver_it();
    EBAMRCellData& state              = *a_rte_states[idx];
    EBAMRCellData& rhs                = *a_rte_sources[idx];
    solver->advance(a_dt, state, rhs);
  }
}

void CdrPlasmaStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("CdrPlasmaStepper::synchronizeSolverTimes");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;

  m_cdr->setTime(a_step,     a_time, a_dt);
  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_rte->setTime(a_step,     a_time, a_dt);
  m_sigma->setTime(a_step,   a_time, a_dt);
}

Real CdrPlasmaStepper::computeElectrodeCurrent(){
  CH_TIME("CdrPlasmaStepper::computeElectrodeCurrent");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectrodeCurrent" << endl;
  }

  // Need to copy onto temporary storage because 
  EBAMRIVData charge_flux;
  m_amr->allocate(charge_flux, m_realm, m_cdr->getPhase(), 1);
  DataOps::setValue(charge_flux, 0.0);
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    DataOps::incr(charge_flux, solver_flux, spec->getChargeNumber()*Units::Qe);
  }

  this->resetDielectricCells(charge_flux);
  m_amr->conservativeAverage(charge_flux, m_realm, m_cdr->getPhase());

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->getDx()[compute_lvl];
  for (DataIterator dit = m_amr->getGrids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const BaseIVFAB<Real>& flx = (*charge_flux[compute_lvl])[dit()];

    const IntVectSet ivs = flx.getIVS() & m_amr->getGrids(m_realm)[compute_lvl].get(dit());
    for (VoFIterator vofit(ivs, flx.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();
      const Real& bndryFrac = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[compute_lvl][dit()].bndryArea(vof);
      const Real& flux = flx(vof, 0);
      sum += flux*bndryFrac;
    }
  }

  sum *= pow(dx, SpaceDim-1);


#ifdef CH_MPI
  const Real sum1 = sum;
  sum = EBLevelDataOps::parallelSum(sum1);  
#endif

  return sum;
}

Real CdrPlasmaStepper::computeDielectricCurrent(){
  CH_TIME("CdrPlasmaStepper::computeDielectricCurrent");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeDielectricCurrent" << endl;
  }

  // Need to copy onto temporary storage because 
  EBAMRIVData charge_flux;
  m_amr->allocate(charge_flux, m_realm, m_cdr->getPhase(), 1);
  DataOps::setValue(charge_flux, 0.0);
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    DataOps::incr(charge_flux, solver_flux, spec->getChargeNumber()*Units::Qe);
  }

  m_sigma->resetCells(charge_flux);
  m_amr->conservativeAverage(charge_flux, m_realm, m_cdr->getPhase());

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->getDx()[compute_lvl];
  for (DataIterator dit = m_amr->getGrids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const BaseIVFAB<Real>& flx = (*charge_flux[compute_lvl])[dit()];

    const IntVectSet ivs = flx.getIVS() & m_amr->getGrids(m_realm)[compute_lvl].get(dit());
    for (VoFIterator vofit(ivs, flx.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof   = vofit();
      const Real& bndryFrac = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[compute_lvl][dit()].bndryArea(vof);
      const Real& flux = flx(vof, 0);
      sum += flux*bndryFrac;
    }
  }

  sum *= pow(dx, SpaceDim-1);


#ifdef CH_MPI
  const Real sum1 = sum;
  sum = EBLevelDataOps::parallelSum(sum1);  
#endif

  return sum;
}

Real CdrPlasmaStepper::computeDomainCurrent(){
  CH_TIME("CdrPlasmaStepper::computeDomainCurrent");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeDomainCurrent" << endl;
  }

  const int comp = 0;

  // Need to copy onto temporary storage because 
  EBAMRIFData charge_flux;
  m_amr->allocate(charge_flux, m_realm, m_cdr->getPhase(), 1);
  DataOps::setValue(charge_flux, 0.0);
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec      = solver_it.getSpecies();
    const EBAMRIFData& solver_flux          = solver->getDomainFlux();

    DataOps::incr(charge_flux, solver_flux, spec->getChargeNumber()*Units::Qe);
  }

  const int compute_lvl = 0;
  Real sum = 0.0;
  const Real dx = m_amr->getDx()[compute_lvl];
  for (DataIterator dit = m_amr->getGrids(m_realm)[compute_lvl].dataIterator(); dit.ok(); ++dit){
    const DomainFluxIFFAB& flux = (*charge_flux[compute_lvl])[dit()];

    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	const BaseIFFAB<Real>& fluxdir = flux(dir, sit());

	FaceStop::WhichFaces stopcrit = FaceStop::AllBoundaryOnly;
	const IntVectSet& ivs  = fluxdir.getIVS();
	const EBGraph& ebgraph = fluxdir.getEBGraph();

	for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  sum += sign(sit())*fluxdir(faceit(), comp);
	}
      }
    }
  }

  sum *= pow(dx, SpaceDim-1);


#ifdef CH_MPI
  const Real sum1 = sum;
  sum = EBLevelDataOps::parallelSum(sum1);  
#endif

  return sum;
}

Real CdrPlasmaStepper::computeOhmicInductionCurrent(){
  CH_TIME("CdrPlasmaStepper::computeRelaxationTime");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeRelaxationTime" << endl;
  }

  Real current = 0.0;

  EBAMRCellData J, E, JdotE;
  m_amr->allocate(J,     m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(E,     m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(JdotE, m_realm, m_cdr->getPhase(), SpaceDim);

  this->computeElectricField(E, m_cdr->getPhase(), m_fieldSolver->getPotential());
  this->computeJ(J);

  // Compute J.dot.E 
  DataOps::dotProduct(JdotE, J,E);
  m_amr->averageDown(JdotE, m_realm, m_cdr->getPhase());

  // Only compue on coarsest level
  const int coar = 0;
  const Real dx  = m_amr->getDx()[coar];
  DataOps::kappaSum(current, *JdotE[coar]);
  current *= pow(dx, SpaceDim);

  return current;
}

Real CdrPlasmaStepper::computeRelaxationTime(){
  CH_TIME("CdrPlasmaStepper::computeRelaxationTime");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeRelaxationTime" << endl;
  }

  // TLDR: This computes the relaxation time as t = eps0/conductivity. Simple as that. 

  EBAMRCellData relaxTime;
  EBAMRCellData conductivity;

  m_amr->allocate(relaxTime,    m_realm, phase::gas, 1);
  m_amr->allocate(conductivity, m_realm, phase::gas, 1);

  this->computeCellConductivity(conductivity);

  m_amr->averageDown  (conductivity, m_realm, phase::gas);
  m_amr->interpGhostMG(conductivity, m_realm, phase::gas);
  
  m_amr->interpToCentroids(conductivity, m_realm, phase::gas);  

  DataOps::setValue(relaxTime, Units::eps0);
  DataOps::divideByScalar(relaxTime, conductivity);

  Real maxVal;
  Real minVal;

  DataOps::getMaxMinNorm(maxVal, minVal, relaxTime);
  
  return minVal;
}

Real CdrPlasmaStepper::getTime() const {
  return m_time;
}

Real CdrPlasmaStepper::getDt(){
  return m_dt;
}

Real CdrPlasmaStepper::getCflDt(){
  return m_dt_cfl;
}

RefCountedPtr<CdrLayout<CdrSolver>>& CdrPlasmaStepper::getCdrSolvers(){
  return m_cdr;
}

RefCountedPtr<FieldSolver>& CdrPlasmaStepper::getFieldSolver(){
  return m_fieldSolver;
}

RefCountedPtr<RtLayout<RtSolver>>& CdrPlasmaStepper::getRadiativeTransferSolvers(){
  return m_rte;
}

RefCountedPtr<SigmaSolver>& CdrPlasmaStepper::getSigmaSolver(){
  return m_sigma;
}

void CdrPlasmaStepper::deallocate() {
  CH_TIME("CdrPlasmaStepper::deallocate");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::deallocate" << endl;
  }

  this->deallocateInternals();
  this->deallocateSolverInternals();  
}

#ifdef CH_USE_HDF5
void CdrPlasmaStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("CdrPlasmaStepper::writeCheckpointData");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::writeCheckpointData" << endl;
  }

  // CDR solvers checkpoint their data
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  // RTE solvers checkpoint their data
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<RtSolver>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->writeCheckpointLevel(a_handle, a_lvl);
  m_sigma->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
void CdrPlasmaStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("CdrPlasmaStepper::readCheckpointData");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::readCheckpointData" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->readCheckpointLevel(a_handle, a_lvl);
  m_sigma->readCheckpointLevel(a_handle, a_lvl);
}
#endif

int CdrPlasmaStepper::getNumberOfPlotVariables() const{
  CH_TIME("CdrPlasmaStepper::getNumberOfPlotVariables");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::getNumberOfPlotVariables" << endl;
  }
  int ncomp = 0;
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }
  
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }

  ncomp += m_fieldSolver->getNumberOfPlotVariables();
  ncomp += m_sigma->getNumberOfPlotVariables();
  ncomp += SpaceDim; // For plotting the current density

  // CdrPlasmaPhysics plotting
  ncomp += m_physics->getNumberOfPlotVariables();

  return ncomp;
}

void CdrPlasmaStepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("CdrPlasmaStepper::writePlotData");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::writePlotData" << endl;
  }

  // Poisson solver copies over its output data
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  m_fieldSolver->writePlotData(a_output, a_icomp);

  // Surface charge solver writes
  a_plotVariableNames.append(m_sigma->getPlotVariableNames());
  m_sigma->writePlotData(a_output, a_icomp);

  // CDR solvers copy their output data
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // Write the current to the output
  this->writeJ(a_output, a_icomp);
  a_plotVariableNames.push_back("x-J");
  a_plotVariableNames.push_back("y-J");
  if(SpaceDim == 3){
    a_plotVariableNames.push_back("z-J");
  }

  // Physics writes to output as well.
  a_plotVariableNames.append(m_physics->getPlotVariableNames());
  this->writePhysics(a_output, a_icomp);
}

void CdrPlasmaStepper::writeJ(EBAMRCellData& a_output, int& a_icomp) const{
  CH_TIME("CdrPlasmaStepper::writeJ");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::writeJ" << endl;
  }

  // Allocates storage and computes J
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, phase::gas, SpaceDim);
  this->computeJ(scratch);

  const Interval src_interv(0, SpaceDim-1);
  const Interval dst_interv(a_icomp, a_icomp + SpaceDim -1);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }
  a_icomp += SpaceDim;
}

void CdrPlasmaStepper::writePhysics(EBAMRCellData& a_output, int& a_icomp) const {
  CH_TIME("CdrPlasmaStepper::writePhysics");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::writePhysics" << endl;
  }

  const int numVars = m_physics->getNumberOfPlotVariables();

  if(numVars > 0){

    const int numCdrSpecies = m_physics->getNumCdrSpecies();    
    const int numRteSpecies = m_physics->getNumRtSpecies ();

    // Compute the electric field
    EBAMRCellData E;
    m_amr->allocate(E, m_realm, phase::gas, SpaceDim);
    this->computeElectricField(E, m_cdr->getPhase(), m_fieldSolver->getPotential());

    // Scratch data. 
    EBAMRCellData scratch;
    m_amr->allocate(scratch, m_realm, phase::gas, 1);

    // CDR and RTE densities
    const Vector<EBAMRCellData*> cdrDensities = m_cdr->getPhis();
    const Vector<EBAMRCellData*> rteDensities = m_rte->getPhis();

    // Compute the gradient of each species. 
    std::vector<std::shared_ptr<EBAMRCellData> > cdrGradients(numCdrSpecies);
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // This storage must be deleted later. 
      cdrGradients[idx] = std::make_shared<EBAMRCellData>();

      // Allocate cell-centered data
      m_amr->allocate(*cdrGradients[idx], m_realm, m_cdr->getPhase(), SpaceDim);

      // Copy the densities to a scratch data holder so we can compute the gradient. Must do this because
      // the gradient is a two-level AMR operator. 
      DataOps::copy(scratch, *cdrDensities[idx]);
      m_amr->interpGhostMG(scratch, m_realm, m_cdr->getPhase());

      // Now compute the gradient and coarsen/interpolate the invalid regions. 
      m_amr->computeGradient(*cdrGradients[idx], scratch, m_realm, m_cdr->getPhase());
      m_amr->averageDown    (*cdrGradients[idx],          m_realm, m_cdr->getPhase());   
      m_amr->interpGhost    (*cdrGradients[idx],          m_realm, m_cdr->getPhase());   
    }

    // This is stuff that is on a per-cell basis. We visit each cell and populate these fields and then pass them to our
    // nifty plasma physics object. 
    Vector<Real>     localCdrDensities(cdrDensities.size(), 0.0           );
    Vector<RealVect> localCdrGradients(cdrDensities.size(), RealVect::Zero);    
    Vector<Real>     localRteDensities(rteDensities.size(), 0.0           );

    // E and the gradients can be put on the centroids immediately because their lifetimes are limited to
    // this function.
    m_amr->interpToCentroids(E, m_realm, phase::gas);
    for (auto& grad : cdrGradients){
      m_amr->interpToCentroids(*grad, m_realm, phase::gas);
    }

    // Level loop.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl   = m_amr->getGrids     (m_realm            )[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, phase::gas)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];

      // Patch loop.
      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box&     cellBox = dbl  [dit()]; // <--- regular region. 
	const EBISBox& ebisBox = ebisl[dit()];

	// Irregular region. 
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];					     

	// Output data holder. 
	EBCellFAB& output = (*a_output[lvl])[dit()];
	BaseFab<Real>& regularData = (*a_output[lvl])[dit()].getSingleValuedFAB();
      
	// Regular kernel. 
	auto regularKernel = [&] (const IntVect& iv) -> void {
	  const RealVect position = m_amr->getProbLo() + (0.5*RealVect::Unit + RealVect(iv))*dx;
	  const RealVect localE   = RealVect(D_DECL((*E[lvl])[dit()].getSingleValuedFAB()(iv, 0),
						    (*E[lvl])[dit()].getSingleValuedFAB()(iv, 1),
						    (*E[lvl])[dit()].getSingleValuedFAB()(iv, 2)));

	  for (int i = 0; i < cdrDensities.size(); i++){
	    localCdrDensities[i] = (*(*cdrDensities[i])[lvl])[dit()].getSingleValuedFAB()(iv, 0);
	  }

	  for (int i = 0; i < numCdrSpecies; i++){
	    localCdrGradients[i] = RealVect(D_DECL((*(*cdrGradients[i])[lvl])[dit()].getSingleValuedFAB()(iv, 0),
						   (*(*cdrGradients[i])[lvl])[dit()].getSingleValuedFAB()(iv, 1),
						   (*(*cdrGradients[i])[lvl])[dit()].getSingleValuedFAB()(iv, 2)));
	  }

	  for (int i = 0; i < rteDensities.size(); i++){
	    localRteDensities[i] = (*(*rteDensities[i])[lvl])[dit()].getSingleValuedFAB()(iv, 0);
	  }

	  // Get plot variables from plasma physics. 
	  const Vector<Real> plotVars = m_physics->getPlotVariables(localCdrDensities,
								    localCdrGradients,
								    localRteDensities,
								    localE,
								    position,
								    dx,
								    m_dt,
								    m_time,
								    1.0);

	  for (int icomp = 0; icomp < numVars; icomp++){
	    regularData(iv, a_icomp + icomp) = plotVars[icomp];
	  }
	};

	auto irregularKernel = [&](const VolIndex& vof) -> void {
	  const RealVect position = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisBox, dx);
	  const RealVect localE   = RealVect(D_DECL((*E[lvl])[dit()](vof, 0),
						    (*E[lvl])[dit()](vof, 1),
						    (*E[lvl])[dit()](vof, 2)));

	  for (int i = 0; i < cdrDensities.size(); i++){
	    localCdrDensities[i] = 0.0;(*(*cdrDensities[i])[lvl])[dit()](vof, 0);
	  }

	  for (int i = 0; i < numCdrSpecies; i++){
	    localCdrGradients[i] = RealVect(D_DECL((*(*cdrGradients[i])[lvl])[dit()](vof, 0),
						   (*(*cdrGradients[i])[lvl])[dit()](vof, 1),
						   (*(*cdrGradients[i])[lvl])[dit()](vof, 2)));
	  }

	  for (int i = 0; i < rteDensities.size(); i++){
	    localRteDensities[i] = (*(*rteDensities[i])[lvl])[dit()](vof, 0);
	  }

	  // Get plot variables from plasma physics. 
	  const Vector<Real> plotVars = m_physics->getPlotVariables(localCdrDensities,
								    localCdrGradients,
								    localRteDensities,
								    localE,
								    position,
								    dx,
								    m_dt,
								    m_time,
								    ebisBox.volFrac(vof));

	  for (int icomp = 0; icomp < numVars; icomp++){
	    output(vof, a_icomp + icomp) = plotVars[icomp];
	  }	  
	};

	// Execute the kernels.
	BoxLoops::loop(cellBox,   regularKernel);
	BoxLoops::loop(vofit,   irregularKernel);
      }
    }

    // Need to let the outside world know that we've written to some of the variables. 
    a_icomp += numVars;
  }
}

void CdrPlasmaStepper::postCheckpointSetup(){
  CH_TIME("CdrPlasmaStepper::postCheckpointSetup");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaStepper::postCheckpointSetup" << endl;
  }

  this->solvePoisson();       // Solve Poisson equation by 
  if(this->stationaryRTE()){  // Solve RTE equations if stationary solvers
    const Real dummy_dt = 0.0;
    this->solveRadiativeTransfer(dummy_dt); // Argument does not matter, it's a stationary solver.
  }
  this->allocateInternals();  // Prepare internal storage for time stepper

  // Fill solvers with important stuff
  this->computeCdrDriftVelocities();
  this->computeCdrDiffusion();
}

void CdrPlasmaStepper::printStepReport(){
  CH_TIME("CdrPlasmaStepper::printStepReport");
  if(m_verbosity > 4){
    pout() << "CdrPlasmaStepper::printStepReport" << endl;
  }

  // Compute the maximum electric field
  Real Emax;
  this->computeMaxElectricField(Emax, phase::gas);

  //
  Real nmax;
  std::string solver_max;
  this->getCdrMax(nmax, solver_max);

  const Real cfl_dt = this->getCflDt();

  std::string str;
  if(m_timeCode == TimeCode::Advection){
    str = " (Restricted by advection)";
  }
  else if(m_timeCode == TimeCode::AdvectionDiffusion){
    str = " (Restricted by advection-diffusion)";
  }
  else if(m_timeCode == TimeCode::Error){
    str = " (Restricted by error)";
  }
  else if(m_timeCode == TimeCode::Diffusion){
    str = " (Restricted by diffusion)";
  }
  else if(m_timeCode == TimeCode::Source){
    MayDay::Error("CdrPlasmaStepper::stepReport - shouldn't happen, source term has been taken out of the design");
    str = " (Restricted by source term)";
  }
  else if(m_timeCode == TimeCode::RelaxationTime){
    str = " (Restricted by relaxation time)";
  }
  else if(m_timeCode == TimeCode::Restricted){
    str = " (Restricted by time stepper)";
  }
  else if(m_timeCode == TimeCode::Hardcap){
    str = " (Restricted by a hardcap)";
  }
  pout() << "                                   mode  = " << str << endl
	 << "                                   cfl   = " << m_dt/cfl_dt << endl
	 << "                                   Emax  = " << Emax << endl
	 << "                                   n_max = " << nmax << "(" + solver_max + ")" << endl;
}

#include <CD_NamespaceFooter.H>
