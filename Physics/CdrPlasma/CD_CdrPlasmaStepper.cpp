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
  CH_TIME("CdrPlasmaStepper::postInitialize()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::postInitialize()" << endl;
  }
}

void CdrPlasmaStepper::registerRealms(){
  CH_TIME("CdrPlasmaStepper::registerRealms()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::registerRealms()" << endl;
  }

  CH_assert(!m_amr.isNull());
  
  m_amr->registerRealm(m_realm);
}

void CdrPlasmaStepper::postRegrid(){
  CH_TIME("CdrPlasmaStepper::postRegrid()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::postRegrid()" << endl;
  }
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
						 const Vector<EBAMRCellData*>&  a_cdrDensities) {
  CH_TIME("CdrPlasmaStepper::computeSpaceChargeDensity(MFAMRCellData, Vector(EBAMRCellData))");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeSpaceChargeDensity(MFAMRCellData, Vector(EBAMRCellData))" << endl;
  }

  CH_assert(a_rho[0]->nComp() == 1);

  // Set space charge to zero everywhere. We will then compute the space charge on
  // the gas phase (solid phase space charge is always zero in the CdrPlasma module). 
  DataOps::setValue(a_rho, 0.0);

  // Get the data holder for the gas phase (the MFAMRCellData holds data on both phases)
  EBAMRCellData rhoGas = m_amr->alias(phase::gas, a_rho);


  // Go through the CDR solvers and add their space charge density. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const EBAMRCellData&             density = *(a_cdrDensities[solverIt.index()]);
    const RefCountedPtr<CdrSpecies>& spec    = solverIt.getSpecies();
    const int                        Z       = spec->getChargeNumber();

    CH_assert(density[0]->nComp() == 1);
      
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

  CH_assert(a_cellConductivity[0]->nComp() == 1);

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
  
  CH_assert(a_conductivityFace[0]->nComp() == 1);
  CH_assert(a_conductivityEB  [0]->nComp() == 1);
  CH_assert(a_conductivityCell[0]->nComp() == 1);  

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

  CH_assert(a_dt >= 0.0);

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

  CH_assert(a_conductivityFace[0]->nComp() == 1  );
  CH_assert(a_conductivityEB  [0]->nComp() == 1  );
  CH_assert(a_factor                       >= 0.0);

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

  CH_assert(a_potential[0]->nComp() == 1);
  CH_assert(a_rho      [0]->nComp() == 1);
  CH_assert(a_sigma    [0]->nComp() == 1);    

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
  CH_TIME("CdrPlasmaStepper::advanceReactionNetwork(Real, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetwork(Real, Real)" << endl;
  }

  CH_assert(a_dt >= 0.0);

  // Compute the electric field on the compute phase. 
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

  CH_assert(a_E[0]->nComp() == SpaceDim);
  CH_assert(a_dt            >= 0.0     );

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
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
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

  CH_assert(a_E[0]->nComp() == SpaceDim);
  CH_assert(a_dt            >= 0.0     );  

  // Number of CDR and RTE solvers.
  const int numCdrSpecies = m_physics->getNumCdrSpecies();
  const int numRteSpecies = m_physics->getNumRtSpecies ();

  CH_assert(a_cdrSources.  size() == numCdrSpecies);
  CH_assert(a_rteSources.  size() == numRteSpecies);  
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  
  CH_assert(a_cdrGradients.size() == numCdrSpecies);
  CH_assert(a_rteDensities.size() == numRteSpecies);      

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

  // TLDR: This is the level version of advanceReactionNetwork. It's purpose is to compute the source terms for the CDR and RTE equations. It will expose
  //       all the input parameters on a per-patch basis and then call the patch version.   

  CH_assert(a_E.nComp() == SpaceDim);
  CH_assert(a_dt        >= 0.0     );

  constexpr int comp = 0;

  // Number of species involved. 
  const int numRteSpecies  = m_physics->getNumRtSpecies();
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrSources.  size() == numCdrSpecies);
  CH_assert(a_rteSources.  size() == numRteSpecies);  
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  
  CH_assert(a_cdrGradients.size() == numCdrSpecies);
  CH_assert(a_rteDensities.size() == numRteSpecies);          

  // Stencils for putting cell-centered data on cell centroids. 
  const IrregAmrStencil<CentroidInterpolationStencil>& interpStencils = m_amr->getCentroidInterpolationStencils(m_realm, m_cdr->getPhase());

  // Grids and EBIS information for this grid level. 
  const DisjointBoxLayout& dbl   = m_amr->getGrids     (m_realm                   )[a_lvl];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real               dx    = m_amr->getDx()                                  [a_lvl];

  // Grid loop
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box cellBox = dbl[dit()];

    // Things that are passed into the kernels. The stuff with EBCellFAB* signatures is for the irregular kernel. The vectors with FArrayBox* signatures
    // are for the regular kernels. 
    Vector<EBCellFAB*> cdrSources      (numCdrSpecies, nullptr);
    Vector<EBCellFAB*> cdrDensities    (numCdrSpecies, nullptr);
    Vector<EBCellFAB*> cdrGradients    (numCdrSpecies, nullptr);
    Vector<EBCellFAB*> cdrVelocities   (numCdrSpecies, nullptr);
    Vector<EBCellFAB*> rteSources      (numRteSpecies, nullptr);
    Vector<EBCellFAB*> rteDensities    (numRteSpecies, nullptr);

    Vector<FArrayBox*> cdrSourcesFAB   (numCdrSpecies, nullptr);
    Vector<FArrayBox*> cdrDensitiesFAB (numCdrSpecies, nullptr);
    Vector<FArrayBox*> cdrGradientsFAB (numCdrSpecies, nullptr);
    Vector<FArrayBox*> cdrVelocitiesFAB(numCdrSpecies, nullptr);
    Vector<FArrayBox*> rteSourcesFAB   (numRteSpecies, nullptr);
    Vector<FArrayBox*> rteDensitiesFAB (numRteSpecies, nullptr);    

    // Fetch things for the CDR solvres. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const RefCountedPtr<CdrSolver>& solver = solverIt();
      
      const int    idx  = solverIt.index();
      
      cdrSources     [idx] = &(*a_cdrSources  [idx])[dit()];
      cdrDensities   [idx] = &(*a_cdrDensities[idx])[dit()];
      cdrGradients   [idx] = &(*a_cdrGradients[idx])[dit()];

      cdrSourcesFAB  [idx] = &(cdrSources  [idx]->getFArrayBox());
      cdrDensitiesFAB[idx] = &(cdrDensities[idx]->getFArrayBox());      
      cdrGradientsFAB[idx] = &(cdrGradients[idx]->getFArrayBox());

      // If the solver is mobile, fetch it's (cell-centered) velocity. 
      if(solver->isMobile()){
	const EBAMRCellData& velo = solver->getCellCenteredVelocity();
	
	cdrVelocities   [idx] = &((*velo[a_lvl])[dit()]);
	
	cdrVelocitiesFAB[idx] = &(cdrVelocities[idx]->getFArrayBox());
      }
    }

    // Fetch things from the RTE solvers. 
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int    idx  = solverIt.index();
      
      rteSources  [idx] = &(*a_rteSources  [idx])[dit()];
      rteDensities[idx] = &(*a_rteDensities[idx])[dit()];

      rteSourcesFAB  [idx] = &(rteSources  [idx]->getFArrayBox());
      rteDensitiesFAB[idx] = &(rteDensities[idx]->getFArrayBox());
    }
    
    // Do regular cells (actually also does irregular cells but those are redone using correct arithmetic in the call below). 
    this->advanceReactionNetworkRegularCells(cdrSourcesFAB,
					     rteSourcesFAB,
					     cdrDensitiesFAB,
					     cdrGradientsFAB,
					     rteDensitiesFAB,
					     a_E[dit()].getFArrayBox(),
					     a_time,
					     a_dt,
					     dx,
					     cellBox);

    // The regular cell loop will also do the irregular cells, so we redo the irregular cells here. 
    this->advanceReactionNetworkIrreg(cdrSources,
				      rteSources,
				      cdrDensities,
				      cdrGradients,
				      cdrVelocities,
				      rteDensities,
				      interpStencils[a_lvl][dit()],
				      a_E[dit()],
				      a_time,
				      a_dt,
				      dx,
				      cellBox,
				      a_lvl,
				      dit());


    // Covered cells are bogus.
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx  = solverIt.index();
      
      cdrSources[idx]->setCoveredCellVal(0.0, comp);
    }
    
    // Covered cells are bogus.
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx  = solverIt.index();
      
      rteSources[idx]->setCoveredCellVal(0.0, comp);
    };
  }
}

void CdrPlasmaStepper::advanceReactionNetworkRegularCells(Vector<FArrayBox*>&       a_cdrSources,
							  Vector<FArrayBox*>&       a_rteSources,
							  const Vector<FArrayBox*>& a_cdrDensities,
							  const Vector<FArrayBox*>& a_cdrGradients,
							  const Vector<FArrayBox*>& a_rteDensities,
							  const FArrayBox&          a_E,
							  const Real&               a_time,
							  const Real&               a_dt,
							  const Real&               a_dx,
							  const Box&                a_cellBox){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkRegularCells(Vector<FArrayBox*>x5, FArrayBox, Realx3, Box)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkRegularCells(Vector<FArrayBox*>x5, FArrayBox, Realx3, Box)" << endl;
  }

  CH_assert(a_E.nComp() == SpaceDim );
  CH_assert(a_dt        >= 0.0      );
  CH_assert(a_dx        >= 0.0      );
  CH_assert(a_cellBox.cellCentered());

  constexpr int  comp  = 0;

  constexpr Real zero  = 0.0;
  constexpr Real kappa = 1.0;
    
  // Number of CDR and RTE solvers. 
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();
  const int numRteSpecies  = m_physics->getNumRtSpecies();

  CH_assert(a_cdrSources.  size() == numCdrSpecies);
  CH_assert(a_rteSources.  size() == numRteSpecies);  
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  
  CH_assert(a_cdrGradients.size() == numCdrSpecies);
  CH_assert(a_rteDensities.size() == numRteSpecies);        

  // Lower-left corner -- physical coordinates. 
  const RealVect probLo  = m_amr->getProbLo();

  // Create some transient storage on which we store the source terms for this grid patch. 
  FArrayBox cdrSrc(a_cellBox, numCdrSpecies);
  FArrayBox rteSrc(a_cellBox, numRteSpecies);

  cdrSrc.setVal(0.0);
  rteSrc.setVal(0.0);

  // These things that are passed into CdrPlasmaPhysics
  Vector<Real    > cdrSources  (numCdrSpecies, 0.0           );
  Vector<Real    > rteSources  (numRteSpecies, 0.0           );  
  Vector<Real    > cdrDensities(numCdrSpecies, 0.0           );
  Vector<RealVect> cdrGradients(numCdrSpecies, RealVect::Zero);
  Vector<Real    > rteDensities(numRteSpecies, 0.0           );  
  
  // Regular kernel. We reconstructor the various cell-centered quantities and put them in the data structure required
  // by CdrPlasmaPhysics. 
  auto regularKernel = [&](const IntVect& iv) -> void {
    // Create the position and electric field. 
    const RealVect pos = probLo + (0.5*RealVect::Unit + RealVect(iv)) * a_dx;
    const RealVect E   = RealVect(D_DECL(a_E(iv, 0), a_E(iv, 1), a_E(iv, 2)));

    // Get the cell-centered CDR density and gradient. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int  idx = solverIt.index();
      const Real phi = (*a_cdrDensities[idx])(iv, comp);
	
      cdrDensities[idx] = std::max(zero, phi);
      cdrGradients[idx] = RealVect(D_DECL((*a_cdrGradients[idx])(iv, 0), (*a_cdrGradients[idx])(iv, 1), (*a_cdrGradients[idx])(iv, 2)));
    }

    // Get the cell-centered radiative transfer densities. 
    for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int  idx = solverIt.index();
      const Real phi = (*a_rteDensities[idx])(iv, comp);
	
      rteDensities[idx] = std::max(zero, phi);
    }

    // Physics now solves for the source terms. 
    m_physics->advanceReactionNetwork(cdrSources,
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

    // Put CDR source terms into temporary data holders. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
	
      cdrSrc(iv, idx) = cdrSources[idx];
    }
    
    // Put RTE source terms into temporary dataholders. 
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
	
      rteSrc(iv, idx) = rteSources[idx];
    }
  };

  // Execute the kernel. This puts the source terms in cdrSrc and rteSrc, but our target data holders are the input data holders
  // with single components. So, copy the result back to these.
  BoxLoops::loop(a_cellBox, regularKernel);  

  // Do it for the CDR solvers. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    (*a_cdrSources[idx]).setVal(0.0);
    (*a_cdrSources[idx]).plus(cdrSrc, idx, comp);
  }

  // Do it for the RTE solvers. 
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    (*a_rteSources[idx]).setVal(0.0);
    (*a_rteSources[idx]).plus(rteSrc, idx, comp);
  }
}

void CdrPlasmaStepper::advanceReactionNetworkIrreg(Vector<EBCellFAB*>&          a_cdrSources,
						   Vector<EBCellFAB*>&          a_rteSources,
						   const Vector<EBCellFAB*>&    a_cdrDensities,
						   const Vector<EBCellFAB*>&    a_cdrGradients,
						   const Vector<EBCellFAB*>&    a_cdrVelocities,
						   const Vector<EBCellFAB*>&    a_rteDensities,
						   const BaseIVFAB<VoFStencil>& a_interpStencils,
						   const EBCellFAB&             a_E,
						   const Real&                  a_time,
						   const Real&                  a_dt,
						   const Real&                  a_dx,
						   const Box&                   a_cellBox,
						   const int                    a_lvl,
						   const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkIrreg(Vector<EBCellFAB*>x6, ...)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkRegularCells(Vector<EBCellFAB*>x6, ...)" << endl;
  }

  CH_assert(a_E.nComp() == SpaceDim );
  CH_assert(a_dt        >= 0.0      );
  CH_assert(a_dx        >= 0.0      );
  CH_assert(a_cellBox.cellCentered());
  
  switch(m_whichSourceTermComputation){
  case SourceTermComputation::Interpolated:
    {
      this->advanceReactionNetworkIrregInterp(a_cdrSources,
					      a_rteSources,
					      a_cdrDensities,
					      a_cdrGradients,
					      a_cdrVelocities,
					      a_rteDensities,
					      a_interpStencils,
					      a_E,
					      a_time,
					      a_dt,
					      a_dx,
					      a_cellBox,
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
					      a_cdrVelocities,					      
					      a_rteDensities,
					      a_interpStencils,
					      a_E,
					      a_time,
					      a_dt,
					      a_dx,
					      a_cellBox,
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
					     a_interpStencils,

					     a_E,
					     a_time,
					     a_dt,
					     a_dx,
					     a_cellBox,
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
					      a_cdrVelocities,					      
					      a_rteDensities,
					      a_interpStencils,
					      a_E,
					      a_time,
					      a_dt,
					      a_dx,
					      a_cellBox,
					      a_lvl,
					      a_dit);      
      break;	
    }
  default:
    {
      MayDay::Error("CdrPlasmaStepper::advanceReactionNetworkIrreg - logic bust");
      
      break;
    }
  }
}

void CdrPlasmaStepper::advanceReactionNetworkIrregInterp(Vector<EBCellFAB*>&          a_cdrSources,
							 Vector<EBCellFAB*>&          a_rteSources,
							 const Vector<EBCellFAB*>&    a_cdrDensities,
							 const Vector<EBCellFAB*>&    a_cdrGradients,
							 const Vector<EBCellFAB*>&    a_cdrVelocities,
							 const Vector<EBCellFAB*>&    a_rteDensities,
							 const BaseIVFAB<VoFStencil>& a_interpStencils,
							 const EBCellFAB&             a_E,
							 const Real&                  a_time,
							 const Real&                  a_dt,
							 const Real&                  a_dx,
							 const Box&                   a_cellBox,
							 const int                    a_lvl,
							 const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkIrregInterp(Vector<EBCellFAB*>x6, ...)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkIrregInterp(Vector<EBCellFAB*>x6, ...)" << endl;
  }

  CH_assert(a_E.nComp() == SpaceDim );
  CH_assert(a_dt        >= 0.0      );
  CH_assert(a_dx        >= 0.0      );
  CH_assert(a_cellBox.cellCentered());  

  constexpr int  comp = 0  ;
  constexpr Real zero = 0.0;

  // Number of CDR and RTE solvers.
  const int numRteSpecies  = m_physics->getNumRtSpecies ();
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrSources.  size() == numCdrSpecies);
  CH_assert(a_rteSources.  size() == numRteSpecies);  
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  
  CH_assert(a_cdrGradients.size() == numCdrSpecies);
  CH_assert(a_rteDensities.size() == numRteSpecies);        

  // EB box. 
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];

  // Lower-left corner in physical coordinates
  const RealVect probLo = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  Vector<Real>     cdrSources  (numCdrSpecies, 0.            );
  Vector<Real>     cdrDensities(numCdrSpecies, 0.            );
  Vector<Real>     rteSources  (numRteSpecies, 0.            );  
  Vector<Real>     rteDensities(numRteSpecies, 0.            );
  Vector<RealVect> cdrGradients(numCdrSpecies, RealVect::Zero);  

  // Irregular grid kernel. 
  auto irregularKernel = [&](const VolIndex& vof) -> void{
    const Real        kappa   = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interpStencils(vof, comp);
    const RealVect    pos     = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, a_dx);
    const RealVect    normal  = ebisbox.normal(vof);

    // I assume that the electric field is on the cell center (which is what FieldSolverMultigrid currently does), but we do want the field on the centroid. 
    RealVect E = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){

      // Do all components. 
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += stencil.weight(i) * a_E(stencil.vof(i), dir);
      }
    }

    // Compute RTE densities on the centroids. Again, I sort of assume that the solver being used
    // is a cell-centered solver so we interpolate the isotropic term to the centroid. 
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      rteDensities[idx] = 0.0;
      
      for (int i = 0; i < stencil.size(); i++){
	rteDensities[idx] += stencil.weight(i) * (*a_rteDensities[idx])(stencil.vof(i), comp);
      }
      
      rteDensities[idx] = std::max(zero, rteDensities[idx]);
    }

    // Compute plasma species densities on the cell centroid. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      cdrDensities[idx] = 0.0;
      
      for (int i = 0; i < stencil.size(); i++){
	cdrDensities[idx] += stencil.weight(i) * (*a_cdrDensities[idx])(stencil.vof(i), comp);
      }
      
      cdrDensities[idx] = std::max(cdrDensities[idx], zero);
    }

    // Compute plasma species gradients on the cell centroid. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();    

      cdrGradients[idx] = RealVect::Zero;
      
      for (int i = 0; i < stencil.size(); i++){

	// Do all components. 
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
				      rteDensities,
				      E,
				      pos,
				      a_dx,
				      a_dt,
				      a_time,
				      kappa);

    // Iterate through the CDR solvers and set the source terms. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      (*a_cdrSources[idx])(vof, comp) = cdrSources[idx];
    }

    // Iterate through the RTE solvers and set the source terms.     
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      (*a_rteSources[idx])(vof, comp) = rteSources[idx];
    }
  };

  // Region for kernel. 
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];

  // Execute the kernel. 
  BoxLoops::loop(vofit, irregularKernel);
}

void CdrPlasmaStepper::advanceReactionNetworkIrregKappa(Vector<EBCellFAB*>&          a_cdrSources,
							Vector<EBCellFAB*>&          a_rteSources,
							const Vector<EBCellFAB*>&    a_cdrDensities,
							const Vector<EBCellFAB*>&    a_cdrGradients,
							const Vector<EBCellFAB*>&    a_rteDensities,
							const BaseIVFAB<VoFStencil>& a_interpStencils,
							const EBCellFAB&             a_E,
							const Real&                  a_time,
							const Real&                  a_dt,
							const Real&                  a_dx,
							const Box&                   a_cellBox,
							const int                    a_lvl,
							const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkIrregKappa(Vector<EBCellFAB*>x6, ...)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkIrregKappa(Vector<EBCellFAB*>x6, ...)" << endl;
  }

  CH_assert(a_E.nComp() == SpaceDim );
  CH_assert(a_dt        >= 0.0      );
  CH_assert(a_dx        >= 0.0      );
  CH_assert(a_cellBox.cellCentered());  

  constexpr int  comp = 0  ;
  constexpr Real zero = 0.0;

  // Number of CDR and RTE solvers.
  const int numRteSpecies  = m_physics->getNumRtSpecies ();
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrSources.  size() == numCdrSpecies);
  CH_assert(a_rteSources.  size() == numRteSpecies);  
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  
  CH_assert(a_cdrGradients.size() == numCdrSpecies);
  CH_assert(a_rteDensities.size() == numRteSpecies);    

  // EB box. 
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];    

  // Lower-left corner in physical coordinates.
  const RealVect probLo  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  Vector<Real>     cdrSources  (numCdrSpecies, 0.            );
  Vector<Real>     cdrDensities(numCdrSpecies, 0.            );
  Vector<Real>     rteSources  (numRteSpecies, 0.            );  
  Vector<Real>     rteDensities(numRteSpecies, 0.            );
  Vector<RealVect> cdrGradients(numCdrSpecies, RealVect::Zero);    

  // Definition of the irregular kernel. 
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const Real        kappa   = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interpStencils(vof, comp);
    const RealVect    pos     = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, a_dx);

    // Input E is on cell center but we need centroid. 
    RealVect E = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += a_E(stencil.vof(i), dir) * stencil.weight(i);
      }
    }

    // On input the CDR densities are on the centroid so just make sure they are positive. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      cdrDensities[idx] = std::max(zero, (*a_cdrDensities[idx])(vof, comp));      
    }

    // On input we assume that the CDR gradients are on the cell center so we interpolate to centroid. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      cdrGradients[idx] = RealVect::Zero;
      for (int i = 0; i < stencil.size(); i++){
	for (int dir = 0; dir < SpaceDim; dir++){
	  cdrGradients[idx][dir] += (*a_cdrGradients[idx])(stencil.vof(i), dir) * stencil.weight(i);
	}
      }
    }

    // On input the RTE densities are on the centroid, so just make sure it's positive. 
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      rteDensities[idx] = std::max(zero, (*a_rteDensities[idx])(vof, comp));
    }

    // Compute source terms
    m_physics->advanceReactionNetwork(cdrSources,
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

    // Put the CDR source term where it belongs. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      (*a_cdrSources[idx])(vof, comp) = cdrSources[idx];
    }

    // Put the RTE source term where it belongs.     
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      (*a_rteSources[idx])(vof, comp) = rteSources[idx];
    }
  };

  // Kernel region. 
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];

  // Run the kernel
  BoxLoops::loop(vofit, irregularKernel);
}

void CdrPlasmaStepper::advanceReactionNetworkIrregUpwind(Vector<EBCellFAB*>&          a_cdrSources,
							 Vector<EBCellFAB*>&          a_rteSources,
							 const Vector<EBCellFAB*>&    a_cdrDensities,
							 const Vector<EBCellFAB*>&    a_cdrGradients,
							 const Vector<EBCellFAB*>&    a_cdrVelocities,
							 const Vector<EBCellFAB*>&    a_rteDensities,
							 const BaseIVFAB<VoFStencil>& a_interpStencils,
							 const EBCellFAB&             a_E,
							 const Real&                  a_time,
							 const Real&                  a_dt,
							 const Real&                  a_dx,
							 const Box&                   a_cellBox,
							 const int                    a_lvl,
							 const DataIndex&             a_dit){
  CH_TIME("CdrPlasmaStepper::advanceReactionNetworkIrregUpwind(Vector<EBCellFAB*>x6, ...)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::advanceReactionNetworkIrregUpwind(Vector<EBCellFAB*>x6, ...)" << endl;
  }

  constexpr int  comp = 0  ;
  constexpr Real zero = 0.0;

  // Number of CDR and RTE solvers.
  const int numRteSpecies  = m_physics->getNumRtSpecies ();
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrSources.  size() == numCdrSpecies);
  CH_assert(a_rteSources.  size() == numRteSpecies);  
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  
  CH_assert(a_cdrGradients.size() == numCdrSpecies);
  CH_assert(a_rteDensities.size() == numRteSpecies);        

  // EB box. 
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];    

  // Lower-left corner in physical coordinates.
  const RealVect probLo  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  Vector<Real>     cdrSources  (numCdrSpecies, 0.            );
  Vector<Real>     cdrDensities(numCdrSpecies, 0.            );
  Vector<Real>     rteSources  (numRteSpecies, 0.            );  
  Vector<Real>     rteDensities(numRteSpecies, 0.            );
  Vector<RealVect> cdrGradients(numCdrSpecies, RealVect::Zero);    

  // This is the kernel that we run. 
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const Real        kappa   = ebisbox.volFrac(vof);
    const VoFStencil& stencil = a_interpStencils(vof, comp);
    const RealVect    pos     = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, a_dx);

    // The input electric field is on the cell center but we need it on the centroid. 
    RealVect E = RealVect::Zero;
    for (int i = 0; i < stencil.size(); i++){
      for (int dir = 0; dir < SpaceDim; dir++){
	E[dir] += stencil.weight(i) * a_E(stencil.vof(i), dir);
      }
    }

    // Project the electric field on the EB normal (recall, it points inwards into the plasma)
    const Real EdotN   = E.dotProduct(ebisbox.normal(vof));

    // Compute the CDR densities that go into the reaction kernel. Note that since we upwind, we check if the
    // flow is into or away from the boundary. If the flow is away from the boundary we don't really have an upwind
    // side so we set the density to zero in that case. This is not really captured by CdrSolver::weightedUpwind because
    // the fallback option is to use the cell-centered value (that is the correct design, when we don't have reactive plasmas). We fix
    // that here.
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      // Get the charge number. 
      const RefCountedPtr<CdrSpecies> species = solverIt.getSpecies();      
      const int Z = species->getChargeNumber();

      // Check if the EB is an "inflow" face. If it is, we will turn off this species in the reaction network. 
      bool inflow = false;
      if(solverIt()->isMobile() && Real(Z) * EdotN > 0.0){
	inflow = true;
      }

      // If it's not an inflow face, interpolate the data to the cell centroid. Otherwise, set the
      // species density to zero. 
      if(!inflow){
	cdrDensities[idx] = zero;
	
	for (int i = 0; i < stencil.size(); i++){
	  cdrDensities[idx] += stencil.weight(i) * (*a_cdrDensities[idx])(stencil.vof(i), comp);
	}

	// Enforce positivity in reaction kernel (interpolation stencils might have negative weights). 
	cdrDensities[idx] = std::max(zero, cdrDensities[idx]);
      }
      else{
	cdrDensities[idx] = zero;
      }
    }

    // Interpolate cell-centered gradients to centroid. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      cdrGradients[idx] = RealVect::Zero;      
      for (int i = 0; i < stencil.size(); i++){
	for (int dir = 0; dir < SpaceDim; dir++){
	  cdrGradients[idx][dir] += stencil.weight(i) * (*a_cdrGradients[idx])(stencil.vof(i), dir);
	}
      }
    }

    // Interpolate the RTE densities to the centroid. 
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // Apply the stencil but always leave a positive number of photons (interpolators might not always have non-negative weights). 
      rteDensities[idx] = 0.0;
      
      for (int i = 0; i < stencil.size(); i++){
	rteDensities[idx] += stencil.weight(i) * (*a_rteDensities[idx])(stencil.vof(i), comp);
      }
      
      rteDensities[idx] = std::max(zero, rteDensities[idx]);
    }    

    // Call the plasma-kinetics framework and have it fill the source terms over a time step a_dt. 
    m_physics->advanceReactionNetwork(cdrSources,
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

    // Iterate through the CDR solvers and set the source terms. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      (*a_cdrSources[idx])(vof, comp) = cdrSources[idx];
    }

    // Iterate through the RTE solvers and set the source terms.     
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      (*a_rteSources[idx])(vof, comp) = rteSources[idx];
    }
  };

  // Kernel region. 
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];

  // Run the kernel
  BoxLoops::loop(vofit, irregularKernel);  
}

void CdrPlasmaStepper::computeCdrDiffusion(){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusion()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusion()" << endl;
  }

  // Allcoate storage for the electric field. 
  EBAMRCellData electricFieldCell;
  EBAMRIVData   electricFieldEB  ;
  
  m_amr->allocate(electricFieldCell, m_realm, m_cdr->getPhase(), SpaceDim);
  m_amr->allocate(electricFieldEB,   m_realm, m_cdr->getPhase(), SpaceDim);

  // Compute field on cell center and on EB centroid. 
  this->computeElectricField(electricFieldCell, m_cdr->getPhase(), m_fieldSolver->getPotential());
  this->computeElectricField(electricFieldEB,   m_cdr->getPhase(), electricFieldCell            );

  // Call the other version. 
  this->computeCdrDiffusion(electricFieldCell, electricFieldEB);
}

void CdrPlasmaStepper::computeCdrDiffusion(const EBAMRCellData& a_electricFieldCell, const EBAMRIVData& a_electricFieldEB) {
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusion(EBAMRCellData, EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusion(EBAMRCellData, EBAMRIVData)" << endl;
  }

  // TLDR: We want to compute the CDR solver diffusion coefficients on the face centers (very important) and on the EB (less important). When we do the faces
  //       we actually compute the diffusion coefficient on the cell center, which we later average to faces. This is what computeCdrDiffusionFace does. On the EB,
  //       we instead extrapolate the things that we need directly to it, and call the same physics coupling function as we did for the cells.

  CH_assert(a_electricFieldCell[0]->nComp() == SpaceDim);
  CH_assert(a_electricFieldEB  [0]->nComp() == SpaceDim);    

  constexpr int numComp       = 1;

  // First, compute the diffusion coefficients on face centers.

  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  // We will fill cdrDcoFace (face-centered diffusion coefficients) and cdrDcoEB. cdrDensities are the cell-centered densities from
  // the CDR solver. 
  Vector<EBAMRFluxData*>       cdrDcoFace   = m_cdr->getFaceCenteredDiffusionCoefficient();
  Vector<EBAMRIVData*>         cdrDcoEB     = m_cdr->getEbCenteredDiffusionCoefficient();
  const Vector<EBAMRCellData*> cdrDensities = m_cdr->getPhis();

  // 1. Compute the diffusion coefficients on faces.
  this->computeCdrDiffusionFace(cdrDcoFace, cdrDensities, a_electricFieldCell, m_time);  

  // 2a. Allocate some storage that allows us to extrapolate the CDR densities to the EB so we can call the same physics on the EB instead of the
  //     grid faces. 
  Vector<EBAMRIVData*> cdrDensitiesExtrap(numCdrSpecies, nullptr);
  
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();

    // Allocate storage. Note that because we use the new operator we must the the storage as well. 
    cdrDensitiesExtrap[idx] = new EBAMRIVData();  // This must be deleted
    m_amr->allocate(*cdrDensitiesExtrap[idx], m_realm, m_cdr->getPhase(), numComp);

    // Extrapolate the cell-centered densities to the EB. 
    const IrregAmrStencil<EbCentroidInterpolationStencil>& stencil = m_amr->getEbCentroidInterpolationStencils(m_realm, m_cdr->getPhase());
    
    stencil.apply(*cdrDensitiesExtrap[idx], *cdrDensities[idx]);
  }

  // 2b. Compute the diffusion coefficeints on the EB.
  this->computeCdrDiffusionEb(cdrDcoEB,  cdrDensitiesExtrap, a_electricFieldEB,   m_time);

  // 2c. Release the extra storage allocate in 2a. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    m_amr->deallocate(*cdrDensitiesExtrap[idx]);
    
    delete cdrDensitiesExtrap[idx];
  }
}

void CdrPlasmaStepper::computeCdrDiffusionCell(Vector<EBAMRCellData>&        a_cdrDcoCell,
					       const Vector<EBAMRCellData*>& a_cdrDensities,
					       const EBAMRCellData&          a_electricFieldCell,
					       const Real&                   a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCell(Vector<EBAMRCellData*>, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCell(Vector<EBAMRCellData*>, Vector<EBAMRCellData*>, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_electricFieldCell[0]->nComp() == SpaceDim);  

  // Number of CDR solvers. 
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrDcoCell.  size() == numCdrSpecies);
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  

  // List of CDR diffusion coefficients and states. Will be populated on each grid level. 
  Vector<LevelData<EBCellFAB>* > cdrDcoCell  (numCdrSpecies, nullptr);
  Vector<LevelData<EBCellFAB>* > cdrDensities(numCdrSpecies, nullptr);  

  // Level loop
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    // Populate the Vectors.     
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      cdrDcoCell  [idx] =   a_cdrDcoCell  [idx] [lvl];
      cdrDensities[idx] = (*a_cdrDensities[idx])[lvl];
    }

    // Call the level version
    this->computeCdrDiffusionCell(cdrDcoCell, cdrDensities, *a_electricFieldCell[lvl], lvl, a_time);
  }
}

void CdrPlasmaStepper::computeCdrDiffusionCell(Vector<LevelData<EBCellFAB>* >&       a_cdrDcoCell,
					       const Vector<LevelData<EBCellFAB>* >& a_cdrDensities,
					       const LevelData<EBCellFAB>&           a_electricFieldCell,
					       const int                             a_lvl,
					       const Real&                           a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCell(Vector<LD<EBCellFAB>* >x2, LD<EBCellFAB, int, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCell(Vector<LD<EBCellFAB>* >x2, LD<EBCellFAB, int, Real)" << endl;
  }

  CH_assert(a_electricFieldCell.nComp() == SpaceDim);

  constexpr int comp = 0;

  // Number of CDR solvers. 
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrDcoCell.  size() == numCdrSpecies);
  CH_assert(a_cdrDensities.size() == numCdrSpecies);    

  // Grids on this level. 
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  
  for (DataIterator dit(dbl); dit.ok(); ++dit){

    // Computational region for the kernels. 
    const Box cellBox = dbl[dit()];

    // Handle to the electric field. 
    const EBCellFAB& electricField    = a_electricFieldCell[dit()];
    const FArrayBox& electricFieldFAB = electricField.getFArrayBox();

    // Get handle to the cell-centered diffusion coefficients and densities for the current patch. 
    Vector<EBCellFAB*> cdrDcoCell  (numCdrSpecies, nullptr);
    Vector<EBCellFAB*> cdrDensities(numCdrSpecies, nullptr);

    Vector<FArrayBox*> cdrDcoCellFAB  (numCdrSpecies, nullptr);
    Vector<FArrayBox*> cdrDensitiesFAB(numCdrSpecies, nullptr);    

    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      cdrDcoCell     [idx] = &(*a_cdrDcoCell  [idx])[dit()];
      cdrDensities   [idx] = &(*a_cdrDensities[idx])[dit()];

      cdrDcoCellFAB  [idx] = &(cdrDcoCell  [idx]->getFArrayBox());
      cdrDensitiesFAB[idx] = &(cdrDensities[idx]->getFArrayBox());
    }

    // Do the regular cells. This will also do irregular cells but we re-do those below. 
    this->computeCdrDiffusionCellRegular(cdrDcoCellFAB,
					 cdrDensitiesFAB,
					 electricFieldFAB,
					 cellBox,
					 m_amr->getDx()[a_lvl],
					 a_time);

    // Re-do the irregular cells
    this->computeCdrDiffusionCellIrregular(cdrDcoCell,
					   cdrDensities,
					   electricField,
					   m_amr->getDx()[a_lvl],
					   a_time,
					   a_lvl,
					   dit());

    // Covered data is bogus. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int  idx = solverIt.index();

      cdrDcoCell[idx]->setCoveredCellVal(0.0, comp);
    }
  }
}

void CdrPlasmaStepper::computeCdrDiffusionCellRegular(Vector<FArrayBox*>&       a_cdrDcoCell,
						      const Vector<FArrayBox*>& a_cdrDensities,
						      const FArrayBox&          a_electricFieldCell,
						      const Box                 a_cellBox,
						      const Real                a_dx,
						      const Real                a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCellRegular(Vector<FArrayBox*>x2, FArrayBox, Box, Real, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCellRegular(Vector<FArrayBox*>x2, FArrayBox, Box, Real, Real)" << endl;
  }

  CH_assert(a_electricFieldCell.nComp() == SpaceDim);
  CH_assert(a_cellBox.cellCentered());

  constexpr Real zero = 0.0;
  constexpr int  comp = 0;

  // Number of CDR solvers. 
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrDcoCell.  size() == numCdrSpecies);
  CH_assert(a_cdrDensities.size() == numCdrSpecies);      

  // Lower-left corner of physical domain. 
  const RealVect probLo  = m_amr->getProbLo();

  // This is populated by m_physics. 
  Vector<Real> cdrDensities(numCdrSpecies, 0.0);

  // Computed coefficients terms go onto here -- then we extract the data later. 
  FArrayBox diffusionCoefficients(a_cellBox, numCdrSpecies);

  // Regular kernel. We just fetch everything on the cell center and then call the physics framework. 
  auto regularKernel = [&](const IntVect iv) -> void {

    // Physical position and electric field. 
    const RealVect pos = probLo + (0.5*RealVect::Unit + RealVect(iv)) * a_dx;
    const RealVect E   = RealVect(D_DECL(a_electricFieldCell(iv, 0), a_electricFieldCell(iv, 1), a_electricFieldCell(iv, 2)));

    // Get the CDR densities in the current cell. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int  idx = solverIt.index();
      
      const Real phi = (*a_cdrDensities[idx])(iv, comp);
      
      cdrDensities[idx] = std::max(zero, phi);
    }

    // Compute the diffusion coefficients. 
    const Vector<Real> Dcos = m_physics->computeCdrDiffusionCoefficients(a_time, pos, E, cdrDensities);

    // Set the diffusion coefficients in the temporary data holder. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      diffusionCoefficients(iv, idx) = Dcos[idx];
    }
  };

  // Launch the kernel over the input box. 
  BoxLoops::loop(a_cellBox, regularKernel);

  // Linearize back -- we had all the diffusion coefficients stored in the temporary data holder -- now put them in the output data holders. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();

    // The check here is important because the input FArrayBox points can be null if the species is not diffusive.     
    if(solverIt()->isDiffusive()) {  
      (*a_cdrDcoCell[idx]).setVal(0.0);
      (*a_cdrDcoCell[idx]).plus(diffusionCoefficients, idx, comp);
    }
  }
}

void CdrPlasmaStepper::computeCdrDiffusionCellIrregular(Vector<EBCellFAB*>&          a_cdrDcoCell,
							const Vector<EBCellFAB*>&    a_cdrDensities,
							const EBCellFAB&             a_electricFieldCell,
							const Real                   a_dx,
							const Real&                  a_time,
							const int                    a_lvl,
							const DataIndex&             a_dit) {
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionCellIrregular(Vector<EBCellFAB*>x2, EBCellFAB, Real, Real, int DataIndex)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionCellIrregular(Vector<EBCellFAB*>x2, EBCellFAB, Real, Real, int DataIndex)" << endl;
  }

  CH_assert(a_electricFieldCell.nComp() == SpaceDim);

  constexpr int comp = 0;

  // Number of CDR solvers
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrDcoCell.  size() == numCdrSpecies);
  CH_assert(a_cdrDensities.size() == numCdrSpecies);  

  // EBISBox and graph
  const EBISBox& ebisbox = a_electricFieldCell.getEBISBox();

  // Lower left corner (physical coordinates)
  const RealVect probLo  = m_amr->getProbLo();

  // Things that are passed into CdrPlasmaPhysics
  Vector<Real> cdrDensities(numCdrSpecies, 0.0);

  // Irreegular kernel. 
  auto irregularKernel = [&](const VolIndex& vof) -> void {

    // Physical position and electric field. 
    const RealVect pos = probLo + Location::position(Location::Cell::Center, vof, ebisbox, a_dx);						     
    const RealVect E   = RealVect(D_DECL(a_electricFieldCell(vof, 0), a_electricFieldCell(vof, 1), a_electricFieldCell(vof, 2)));

    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx  = solverIt.index();

      const Real phi = (*a_cdrDensities[idx])(vof, comp);

      cdrDensities[idx] = std::max(0.0, phi);
    }

    // Let our nifty plasma physics framework take care of the diffusion coefficients. 
    const Vector<Real> Dcos = m_physics->computeCdrDiffusionCoefficients(a_time,
									 pos,
									 E,
									 cdrDensities);

    // Put the result in the correct diffusion coefficient data holder. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // The check here is important because the input EBCellFAB* will be nullptr if the species is not diffusive. 
      if(solverIt()->isDiffusive()){
	(*a_cdrDcoCell[idx])(vof, comp) = Dcos[idx];
      }
    }
  };

  // Kernel region.
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];  

  // Launch the kernel.
  BoxLoops::loop(vofit, irregularKernel);
}

void CdrPlasmaStepper::computeCdrDiffusionFace(Vector<EBAMRFluxData*>&       a_cdrDcoFace,
					       const Vector<EBAMRCellData*>& a_cdrDensities,
					       const EBAMRCellData&          a_electricFieldCell,
					       const Real&                   a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionFace(Vector<EBAMRFluxData*>, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionFace(Vector<EBAMRFluxData*>, Vector<EBAMRCellData*>, EBAMRCellData, Real)" << endl;
  }

  // TLDR: This is a routine for computing the CDR diffusion coefficients on face centers (not centroids!). It does so by first
  //       computing the diffusion coefficient on cell centers and then averaging those to face centers. To do this, we need
  //       to allocate some transient memory (which is automatically released).
  //
  //       The cell-to-face averaging will also fill one ghost face. This is important because the CDR solvers will interpolate
  //       face-centered fluxes to face centroids and will thus reach into one ghost face.

  CH_assert(a_electricFieldCell[0]->nComp() == SpaceDim);

  const int nComp = 1;

  // Number of CDR solvers
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrDcoFace.  size() == numCdrSpecies);
  CH_assert(a_cdrDensities.size() == numCdrSpecies);    

  // Allocate data for cell-centered diffusion coefficients. 
  Vector<EBAMRCellData> cdrDcoCell(numCdrSpecies);
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    m_amr->allocate(cdrDcoCell[idx], m_realm, m_cdr->getPhase(), nComp);
  }

  // Compute the diffusion coefficients on cell centers. 
  this->computeCdrDiffusionCell(cdrDcoCell, a_cdrDensities, a_electricFieldCell, a_time);

  // Now compute face-centered things by taking the average of cell-centered things. Note that this will include
  // one ghost face because the CDR solvers will interpolate the face-centered diffusion coefficients to face centroids. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    
    // The check here is important because some EBAMRFluxData input pointers will be nullptr when the species is not diffusive. 
    if(solver->isDiffusive()){ 
      const int idx = solverIt.index();

      // Set to something stupendously large to spot errors. 
      DataOps::setValue(*a_cdrDcoFace[idx], std::numeric_limits<Real>::max());

      // Coarsen the cell-centered diffusion coefficient before averaging to faces. 
      m_amr->averageDown  (cdrDcoCell[idx], m_realm, m_cdr->getPhase());
      m_amr->interpGhostMG(cdrDcoCell[idx], m_realm, m_cdr->getPhase()); 

      // Average to cell faces. Note that this call also includes one ghost face. 
      DataOps::averageCellScalarToFaceScalar(*a_cdrDcoFace[idx], cdrDcoCell[idx], m_amr->getDomains());
    }
  }
}


void CdrPlasmaStepper::computeCdrDiffusionEb(Vector<EBAMRIVData*>&       a_cdrDcoEB,
					     const Vector<EBAMRIVData*>& a_cdrDensitiesEB,
					     const EBAMRIVData&          a_electricFieldEB,
					     const Real&                 a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionEb(Vector<EBAMRIVData*>x2, EBAMRIVData, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionEb(Vector<EBAMRIVData*>x2, EBAMRIVData, Real)" << endl;
  }

  // TLDR: This is responsible for computing the CDR diffusion coefficients on EB faces. This is necessary because the finite-volume
  //       method needs it. This routine will go through each level and call an equivalent routine.

  CH_assert(a_electricFieldEB[0]->nComp() == SpaceDim);  

  // Numer of CDR species. 
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrDcoEB  .    size() == numCdrSpecies);
  CH_assert(a_cdrDensitiesEB.size() == numCdrSpecies);  

  // Level loop.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    // Populate the level-view of the diffusion coefficients and CDR densities. 
    Vector<LevelData<BaseIVFAB<Real> >* > cdrDcoEB      (numCdrSpecies, nullptr);
    Vector<LevelData<BaseIVFAB<Real> >* > cdrDensitiesEB(numCdrSpecies, nullptr);

    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      cdrDcoEB      [idx] = (*a_cdrDcoEB      [idx])[lvl];
      cdrDensitiesEB[idx] = (*a_cdrDensitiesEB[idx])[lvl];
    }

    // Call the level version.
    this->computeCdrDiffusionEb(cdrDcoEB, cdrDensitiesEB, *a_electricFieldEB[lvl], a_time, lvl);
  }
}


void CdrPlasmaStepper::computeCdrDiffusionEb(Vector<LevelData<BaseIVFAB<Real> >* >&       a_cdrDcoEB,
					     const Vector<LevelData<BaseIVFAB<Real> >* >& a_cdrDensitiesEB,
					     const LevelData<BaseIVFAB<Real> >&           a_electricFieldEB,
					     const Real&                                  a_time,
					     const int                                    a_lvl){
  CH_TIME("CdrPlasmaStepper::computeCdrDiffusionEb(Vector<LD<BaseIVFAB<Real> >* >x2, LD<BaseIVFAB<Real> >, Real, int)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDiffusionEb(Vector<LD<BaseIVFAB<Real> >* >x2, LD<BaseIVFAB<Real> >, Real, int)" << endl;
  }

  // TLDR: This is the level version which computes the CDR diffusion coefficients on EB centroids. It will go through all patches and
  //       call an irregular kernel. 

  CH_assert(a_electricFieldEB.nComp() == SpaceDim);    

  constexpr int  comp = 0  ;
  constexpr Real zero = 0.0;

  // Number of CDR solvers
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_cdrDcoEB.      size() == numCdrSpecies);
  CH_assert(a_cdrDensitiesEB.size() == numCdrSpecies);

  // Physical coordinates of lower left corner 
  const RealVect probLo = m_amr->getProbLo();

  // This is passed into the plasma kinetics framework. 
  Vector<Real> cdrDensitiesEB(numCdrSpecies, 0.0);  

  // Grids on this level. 
  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real               dx    = m_amr->getDx()[a_lvl];

  // Grid patch loop.
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const EBISBox& ebisbox   = ebisl[dit()];

    // Electric field in this grid patch. 
    const BaseIVFAB<Real>& electricFieldPatchEB = a_electricFieldEB[dit()];

    // Irregular kernel.
    auto irregularKernel = [&] (const VolIndex& vof) -> void {

      // Physical coordinates and electric field. 
      const RealVect pos  = probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);							      
      const RealVect E    = RealVect(D_DECL(electricFieldPatchEB(vof,0), electricFieldPatchEB(vof,1), electricFieldPatchEB(vof,2)));

      // Construct the CDR densities on the EB -- it needs to be in a form understandable by CdrPlasmaPhysics
      for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	const int  idx  = solverIt.index();
	const Real phi = (*a_cdrDensitiesEB[idx])[dit()](vof, comp);
	
	cdrDensitiesEB[idx] = std::max(zero, phi);
      }

      // Call our nifty plasma kinetics framework for computing the diffusion coefficients. 
      const Vector<Real> Dcos = m_physics->computeCdrDiffusionCoefficients(a_time,
									     pos,
									     E,
									     cdrDensitiesEB);

      // CdrPlasmaPhysics returns the diffusion coefficients in the definition order -- we need to put the coefficients
      // in the solver data holders now. 
      for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	const int idx = solverIt.index();

	// The check here is important because some of the input data will be nullptr if the species is not diffusive. 
	if(solverIt()->isDiffusive()){
	  (*a_cdrDcoEB[idx])[dit()](vof, comp) = Dcos[idx];
	}
      }
    };

    // Kernel region.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];

    // Launch the kernel.
    BoxLoops::loop(vofit, irregularKernel);
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocities(){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocities()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocities()" << endl;
  }

  // Compute the electric field first. 
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, m_cdr->getPhase(), SpaceDim);
  this->computeElectricField(E, m_cdr->getPhase(), m_fieldSolver->getPotential());

  // Get handle to CDR drift velocities and densities. 
  Vector<EBAMRCellData*>       cdrVelocities = m_cdr->getVelocities();
  const Vector<EBAMRCellData*> cdrDensities  = m_cdr->getPhis      ();  

  // Call the other version. 
  this->computeCdrDriftVelocities(cdrVelocities, cdrDensities, E, m_time);
}

void CdrPlasmaStepper::computeCdrDriftVelocities(Vector<EBAMRCellData*>&       a_cdrVelocities,
						 const Vector<EBAMRCellData*>& a_cdrDensities,
						 const EBAMRCellData&          a_electricField,
						 const Real&                   a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocities(Vector<EBAMRCellData*>x2, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocities(Vector<EBAMRCellData*>x2, EBAMRCellData, Real)" << endl;
  }

  // Number of CDR solvers/species.
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();

  CH_assert(a_electricField[0]->nComp() == SpaceDim     );
  CH_assert(a_cdrVelocities.size()      == numCdrSpecies);  
  CH_assert(a_cdrDensities. size()      == numCdrSpecies);

  // Level loop -- we will just call the level version. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    // Vectors for exposing the data on each level. 
    Vector<LevelData<EBCellFAB>* > cdrVelocities(numCdrSpecies, nullptr);
    Vector<LevelData<EBCellFAB>* > cdrDensities (numCdrSpecies, nullptr);

    // Fetch data on each level. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      cdrVelocities[idx] = (*a_cdrVelocities[idx])[lvl];
      cdrDensities [idx] = (*a_cdrDensities [idx])[lvl];
    }

    // Call the level version.
    this->computeCdrDriftVelocities(cdrVelocities, cdrDensities, *a_electricField[lvl], lvl, a_time);
  }

  // Coarsen and update ghost cells
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    if(solverIt()->isMobile()){
      m_amr->averageDown  (*a_cdrVelocities[idx], m_realm, m_cdr->getPhase()); 
      m_amr->interpGhostMG(*a_cdrVelocities[idx], m_realm, m_cdr->getPhase());
    }
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocities(Vector<LevelData<EBCellFAB> *>&       a_cdrVelocities,
						 const Vector<LevelData<EBCellFAB> *>& a_cdrDensities,
						 const LevelData<EBCellFAB> &          a_electricField,
						 const int                             a_lvl,
						 const Real&                           a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocities(Vector<LD<EBCellFAB>*>x2, LD<EBCellFAB>, int, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocities(Vector<LD<EBCellFAB>*>x2, LD<EBCellFAB>, int, Real)" << endl;
  }

  const int numCdrSpecies = m_physics->getNumCdrSpecies();  

  CH_assert(a_electricField.nComp() == SpaceDim     );
  CH_assert(a_cdrVelocities. size() == numCdrSpecies);
  CH_assert(a_cdrDensities.  size() == numCdrSpecies);    

  // Grid level and resolution
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const Real dx                 = m_amr->getDx()          [a_lvl];

  // Grid loop
  for (DataIterator dit(dbl); dit.ok(); ++dit){

    // Computational region. 
    const Box cellBox = dbl[dit()];

    // Data holding velocities and densities on each grid patch. 
    Vector<EBCellFAB*> cdrVelocities(numCdrSpecies, nullptr);
    Vector<EBCellFAB*> cdrDensities (numCdrSpecies, nullptr);

    Vector<FArrayBox*> cdrVelocitiesFAB(numCdrSpecies, nullptr);
    Vector<FArrayBox*> cdrDensitiesFAB (numCdrSpecies, nullptr);

    // Data holding the electric field. 
    const EBCellFAB& electricField    = a_electricField[dit()];
    const FArrayBox& electricFieldFAB = electricField.getFArrayBox();

    // Populate the patch data. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // If the species is mobile we fetch the data holding it's cell centered velocity. If it's not, the data holder
      // will be nullptr (because the solver did not allocate the data). 
      if(solverIt()->isMobile()){
	cdrVelocities   [idx] = &(*a_cdrVelocities[idx])[dit()];
	cdrVelocitiesFAB[idx] = &(   cdrVelocities[idx]->getFArrayBox());
      }

      // Densities -- this is always allocated. 
      cdrDensities   [idx] = &(*a_cdrDensities[idx])[dit()];
      cdrDensitiesFAB[idx] = &(   cdrDensities[idx]->getFArrayBox());      
    }

    // Call the regular version. This will also do irregular cells but we have to redo those below (because there might be multicells here, also).
    this->computeCdrDriftVelocitiesRegular(cdrVelocitiesFAB,
    					   cdrDensitiesFAB,
    					   electricFieldFAB,
    					   cellBox,
    					   a_time,
    					   dx);

    // Do the irregular and multi-valued cells. 
    this->computeCdrDriftVelocitiesIrregular(cdrVelocities,
					     cdrDensities,
					     electricField,
					     a_time,
					     dx,
					     a_lvl,
					     dit());

    // Velocities in covered cells are always bogus.
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // If it's mobile, set all velocity components to zero in the cut-cells. 
      if(solverIt()->isMobile()){
	for (int dir = 0; dir < SpaceDim; dir++){
	  cdrVelocities[idx]->setCoveredCellVal(0.0, dir);
	}	
      }
    }
  }
}

void CdrPlasmaStepper::computeCdrDriftVelocitiesRegular(Vector<FArrayBox*>&       a_cdrVelocities,
							const Vector<FArrayBox*>& a_cdrDensities,
							const FArrayBox&          a_electricField,
							const Box&                a_cellBox,
							const Real&               a_time,
							const Real&               a_dx){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocitiesRegular(Vector<FArrayBox*>x2, FArrayBox, Box, Real, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocitiesRegular(Vector<FArrayBox*>x2, FArrayBox, Box, Real, Real)" << endl;
  }

  // TLDR: This code will go through all single-valued grid cells and compute the drift velocities from our nifty
  //       little plasma physics framework. 

  constexpr int  comp = 0  ;
  constexpr Real zero = 0.0;  

  const int numCdrSpecies = m_physics->getNumCdrSpecies();  

  CH_assert(a_cdrVelocities. size() == numCdrSpecies);  
  CH_assert(a_cdrDensities.  size() == numCdrSpecies);        
  CH_assert(a_electricField.nComp() == SpaceDim     );

  // Lower-left corner in physical coordinates. 
  const RealVect probLo  = m_amr->getProbLo();

  // To be populated in each grid cell. 
  Vector<Real> cdrDensities(numCdrSpecies, 0.0);

  // Regular kernel. This will fill the velocities in all the input cells. 
  auto regularKernel = [&](const IntVect& iv) -> void {

    // Physical position and electric field
    const RealVect pos = probLo + (0.5*RealVect::Unit + RealVect(iv)) * a_dx;
    const RealVect E   = RealVect(D_DECL(a_electricField(iv, 0), a_electricField(iv, 1), a_electricField(iv, 2)));

    // Get all the densities in this grid cell. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx  = solverIt.index();
      const Real phi = (*a_cdrDensities[idx])(iv, comp);
      
      cdrDensities[idx] = std::max(phi, zero);
    }

    // Compute velocities
    const Vector<RealVect> velocities = m_physics->computeCdrDriftVelocities(a_time, pos, E, cdrDensities);

    // Put velocities in the appropriate place. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();

      // The check here is important because if the solver was not mobile then the memory for the velocity data holder
      // was not allocated either. 
      if(solverIt()->isMobile()){
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_cdrVelocities[idx])(iv, dir) = velocities[idx][dir];
	}
      }
    }
  };

  // Launch the kernel
  BoxLoops::loop(a_cellBox, regularKernel);
}

void CdrPlasmaStepper::computeCdrDriftVelocitiesIrregular(Vector<EBCellFAB*>&       a_cdrVelocities,
							  const Vector<EBCellFAB*>& a_cdrDensities,
							  const EBCellFAB&          a_electricField,
							  const Real&               a_time,
							  const Real&               a_dx,
							  const int                 a_lvl,
							  const DataIndex&          a_dit){
  CH_TIME("CdrPlasmaStepper::computeCdrDriftVelocitiesIrregular(Vector<EBCellFAB*>x2, EBCellFAB, Real, Real, int, DataIndex)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDriftVelocitiesIrregular(Vector<EBCellFAB*>x2, EBCellFAB, Real, Real, int, DataIndex)" << endl;
  }

  // TLDR: This will compute the drift velocities in each irregular grid cell. This also includes multi-valued cells if they are there.

  constexpr int comp  = 0  ;
  constexpr Real zero = 0.0;  

  // Number of CDR solvers. 
  const int numCdrSpecies = m_physics->getNumCdrSpecies();  

  CH_assert(a_cdrVelocities. size() == numCdrSpecies);  
  CH_assert(a_cdrDensities.  size() == numCdrSpecies);        
  CH_assert(a_electricField.nComp() == SpaceDim     );

  // Lower-left corner in physical coordinates. 
  const RealVect probLo  = m_amr->getProbLo();

  // EB box
  const EBISBox& ebisBox = a_electricField.getEBISBox();

  // CDR densities -- will be populated in each grid cell
  Vector<Real> cdrDensities(numCdrSpecies, 0.0);

  // Irregular kernel. 
  auto irregularKernel = [&](const VolIndex& vof) -> void {

    // Position and electric field
    const RealVect pos  = probLo + Location::position(Location::Cell::Center, vof, ebisBox, a_dx);        
    const RealVect E    = RealVect(D_DECL(a_electricField(vof, 0), a_electricField(vof, 1), a_electricField(vof, 2)));

    // Get CDR densities
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int  idx = solverIt.index();
      const Real phi = (*a_cdrDensities[idx])(vof, comp);
      
      cdrDensities[idx] = std::max(zero, phi);
    }
    
    // Plasma physics framework computes the velocities. 
    const Vector<RealVect> velocities = m_physics->computeCdrDriftVelocities(a_time, pos, E, cdrDensities);

    // Put velocities in the appropriate place. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      if(solverIt()->isMobile()){
	const int idx = solverIt.index();

	// Set each velocity component. 
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_cdrVelocities[idx])(vof, dir) = velocities[idx][dir];
	}
      }
    }
  };

  // Kernel region
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];

  // Launch the kernel
  BoxLoops::loop(vofit, irregularKernel);
}

void CdrPlasmaStepper::computeCdrFluxes(Vector<EBAMRIVData*>&       a_cdrFluxesEB,
					const Vector<EBAMRIVData*>& a_extrapCdrFluxes,
					const Vector<EBAMRIVData*>& a_extrapCdrDensities,
					const Vector<EBAMRIVData*>& a_extrapCdrVelocities,
					const Vector<EBAMRIVData*>& a_extrapCdrGradients,
					const Vector<EBAMRIVData*>& a_extrapRteFluxes,
					const EBAMRIVData&          a_electricField,
					const Real&                 a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrFluxes(Vector<EBAMRIVData*>x6, EBAMRIVData, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrFluxes(Vector<EBAMRIVData*>x6, EBAMRIVData, Real)" << endl;
  }

  // TLDR: This is the coupling function that computes the CDR fluxes on EB grid cells. It will call the level version. 

  // Number of CDR and RTE species. 
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();
  const int numRteSpecies  = m_physics->getNumRtSpecies();

  // Basic assertions. 
  CH_assert(a_cdrFluxesEB.        size() == numCdrSpecies);
  CH_assert(a_extrapCdrFluxes.    size() == numCdrSpecies);
  CH_assert(a_extrapCdrDensities. size() == numCdrSpecies);    
  CH_assert(a_extrapCdrVelocities.size() == numCdrSpecies);
  CH_assert(a_extrapCdrGradients. size() == numCdrSpecies);
  CH_assert(a_extrapRteFluxes.    size() == numRteSpecies);
  CH_assert(a_electricField[0]->nComp()  == SpaceDim     );

  // Go through each grid level. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    Vector<LevelData<BaseIVFAB<Real> >* > cdrFluxes          (numCdrSpecies, nullptr);
    Vector<LevelData<BaseIVFAB<Real> >* > extrapCdrFluxes    (numCdrSpecies, nullptr);
    Vector<LevelData<BaseIVFAB<Real> >* > extrapCdrDensities (numCdrSpecies, nullptr);
    Vector<LevelData<BaseIVFAB<Real> >* > extrapCdrVelocities(numCdrSpecies, nullptr);
    Vector<LevelData<BaseIVFAB<Real> >* > extrapCdrGradients (numCdrSpecies, nullptr);
    Vector<LevelData<BaseIVFAB<Real> >* > extrapRteFluxes    (numRteSpecies, nullptr);

    // Populate the level data for the CDR quantities. 
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      cdrFluxes          [idx] = (*a_cdrFluxesEB        [idx])[lvl];
      extrapCdrFluxes    [idx] = (*a_extrapCdrFluxes    [idx])[lvl];
      extrapCdrDensities [idx] = (*a_extrapCdrDensities [idx])[lvl];
      extrapCdrVelocities[idx] = (*a_extrapCdrVelocities[idx])[lvl];
      extrapCdrGradients [idx] = (*a_extrapCdrGradients [idx])[lvl];
    }

    // Populate the level data for the RTE quantities. 
    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      
      extrapRteFluxes[idx] = (*a_extrapRteFluxes[idx])[lvl];
    }

    // Now call the level version
    this->computeCdrFluxes(cdrFluxes,
			   extrapCdrFluxes,
			   extrapCdrDensities,
			   extrapCdrVelocities,
			   extrapCdrGradients,
			   extrapRteFluxes,
			   *a_electricField[lvl],
			   a_time,
			   lvl);
  }
}

void CdrPlasmaStepper::computeCdrFluxes(Vector<LevelData<BaseIVFAB<Real> >*>&       a_cdrFluxesEB,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrapCdrFluxes,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrapCdrDensities,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrapCdrVelocities,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrapCdrGradients,
					const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrapRteFluxes,
					const LevelData<BaseIVFAB<Real> >&          a_electricField,
					const Real&                                 a_time,
					const int                                   a_lvl){
  CH_TIME("CdrPlasmaStepper::computeCdrFluxes(Vector<LD<BaseIVFAB<Real> >*>x6, LD<BaseIVFAB<Real> >, Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrFluxes(Vector<LD<BaseIVFAB<Real> >*>x6, LD<BaseIVFAB<Real> >, Real)" << endl;
  }

  constexpr int comp  = 0;  

  // Number of CDR and RTE species.   
  const int numCdrSpecies  = m_physics->getNumCdrSpecies();
  const int numRteSpecies  = m_physics->getNumRtSpecies();

  // Basic assertions.
  CH_assert(a_cdrFluxesEB.        size() == numCdrSpecies);
  CH_assert(a_extrapCdrFluxes.    size() == numCdrSpecies);
  CH_assert(a_extrapCdrDensities. size() == numCdrSpecies);    
  CH_assert(a_extrapCdrVelocities.size() == numCdrSpecies);
  CH_assert(a_extrapCdrGradients. size() == numCdrSpecies);
  CH_assert(a_extrapRteFluxes.    size() == numRteSpecies);
  CH_assert(a_electricField.     nComp() == SpaceDim     );

  // Things that will be passed into our nifty little plasma physics framework. 
  Vector<Real> extrapCdrFluxes    (numCdrSpecies, 0.0);
  Vector<Real> extrapCdrDensities (numCdrSpecies, 0.0);
  Vector<Real> extrapCdrVelocities(numCdrSpecies, 0.0);
  Vector<Real> extrapCdrGradients (numCdrSpecies, 0.0);
  Vector<Real> extrapRteFluxes    (numRteSpecies, 0.0);

  // Lower left corner (physical coordinates). 
  const RealVect probLo = m_amr->getProbLo();  

  // Grids and EB information on this level. 
  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const Real               dx    = m_amr->getDx()[a_lvl];
  const MFLevelGrid&       mflg  = *(m_amr->getMFLevelGrid(m_realm)[a_lvl]);

  // Patch loop
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box&     cellBox = dbl  [dit()];
    const EBISBox& ebisBox = ebisl[dit()];
    const EBGraph& ebgraph = ebisBox.getEBGraph();

    // These are the dielectric and electrode cells in the current grid patch. They form a complete set of
    // cut-cells. 
    const IntVectSet& dielectricCells = mflg.interfaceRegion(cellBox, dit());
    const IntVectSet& electrodeCells  = ebisBox.getIrregIVS(cellBox) - dielectricCells;

    // Electric field in this patch. 
    const BaseIVFAB<Real>& electricField = a_electricField[dit()];    

    // Kernel for generating boundary conditions on an electrode cut-cell. This kernel will call the CdrPlasmaPhysics routine
    // which computes the electrode EB boundary fluxes. 
    auto electrodeKernel = [&](const VolIndex& vof) -> void {

      // Position, normal vector, and electric field.
      const RealVect pos    = probLo + Location::position(Location::Cell::Boundary, vof, ebisBox, dx);      
      const RealVect normal = ebisBox.normal(vof);
      const RealVect E      = RealVect(D_DECL(electricField(vof,0), electricField(vof,1), electricField(vof,2)));


      // Populate the CDR related quantities. 
      for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	const int idx = solverIt.index();
	
	extrapCdrFluxes    [idx] = (*a_extrapCdrFluxes    [idx])[dit()](vof,comp);
	extrapCdrDensities [idx] = (*a_extrapCdrDensities [idx])[dit()](vof,comp);
	extrapCdrVelocities[idx] = (*a_extrapCdrVelocities[idx])[dit()](vof,comp);
	extrapCdrGradients [idx] = (*a_extrapCdrGradients [idx])[dit()](vof,comp);
      }

      // Populate the RTE related quantities. 
      for (RtIterator<RtSolver> solverIt(*m_rte); solverIt.ok(); ++solverIt){
	const int idx = solverIt.index();
	
	extrapRteFluxes[idx] = (*a_extrapRteFluxes[idx])[dit()](vof,comp);
      }

      // Compute the CDR fluxes using our framework. 
      const Vector<Real> cdrFluxes = m_physics->computeCdrElectrodeFluxes(a_time,
									  pos,
									  normal,
									  E,
									  extrapCdrDensities,
									  extrapCdrVelocities,
									  extrapCdrGradients,
									  extrapRteFluxes,
									  extrapCdrFluxes);
	
      // Put the fluxes in their respective place
      for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	const int idx = solverIt.index();
	(*a_cdrFluxesEB[idx])[dit()](vof, comp) = cdrFluxes[idx];
      }
    };

    // Kernel for generating boundary conditions on a dielectric cut-cell. This kernel will call the CdrPlasmaPhysics routine
    // which computes the dielectric EB boundary fluxes. 
    auto dielectricKernel = [&](const VolIndex& vof) -> void {

      // Position, normal vector, and electric field.
      const RealVect pos    = probLo + Location::position(Location::Cell::Boundary, vof, ebisBox, dx);      
      const RealVect normal = ebisBox.normal(vof);
      const RealVect E      = RealVect(D_DECL(electricField(vof,0), electricField(vof,1), electricField(vof,2)));


      // Populate the CDR related quantities. 
      for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	const int idx = solverIt.index();
	
	extrapCdrFluxes    [idx] = (*a_extrapCdrFluxes    [idx])[dit()](vof,comp);
	extrapCdrDensities [idx] = (*a_extrapCdrDensities [idx])[dit()](vof,comp);
	extrapCdrVelocities[idx] = (*a_extrapCdrVelocities[idx])[dit()](vof,comp);
	extrapCdrGradients [idx] = (*a_extrapCdrGradients [idx])[dit()](vof,comp);
      }

      // Populate the RTE related quantities. 
      for (RtIterator<RtSolver> solverIt(*m_rte); solverIt.ok(); ++solverIt){
	const int idx = solverIt.index();
	
	extrapRteFluxes[idx] = (*a_extrapRteFluxes[idx])[dit()](vof,comp);
      }

      // Compute the CDR fluxes using our framework. 
      const Vector<Real> cdrFluxes = m_physics->computeCdrDielectricFluxes(a_time,
									   pos,
									   normal,
									   E,
									   extrapCdrDensities,
									   extrapCdrVelocities,
									   extrapCdrGradients,
									   extrapRteFluxes,
									   extrapCdrFluxes);
	
      // Put the fluxes in their respective place
      for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	const int idx = solverIt.index();
	(*a_cdrFluxesEB[idx])[dit()](vof, comp) = cdrFluxes[idx];
      }
    };    

    // Kernel regions.
    VoFIterator vofitElectrode (electrodeCells,  ebgraph);
    VoFIterator vofitDielectric(dielectricCells, ebgraph);        

    // Launch the kernels
    BoxLoops::loop(vofitElectrode,  electrodeKernel );
    BoxLoops::loop(vofitDielectric, dielectricKernel);
  }
}



void CdrPlasmaStepper::computeCdrDomainFluxes(Vector<EBAMRIFData*>&       a_fluxes,
					      const Vector<EBAMRIFData*>& a_extrapCdrFluxes,
					      const Vector<EBAMRIFData*>& a_extrapCdrDensities,
					      const Vector<EBAMRIFData*>& a_extrapCdrVelocities,
					      const Vector<EBAMRIFData*>& a_extrapCdrGradients,
					      const Vector<EBAMRIFData*>& a_extrapRteFluxes,
					      const EBAMRIFData&          a_E,
					      const Real&                 a_time){
  CH_TIME("CdrPlasmaStepper::computeCdrDomainFluxes(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDomainFluxes(level)" << endl;
  }

  const int numCdrSpecies  = m_physics->getNumCdrSpecies();
  const int numRteSpecies  = m_physics->getNumRtSpecies();
  const int finest_level = m_amr->getFinestLevel();

  // Things that will be passed into physics
  Vector<Real> extrapCdrFluxes(numCdrSpecies);
  Vector<Real> extrapCdrDensities(numCdrSpecies);
  Vector<Real> extrapCdrVelocities(numCdrSpecies);
  Vector<Real> extrapCdrGradients(numCdrSpecies);
  Vector<Real> extrapRteFluxes(numRteSpecies);

  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<LevelData<DomainFluxIFFAB>* > fluxes(numCdrSpecies);
    Vector<LevelData<DomainFluxIFFAB>* > extrapCdrFluxes(numCdrSpecies);
    Vector<LevelData<DomainFluxIFFAB>* > extrapCdrDensities(numCdrSpecies);
    Vector<LevelData<DomainFluxIFFAB>* > extrapCdrVelocities(numCdrSpecies);
    Vector<LevelData<DomainFluxIFFAB>* > extrapCdrGradients(numCdrSpecies);
    Vector<LevelData<DomainFluxIFFAB>* > extrapRteFluxes(numRteSpecies);
    
    for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      fluxes[idx]                = (*a_fluxes[idx])[lvl];
      extrapCdrFluxes[idx]     = (*a_extrapCdrFluxes[idx])[lvl];
      extrapCdrDensities[idx]  = (*a_extrapCdrDensities[idx])[lvl];
      extrapCdrVelocities[idx] = (*a_extrapCdrVelocities[idx])[lvl];
      extrapCdrGradients[idx]  = (*a_extrapCdrGradients[idx])[lvl];
    }

    for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
      const int idx = solverIt.index();
      extrapRteFluxes[idx] = (*a_extrapRteFluxes[idx])[lvl];
    }

    // Call the level version
    computeCdrDomainFluxes(fluxes, extrapCdrFluxes, extrapCdrDensities, extrapCdrVelocities,
			   extrapCdrGradients, extrapRteFluxes, *a_E[lvl], a_time, lvl);
    
  }
}

void CdrPlasmaStepper::computeCdrDomainFluxes(Vector<LevelData<DomainFluxIFFAB>*>        a_fluxes,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrFluxes,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrDensities,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrVelocities,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrGradients,
					      const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapRteFluxes,
					      const LevelData<DomainFluxIFFAB>&          a_E,
					      const Real&                                a_time,
					      const int                                  a_lvl){
  CH_TIME("CdrPlasmaStepper::computeCdrDomainFluxes(level)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeCdrDomainFluxes(level)" << endl;
  }

  const int numCdrSpecies  = m_physics->getNumCdrSpecies();
  const int numRteSpecies  = m_physics->getNumRtSpecies();
  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();

  // Things that will be passed into physics
  Vector<Real> extrapCdrFluxes(numCdrSpecies);
  Vector<Real> extrapCdrDensities(numCdrSpecies);
  Vector<Real> extrapCdrVelocities(numCdrSpecies);
  Vector<Real> extrapCdrGradients(numCdrSpecies);
  Vector<Real> extrapRteFluxes(numRteSpecies);

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl       = m_amr->getEBISLayout(m_realm, m_cdr->getPhase())[a_lvl];
  const ProblemDomain& domain   = m_amr->getDomains()[a_lvl];
  const Real dx                 = m_amr->getDx()[a_lvl];
  const RealVect probLo         = m_amr->getProbLo();

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
	  const RealVect pos    = EBArith::getFaceLocation(face, dx*RealVect::Unit, probLo);
	    
	  // Define the electric field
	  const RealVect E = RealVect(D_DECL((a_E)[dit()](dir, sit())(face,0),
					     (a_E)[dit()](dir, sit())(face,1),
					     (a_E)[dit()](dir, sit())(face,2)));

	  // Ion densities. 
	  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	    const int idx = solverIt.index();
	    extrapCdrFluxes[idx]     = (*a_extrapCdrFluxes[idx])[dit()](dir, sit())(face,comp);
	    extrapCdrDensities[idx]  = (*a_extrapCdrDensities[idx])[dit()](dir, sit())(face,comp);
	    extrapCdrVelocities[idx] = (*a_extrapCdrVelocities[idx])[dit()](dir, sit())(face,comp);
	    extrapCdrGradients[idx]  = (*a_extrapCdrGradients[idx])[dit()](dir, sit())(face,comp);
	  }

	  // Photon fluxes
	  for (RtIterator<RtSolver> solverIt(*m_rte); solverIt.ok(); ++solverIt){
	    const int idx = solverIt.index();
	    extrapRteFluxes[idx] = (*a_extrapRteFluxes[idx])[dit()](dir, sit())(face,comp);
	  }

	  // Call CdrPlasmaPhysics
	  const Vector<Real> fluxes = m_physics->computeCdrDomainFluxes(a_time,
									pos,
									dir,
									sit(),
									E,
									extrapCdrDensities,
									extrapCdrVelocities,
									extrapCdrGradients,
									extrapRteFluxes,
									extrapCdrFluxes);

	  // Put fluxes where they belong
	  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
	    const int idx = solverIt.index();
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

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const int idx = solverIt.index();
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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const int idx = solverIt.index();
    if(solver->isMobile()){
      extrapolateToVectorDomainFaces(*a_extrap[idx], a_phase, *a_velocities[idx]);
    }
    else{
      DataOps::setValue(*a_extrap[idx], 0.0);
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


void CdrPlasmaStepper::computeJ(EBAMRCellData& a_J) const{
  CH_TIME("CdrPlasmaStepper::computeJ(amr)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeJ(amr)" << endl;
  }

  CH_assert(a_J[0]->nComp() == SpaceDim);

  DataOps::setValue(a_J, 0.0);

  // Allocate scratch storage for holding the cell-centered diffusion coefficient and
  // D * grad(phi)
  EBAMRCellData scratchONE;
  EBAMRCellData scratchDIM;

  m_amr->allocate(scratchONE, m_realm, m_phase, 1       );
  m_amr->allocate(scratchDIM, m_realm, m_phase, SpaceDim);  
  
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solverIt.getSpecies();
    
    const int  Z           = species->getChargeNumber();
    const bool isMobile    = solver->isMobile();
    const bool isDiffusive = solver->isDiffusive();

    // Cell-centered density
    const EBAMRCellData& phi = solver->getPhi();

    // Add the drift contribution to the current. We compute v*n*Z and add it to J
    if(Z != 0 && isMobile){

      // Cell-centered velocity
      const EBAMRCellData& vel = solver->getCellCenteredVelocity();

      // Now compute vel*phi and add the contribution to the current density. 
      DataOps::copy(scratchDIM, vel);
      DataOps::multiplyScalar(scratchDIM, phi);
      DataOps::incr(a_J, scratchDIM, Real(Z));
    }

    // Add the diffusive contribution to the current.
    if(Z != 0 && isDiffusive){

      // We need updated ghost cells when computing the gradient so we use scratchONE as a scratch data holder when we compute grad(phi)
      DataOps::copy(scratchONE, phi);

      m_amr->averageDown    (scratchONE,             m_realm, m_phase);
      m_amr->interpGhostMG  (scratchONE,             m_realm, m_phase);
      m_amr->computeGradient(scratchDIM, scratchONE, m_realm, m_phase);

      // scratchONE now holds grad(phi). We need to put the diffusion coefficient on the cell center now.
      const EBAMRFluxData& diffusionCoefficientFace = solver->getFaceCenteredDiffusionCoefficient();

      DataOps::averageFaceToCell(scratchONE, diffusionCoefficientFace, m_amr->getDomains());

      // Now make DgradPhi = D * grad(Phi)
      DataOps::multiplyScalar(scratchDIM, scratchONE);

      // Add contribution -Z * D * grad(Phi) to the current
      DataOps::incr(a_J, scratchDIM, -Real(Z));
    }
  }

  // Now scale by electron charge.
  DataOps::scale(a_J, Units::Qe);

  // Coarsen and update ghost cells
  m_amr->averageDown(a_J, m_realm, m_phase);
  m_amr->interpGhost(a_J, m_realm, m_phase);  
}

void CdrPlasmaStepper::computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(MFAMRCellData, MFAMRCellData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(MFAMRCellData, MFAMRCellData)" << endl;
  }


  CH_assert(a_E[0]        ->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1       );  

  // Field solver computes. 
  m_fieldSolver->computeElectricField(a_E, a_potential);

  // Update ghost cells. 
  m_amr->averageDown  (a_E, m_realm);
  m_amr->interpGhostMG(a_E, m_realm);
}

void CdrPlasmaStepper::computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(EBAMRCellData, phase)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(EBAMRCellData, phase)" << endl;
  }

  CH_assert(a_E[0]->nComp() == SpaceDim);

  this->computeElectricField(a_E, a_phase, m_fieldSolver->getPotential());
}

void CdrPlasmaStepper::computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(EBAMRCellData, phase, MFAMRCellData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(EBAMRCellData, phase, MFAMRCellData)" << endl;
  }

  CH_assert(a_E[0]        ->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1       );    

  m_fieldSolver->computeElectricField(a_E, a_phase, a_potential);

  m_amr->averageDown  (a_E, m_realm, a_phase);
  m_amr->interpGhostMG(a_E, m_realm, a_phase);
}

void CdrPlasmaStepper::computeElectricField(EBAMRFluxData& a_electricFieldFace, const phase::which_phase a_phase, const EBAMRCellData& a_electricFieldCell) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(EBAMRFluxData, phase, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(EBAMRFluxData, phase, EBAMRCellData)" << endl;
  }

  // TLDR: This will compute the electric field on grid faces as an arithmetic average
  //       of the cell-centered field. 

  CH_assert(a_electricFieldFace[0]->nComp() == SpaceDim);
  CH_assert(a_electricFieldCell[0]->nComp() == SpaceDim);

  // Grid loop.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    // Patch distribution and EB information on this level.
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, a_phase)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];

    // Patch loop
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const EBCellFAB& electricFieldCell = (*a_electricFieldCell[lvl])[dit()];
      const EBISBox&   ebisbox  = ebisl[dit()];
      const EBGraph&   ebgraph  = ebisbox.getEBGraph();
      const Box&       cellBox          = dbl.get(dit());

      // Do faces in all directions. 
      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& electricFieldFace = (*a_electricFieldFace[lvl])[dit()][dir];
	
	electricFieldFace.setVal(0.0);

	// Average to faces. 
	EBLevelDataOps::averageCellToFace(electricFieldFace,
					  electricFieldCell,
					  ebgraph,
					  cellBox,
					  0,
					  dir,
					  domain,
					  dir,
					  dir);
      }

    }
    //    DataOps::averageCellToFace(*a_electricFieldFace[lvl], *a_electricFieldCell[lvl], m_amr->getDomains()[lvl]);
    a_electricFieldFace[lvl]->exchange();
  }
}

void CdrPlasmaStepper::computeElectricField(EBAMRIVData& a_electricFieldEB, const phase::which_phase a_phase, const EBAMRCellData& a_electricFieldCell) const {
  CH_TIME("CdrPlasmaStepper::computeElectricField(EBAMRIVData, phase, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeElectricField(EBAMRIVData, phase, EBAMRCellData)" << endl;
  }

  CH_assert(a_electricFieldEB  [0]->nComp() == SpaceDim);
  CH_assert(a_electricFieldCell[0]->nComp() == SpaceDim);

  // Interpolate to the EB centroid
  const IrregAmrStencil<EbCentroidInterpolationStencil>& interpStencils = m_amr->getEbCentroidInterpolationStencils(m_realm, a_phase);
  interpStencils.apply(a_electricFieldEB, a_electricFieldCell);
}

void CdrPlasmaStepper::computeMaxElectricField(Real& a_maximumElectricField, const phase::which_phase a_phase){
  CH_TIME("CdrPlasmaStepper::computeMaxElectricField(Real, phase)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeMaxElectricField(Real, phase)" << endl;
  }

  // Compute the electric field on the input phase and interpolate it to centroids. 
  EBAMRCellData E;
  m_amr->allocate(E, m_realm, a_phase, SpaceDim);

  this->computeElectricField(E, a_phase, m_fieldSolver->getPotential());
  m_amr->interpToCentroids(E, m_realm, a_phase);

  // Get the maximum and minimum values of |E|. 
  Real max = -std::numeric_limits<Real>::max();
  Real min =  std::numeric_limits<Real>::max();
  
  DataOps::getMaxMinNorm(max, min, E);

  a_maximumElectricField = max;
}

void CdrPlasmaStepper::computeChargeFlux(EBAMRIVData& a_flux, Vector<EBAMRIVData*>& a_cdr_fluxes){
  CH_TIME("CdrPlasmaStepper::computeChargeFlux");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::computeChargeFlux" << endl;
  }

  MayDay::Abort("CdrPlasmaStepper::computeChargeFlux - I'm suspecting that this is deprecated code which is no longer used");

  DataOps::setValue(a_flux, 0.0);

  for (CdrIterator<CdrSolver> solverIt(*m_cdr); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const RefCountedPtr<CdrSpecies>& spec      = solverIt.getSpecies();
    const EBAMRIVData& solver_flux          = *a_cdr_fluxes[solverIt.index()];

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

  const IrregAmrStencil<EbCentroidInterpolationStencil>& interpStencils = m_amr->getEbCentroidInterpolationStencils(m_realm, a_phase);

  //  for (int i = 0; i < a_fluxes.size(); i++){
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver = solverIt();
    const int idx = solverIt.index();
    if(solver->isMobile()){
      interpStencils.apply(eb_vel, *a_velocities[idx]);
      interpStencils.apply(eb_phi, *a_densities[idx]);

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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const int idx = solverIt.index();
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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const int idx = solverIt.index();
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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver = solverIt();
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

  const RealVect probLo  = m_amr->getProbLo();
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
	const RealVect pos  = probLo + vof.gridIndex()*dx + ebisbox.bndryCentroid(vof)*dx;
	
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

  this->solveRadiativeTransfer(states, rhs, cdr_states, E, m_time, a_dt);
}

void CdrPlasmaStepper::solveRadiativeTransfer(Vector<EBAMRCellData*>&       a_rte_states,
					      Vector<EBAMRCellData*>&       a_rte_sources,
					      const Vector<EBAMRCellData*>& a_cdr_states,
					      const EBAMRCellData&          a_E,
					      const Real                    a_time,
					      const Real                    a_dt) {
  CH_TIME("CdrPlasmaStepper::solveRadiativeTransfer(full)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaStepper::solveRadiativeTransfer(full)" << endl;
  }

  //  this->compute_rte_sources(a_rte_sources, a_cdr_states, a_E, a_time, a_centering);
  //  advanceReactionNetwork(a_time, a_dt);

  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    RefCountedPtr<RtSolver>& solver = solverIt();
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
  
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const RefCountedPtr<CdrSpecies>& spec      = solverIt.getSpecies();
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
  
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const RefCountedPtr<CdrSpecies>& spec      = solverIt.getSpecies();
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
  
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const RefCountedPtr<CdrSpecies>& spec      = solverIt.getSpecies();
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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  // RTE solvers checkpoint their data
  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<RtSolver>& solver = solverIt();
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

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver = solverIt();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<RtSolver>& solver = solverIt();
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
  
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver = solverIt();
    ncomp += solver->getNumberOfPlotVariables();
  }
  
  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<RtSolver>& solver = solverIt();
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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver = solverIt();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<RtSolver>& solver = solverIt();
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
