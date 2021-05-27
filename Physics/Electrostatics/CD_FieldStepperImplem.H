/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldStepperImplem.H
  @brief  Implementation of CD_FieldStepper.H
  @author Robert Marskar
*/

#ifndef CD_FieldStepperImplem_H
#define CD_FieldStepperImplem_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_FieldStepper.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace Electrostatics {
  
    template <typename T>
    Real FieldStepper<T>::s_constant_one(const RealVect a_pos){
      return 1.0;
    }

    template <typename T>
    Real FieldStepper<T>::s_potential_one(const Real a_time){
      return 1.0;
    }

    template <class T>
    FieldStepper<T>::FieldStepper(){
      m_verbosity = -1;

      ParmParse pp("FieldStepper");

      std::string str;

      pp.get("load_balance", m_loadBalance);
      pp.get("box_sorting",  str);
      pp.get("realm",        m_realm);
      pp.get("verbosity",    m_verbosity);
      
      if( str == "none"){
	m_boxSort = BoxSorting::None;
      }
      if( str == "std"){
	m_boxSort = BoxSorting::Std;
      }
      else if(str == "shuffle"){
	m_boxSort = BoxSorting::Shuffle;
      }
      else if(str == "morton"){
	m_boxSort = BoxSorting::Morton;
      }
      else {
	MayDay::Abort("FieldStepper::FieldStepper - unknown box sorting method requested for argument 'BoxSorting'");
      }
    }

    template <class T>
    FieldStepper<T>::~FieldStepper(){

    }

    template <class T>
    void FieldStepper<T>::setupSolvers(){
      
      // Poisson setup
      m_fieldSolver = RefCountedPtr<FieldSolver> (new T());
      m_fieldSolver->setVerbosity(m_verbosity);
      m_fieldSolver->parseOptions();
      m_fieldSolver->setAmr(m_amr);
      m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
      m_fieldSolver->setVoltage(FieldSolver::s_voltageOne);
      m_fieldSolver->setRealm(m_realm);
      m_fieldSolver->setTime(0, 0.0, 0.0);

      // Setup sigma
      m_sigma = RefCountedPtr<SigmaSolver> (new SigmaSolver());
      m_sigma->setVerbosity(m_verbosity);
      m_sigma->setAmr(m_amr);
      m_sigma->setComputationalGeometry(m_computationalGeometry);
      m_sigma->setRealm(m_realm);

    }

    template <class T>
    void FieldStepper<T>::registerOperators(){
      m_fieldSolver->registerOperators();
    }

    template <class T>
    void FieldStepper<T>::registerRealms(){
      m_amr->registerRealm(m_realm);
    }

    template <class T>
    void FieldStepper<T>::allocate() {
      m_fieldSolver->allocateInternals();
      m_sigma->allocateInternals();
    }

    template <class T>
    void FieldStepper<T>::initialData(){
      ParmParse pp("FieldStepper");
      pp.get("init_rho",   m_init_rho);
      pp.get("init_sigma", m_init_sigma);

      MFAMRCellData& state = m_fieldSolver->getPotential();  
      MFAMRCellData& rho   = m_fieldSolver->getRho();
      EBAMRIVData& sigma   = m_sigma->getPhi();

      DataOps::setValue(state, 0.0);
      DataOps::setValue(rho,   m_init_rho);
      DataOps::setValue(sigma, m_init_sigma);
    }

    template <class T>
    void FieldStepper<T>::solvePoisson(){
      
      MFAMRCellData& state = m_fieldSolver->getPotential();  
      MFAMRCellData& rho   = m_fieldSolver->getRho();
      EBAMRIVData& sigma   = m_sigma->getPhi();

      const bool converged = m_fieldSolver->solve(state, rho, sigma);
      
      m_fieldSolver->computeElectricField();

      if(!converged){
	MayDay::Warning("FieldStepper<T>::advance - did not converge");
      }
    }
    
    template <class T>
    void FieldStepper<T>::postInitialize() {
      this->solvePoisson();
    }

    template <class T>
    void FieldStepper<T>::postRegrid(){
      this->solvePoisson();
    }

    template<class T>
    Real FieldStepper<T>::advance(const Real a_dt){
      MayDay::Abort("FieldStepper<T>::advance - callling this is an error. Please set Driver.max_steps = 0");
      return 1.E99;
    }

    template <class T>
    void FieldStepper<T>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
      MayDay::Abort("FieldStepper<T>::readCheckpointData - checkpointing not supported for this class");
    }

    template <class T>
    int FieldStepper<T>::getNumberOfPlotVariables() const{
      int ncomp = m_fieldSolver->getNumberOfPlotVariables();
      ncomp += m_sigma->getNumberOfPlotVariables();
      return ncomp;
    }

    template <class T>
    void FieldStepper<T>::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const{
      a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
      m_fieldSolver->writePlotData(a_output, a_icomp);

      a_plotVariableNames.append(m_sigma->getPlotVariableNames());
      m_sigma->writePlotData(a_output, a_icomp);
    }

    template <class T>
    void FieldStepper<T>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt){
      m_timeStep = a_step;
      m_time = a_time;
      m_dt   = a_dt;

      m_fieldSolver->setTime(a_step, a_time, a_dt);
    }

    template <class T>
    void FieldStepper<T>::preRegrid(const int a_lbase, const int a_oldFinestLevel){
      m_fieldSolver->preRegrid(a_lbase, a_oldFinestLevel);
      m_sigma->preRegrid(a_lbase, a_oldFinestLevel);
    }

    template <class T>
    void FieldStepper<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
      m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
      m_sigma->regrid(a_lmin,   a_oldFinestLevel, a_newFinestLevel);

      // No solve here. That is done in postRegrid for this class. 
      //      m_fieldSolver->solve();
    }

    template <class T>
    bool FieldStepper<T>::loadBalanceThisRealm(const std::string a_realm) const {
      return m_loadBalance ? true : false;
    }

    template <class T>
    void FieldStepper<T>::loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
					   Vector<Vector<Box> >&            a_boxes,
					   const std::string                a_realm,
					   const Vector<DisjointBoxLayout>& a_grids,
					   const int                        a_lmin,
					   const int                        a_finestLevel){
      if(m_loadBalance){
	a_procs.resize(1+a_finestLevel);
	a_boxes.resize(1+a_finestLevel);

	// We need to make AmrMesh restore some operators that we need. Note that this is overkill because we only need
	// the multigrid ghost cell interpolation operator. But for simplicity, do all of them. 
	m_amr->regridOperators(a_lmin, a_finestLevel, 1);


	// Field solver implementation gets the responsibility of computing loads on each level. 
	for (int lvl = 0; lvl <= a_finestLevel; lvl++){
	  Vector<long long> loads = m_fieldSolver->computeLoads(a_grids[lvl], lvl);

	  // Do the desired sorting and load balancing
	  a_boxes[lvl] = a_grids[lvl].boxArray();
	  
	  LoadBalancing::sort(a_boxes[lvl], loads, m_boxSort);
	  LoadBalancing::makeBalance(a_procs[lvl], loads, a_boxes[lvl]);
	}
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
