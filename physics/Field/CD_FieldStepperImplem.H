/*!
  @file   FieldStepperI.H
  @brief  Implementation of the FieldStepper
  @author Robert Marskar
  @data   March 2020
*/

#include "CD_FieldStepper.H"

#include <ParmParse.H>
#include <EBEllipticLoadBalance.H>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace poisson {
  
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

      pp.get("load_balance", m_load_balance);
      pp.get("box_sorting", str);
      pp.get("realm",       m_realm);
      pp.get("verbosity",   m_verbosity);
      
      if( str == "none"){
	m_boxsort = box_sorting::none;
      }
      if( str == "std"){
	m_boxsort = box_sorting::std;
      }
      else if(str == "shuffle"){
	m_boxsort = box_sorting::shuffle;
      }
      else if(str == "morton"){
	m_boxsort = box_sorting::morton;
      }
      else {
	MayDay::Abort("FieldStepper::FieldStepper - unknown box sorting method requested for argument 'box_sorting'");
      }
    }

    template <class T>
    FieldStepper<T>::~FieldStepper(){

    }

    template <class T>
    void FieldStepper<T>::setup_solvers(){
      // Poisson setup
      m_poisson = RefCountedPtr<FieldSolver> (new T());
      m_poisson->set_verbosity(m_verbosity);
      m_poisson->parse_options();
      m_poisson->set_amr(m_amr);
      m_poisson->set_computational_geometry(m_compgeom);
      m_poisson->set_potential(s_potential_one);
      m_poisson->set_realm(m_realm);
      m_poisson->set_time(0, 0.0, 0.0);

      // Setup sigma
      m_sigma = RefCountedPtr<sigma_solver> (new sigma_solver());
      m_sigma->set_verbosity(m_verbosity);
      m_sigma->set_amr(m_amr);
      m_sigma->set_computational_geometry(m_compgeom);
      m_sigma->set_realm(m_realm);

    }

    template <class T>
    void FieldStepper<T>::register_operators(){
      m_poisson->register_operators();
    }

    template <class T>
    void FieldStepper<T>::register_realms(){
      m_amr->register_realm(m_realm);
    }

    template <class T>
    void FieldStepper<T>::allocate() {
      m_poisson->allocate_internals();
      m_sigma->allocate_internals();
    }

    template <class T>
    void FieldStepper<T>::initial_data(){
      ParmParse pp("FieldStepper");
      pp.get("init_rho",   m_init_rho);
      pp.get("init_sigma", m_init_sigma);

      MFAMRCellData& state = m_poisson->get_state();  
      MFAMRCellData& rho   = m_poisson->get_source();
      EBAMRIVData& sigma   = m_sigma->get_state();

      data_ops::set_value(state, 0.0);
      data_ops::set_value(rho,   m_init_rho);
      data_ops::set_value(sigma, m_init_sigma);
    }

    template <class T>
    void FieldStepper<T>::solve_poisson(){
      MFAMRCellData& state = m_poisson->get_state();  
      MFAMRCellData& rho   = m_poisson->get_source();
      EBAMRIVData& sigma   = m_sigma->get_state();

      const Real t0 = MPI_Wtime();
      const bool converged = m_poisson->solve(state, rho, sigma);
      const Real t1 = MPI_Wtime();
      // if(procID() == 0){
      // 	std::cout << t1 - t0 << std::endl;
      // }
      m_poisson->compute_E();

      if(!converged){
	MayDay::Warning("FieldStepper<T>::advance - did not converge");
      }
    }
    
    template <class T>
    void FieldStepper<T>::post_initialize() {
      this->solve_poisson();
    }

    template <class T>
    void FieldStepper<T>::post_regrid(){
      this->solve_poisson();
    }

    template<class T>
    Real FieldStepper<T>::advance(const Real a_dt){
      MayDay::Abort("FieldStepper<T>::advance - callling this is an error. Please set driver.max_steps = 0");
      return 1.E99;
    }

    template <class T>
    void FieldStepper<T>::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
      MayDay::Abort("FieldStepper<T>::read_checkpoint_data - checkpointing not supported for this class");
    }

    template <class T>
    int FieldStepper<T>::get_num_plot_vars() const{
      int ncomp = m_poisson->get_num_plotvars();
      ncomp += m_sigma->get_num_plotvars();
      return ncomp;
    }

    template <class T>
    void FieldStepper<T>::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const{
      a_plotvar_names.append(m_poisson->get_plotvar_names());
      m_poisson->write_plot_data(a_output, a_icomp);

      a_plotvar_names.append(m_sigma->get_plotvar_names());
      m_sigma->write_plot_data(a_output, a_icomp);
    }

    template <class T>
    void FieldStepper<T>::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
      m_step = a_step;
      m_time = a_time;
      m_dt   = a_dt;
    }

    template <class T>
    void FieldStepper<T>::pre_regrid(const int a_lbase, const int a_old_finest_level){
      m_poisson->pre_regrid(a_lbase, a_old_finest_level);
      m_sigma->pre_regrid(a_lbase, a_old_finest_level);
    }

    template <class T>
    void FieldStepper<T>::deallocate(){
      m_poisson->deallocate_internals();
      m_sigma->deallocate_internals();
    }

    template <class T>
    void FieldStepper<T>::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
      m_poisson->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
      m_sigma->regrid(a_lmin,   a_old_finest_level, a_new_finest_level);

      // No solve here. That is done in post_regrid for this class. 
      //      m_poisson->solve();
    }

    template <class T>
    bool FieldStepper<T>::load_balance_realm(const std::string a_realm) const {
      return m_load_balance ? true : false;
    }

    template <class T>
    void FieldStepper<T>::load_balance_boxes(Vector<Vector<int> >&            a_procs,
					      Vector<Vector<Box> >&            a_boxes,
					      const std::string                a_realm,
					      const Vector<DisjointBoxLayout>& a_grids,
					      const int                        a_lmin,
					      const int                        a_finest_level){
      if(m_load_balance){
	constexpr int num_apply = 50;
	
	a_procs.resize(1+a_finest_level);
	a_boxes.resize(1+a_finest_level);

	FieldSolverMultigrid* solver = (FieldSolverMultigrid*) (&(*m_poisson));

	MFAMRCellData dummy1, dummy2;

	m_amr->allocate(dummy1, m_realm, 1);
	m_amr->allocate(dummy2, m_realm, 1);

	const Vector<ProblemDomain>& domains = m_amr->get_domains();

	// We need to make amr_mesh restore some operators that we need. This is overkill because we only need
	// operators for ghost cell interpolation....

	m_amr->regrid_operators(a_lmin, a_finest_level, 1);

	// Set up the operator factor and fetch it.
	solver->set_coefficients();
	solver->setup_operator_factory();
	RefCountedPtr<mfconductivityopfactory>& factory = solver->get_opfact();

	// Go through each level and load balance. 
	for (int lvl = 0; lvl <= a_finest_level; lvl++){
	  // Make an operator
	  auto oper = RefCountedPtr<mfconductivityop> ((mfconductivityop*) factory->MGnewOp(domains[lvl], 0, false));

	  // Reset time
	  TimedDataIterator dit = a_grids[lvl].timedDataIterator();
	  dit.clearTime();
	  dit.enableTime();

	  // Do some relaxations on each level. This includes BCs.
	  for (int k = 0; k < num_apply; k++){
	    oper->applyOp(*dummy1[lvl], *dummy2[lvl], dit, true);
	  }

	  // Merge times
	  dit.disableTime();
	  dit.mergeTime();

	  // Now do the load balancing. When we do this, the boxes are standard-sorted! I.e. the new mesh ignores the desired sorting from amr_mesh
	  Vector<unsigned long long> loads = dit.getTime();
	  a_boxes[lvl] = dit.getBoxes();

	  // Do the desired sorting
	  load_balance::sort(a_boxes[lvl], loads, m_boxsort);

	  // Bah - need long long for LoadBalance
	  Vector<long long> lloads(loads.size());
	  for (int i = 0; i < loads.size(); i++){
	    lloads[i] = (long long) loads[i];
	  }
	  
	  // Do the friggin load balancing. 
	  load_balance::make_balance(a_procs[lvl], lloads, a_boxes[lvl]);
	}
      }
    }
  }
}
#include "CD_NamespaceFooter.H"
