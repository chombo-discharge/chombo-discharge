/*!
  @file   ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef _ITO_PLASMA_STEPPER_
#define _ITO_PLASMA_STEPPER_

#include <functional>

#include "time_stepper.H"
#include "ito_plasma_physics.H"
#include "ito_layout.H"
#include "point_particle.H"
#include "rte_layout.H"
#include "mc_photo.H"
#include "CD_FieldSolver.H"
#include "sigma_solver.H"

#include "CD_NamespaceHeader.H"

namespace physics {
  namespace ito_plasma {

    // Differentiation between particles. 
    enum class which_particles {
      all,
      all_mobile,
      all_diffusive,
      charged_mobile,
      charged_diffusive,
      all_mobile_or_diffusive,
      charged_and_mobile_or_diffusive,
      stationary,
    };

    class ito_plasma_stepper : public time_stepper {
    public:

      static Real s_constant_one(const RealVect a_pos);

      // Ctor and dtor
      ito_plasma_stepper();
      ito_plasma_stepper(RefCountedPtr<ito_plasma_physics>& a_physics);
      ~ito_plasma_stepper();

      // Setup routines
      virtual void parseOptions() = 0;
      virtual void parseRuntimeOptions() = 0;
      virtual void allocate_internals() = 0;
      virtual void setup_solvers() override;
      virtual void allocate() override;
      virtual void initial_data() override;
      virtual void initial_sigma();
      virtual void post_checkpoint_setup() override;
      virtual void post_checkpoint_poisson();
      virtual void post_initialize() override;
      
      void setup_ito();
      void setup_poisson();
      void setup_rte();
      void setup_sigma();
  
      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const override;
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      virtual void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;
      void write_J(EBAMRCellData& a_output, int& a_icomp) const;
      void write_num_particles_per_patch(EBAMRCellData& a_output, int& a_icomp) const;
      int  get_num_plot_vars() const override;
      Vector<long int> get_checkpoint_loads(const std::string a_realm, const int a_level) const override;
  
      // Advance routines
      virtual Real advance(const Real a_dt) = 0;
      virtual void compute_dt(Real& a_dt, time_code& a_timecode) override;
      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override;

      // Registration routines
      void registerRealms() override;
      void registerOperators() override;
      
      // Regrid routines. 
      void deallocate() override;
      virtual void pre_regrid(const int a_lmin, const int a_old_finest_level) override;
      virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;
      void post_regrid() override;
      bool load_balance_realm(const std::string a_realm) const override;
      void load_balance_boxes(Vector<Vector<int> >&            a_procs,
			      Vector<Vector<Box> >&            a_boxes,
			      const std::string                a_realm,
			      const Vector<DisjointBoxLayout>& a_grids,
			      const int                        a_lmin,
			      const int                        a_finest_level) override;


      // Particle handling routines.
      void intersect_particles(const which_particles a_which, const EB_representation a_representation, const bool a_delete);
      void intersect_particles(const which_particles a_which,
			       const ito_solver::which_container a_particles,
			       const ito_solver::which_container a_eb_particles,
			       const ito_solver::which_container a_domain_particles,
			       const EB_representation           a_representation,
			       const bool                        a_delete);
			       
      void remove_covered_particles(const which_particles a_which, const EB_representation a_representation, const Real a_tolerance);
      void remove_covered_particles(const which_particles             a_which,
				    const ito_solver::which_container a_container,
				    const EB_representation           a_representation,
				    const Real                        a_tolerance);

      void transfer_covered_particles(const which_particles a_which, const EB_representation a_representation, const Real a_tolerance);
      void transfer_covered_particles(const which_particles             a_which,
				      const ito_solver::which_container a_containerFrom,
				      const ito_solver::which_container a_containerTo,
				      const EB_representation           a_representation,
				      const Real                        a_tolerance);

      // Internal routines
      void set_verbosity(const int a_verbosity);

      // Set routines
      void set_ito(RefCountedPtr<ito_layout<ito_solver> >& a_cdr);
      void set_poisson(RefCountedPtr<FieldSolver>& a_poisson);
      void set_rte(RefCountedPtr<rte_layout<mc_photo> >& a_rte);
      void set_potential(const std::function<Real(const Real a_time)>& a_potential);

      // Get routines
      Real get_time() const;
      Real compute_Emax(const phase::which_phase a_phase);

      // E-field computations
      void compute_E(MFAMRCellData& a_E,      const MFAMRCellData& a_potential);
      void compute_E(EBAMRCellData& a_E,      const phase::which_phase a_phase);
      void compute_E(EBAMRCellData& a_E,      const phase::which_phase a_phase, const MFAMRCellData& a_potential);
      void compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);
      void compute_E(EBAMRIVData& a_E_eb,     const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);
      
      // Compute rho routines
      void compute_rho();
      void compute_rho(EBAMRCellData& a_rho, const phase::which_phase a_phase);
      void compute_rho(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>&  a_densities);

      // Compute conductivity routines
      void compute_conductivity(EBAMRCellData& a_conductivity);
      void compute_conductivity(EBAMRCellData& a_conductivity, const Vector<particle_container<ito_particle>* >& a_particles);

      // Compute J routines
      void compute_J(EBAMRCellData& a_J, const Real a_dt);

      // Relaxation time
      Real compute_relaxation_time();
      Real compute_relaxation_time(const int a_level);
      Real compute_relaxation_time(const int a_level, const DataIndex a_dit);

      // Poisson-solving routines
      bool solve_poisson();
      bool solve_poisson(MFAMRCellData&                 a_potential,
			 MFAMRCellData&                 a_rhs,
			 const Vector<EBAMRCellData* >& a_densities,
			 const EBAMRIVData&             a_sigma);

      // Remap and deposition routines for solver subsets.
      void deposit_particles(const which_particles a_which_particles);
      void deposit_particles(const which_particles a_which_particles, const ito_solver::which_container a_container);
      void remap_particles(const which_particles a_which_particles);
      void remap_particles(const which_particles a_which_particles, const ito_solver::which_container a_container);

      // This sets the velocity function for each species. This is sgn(charge)*E
      void set_ito_velocity_funcs();

      // Routines for computing velocities and diffusion coefficients
      void compute_ito_velocities();
      void compute_ito_diffusion();

      // Routines for computing mobilities in the LFA. This computes mesh-defined mobilities and interpolates them to the particle positions
      void compute_ito_mobilities_lfa();
      void compute_ito_mobilities_lfa(Vector<EBAMRCellData*>& a_meshMobilities, const EBAMRCellData& a_E, const Real a_time);
      void compute_ito_mobilities_lfa(Vector<LevelData<EBCellFAB>* >& a_meshMobilities,
				      const LevelData<EBCellFAB>&     a_E,
				      const int                       a_level,
				      const Real                      a_time);
      void compute_ito_mobilities_lfa(Vector<EBCellFAB*>& a_meshMobilities,
				      const EBCellFAB&    a_E,
				      const int           a_level,
				      const DataIndex     a_dit,
				      const Box           a_box,
				      const Real          a_time);

      // Routines for computing mobilities in the LEA. This runs through all the particles and calls a mobility computation function,
      // which you must implement in ito_species
      void compute_ito_mobilities_lea();
      void compute_ito_mobilities_lea(Vector<particle_container<ito_particle>* >& a_particles);
      void compute_ito_mobilities_lea(Vector<ParticleData<ito_particle>* >& a_particles);

      // Routines for computing diffusion coefficients.
      // These routines compute the diffusion coefficients by computing D(E)
      // as a mesh variable, and then interpolating that to the particle position. 
      void compute_ito_diffusion_lfa();
      void compute_ito_diffusion_lfa(Vector<EBAMRCellData*>&       a_diffco,
				     const Vector<EBAMRCellData*>& a_densities,
				     const EBAMRCellData&          a_E,
				     const Real                    a_time);
      void compute_ito_diffusion_lfa(Vector<LevelData<EBCellFAB>* >&       a_diffco,
				     const Vector<LevelData<EBCellFAB>* >& a_densities,
				     const LevelData<EBCellFAB>&           a_E,
				     const int                             a_level,
				     const Real                            a_time);
      void compute_ito_diffusion_lfa(Vector<EBCellFAB*>&       a_diffco,
				     const Vector<EBCellFAB*>& a_densities,
				     const EBCellFAB&          a_E,
				     const int                 a_level,
				     const DataIndex           a_dit,
				     const Box                 a_box,
				     const Real                a_time);

      // Routines for computing diffusion coefficients in the LEA
      void compute_ito_diffusion_lea();

      // Routines for computing the number of particles per cell. This is done on the particle realm. 
      void compute_reactive_particles_per_cell(EBAMRCellData& a_ppc);
      void compute_reactive_particles_per_cell(LevelData<EBCellFAB>& a_ppc, const int a_level);
      void compute_reactive_particles_per_cell(EBCellFAB& a_ppc, const int a_level, const DataIndex a_dit, const Box a_box, const EBISBox& a_ebisbox);

      // Routines for computing the mean energy of the particles in grid cells
      void compute_reactive_mean_energies_per_cell(EBAMRCellData& a_mean_energies);
      void compute_reactive_mean_energies_per_cell(LevelData<EBCellFAB>& a_mean_energies, const int a_level);
      void compute_reactive_mean_energies_per_cell(EBCellFAB& a_mean_energies, const int a_level, const DataIndex a_dit, const Box a_box, const EBISBox& a_ebisbox);

      // New routines for chemistry advances. 
      void advance_reaction_network_nwo(const Real a_dt);
      void advance_reaction_network_nwo(const EBAMRCellData& a_E, const EBAMRCellData& a_EdotJ, const Real a_dt);
      void advance_reaction_network_nwo(LevelData<EBCellFAB>&       a_ppc,
					LevelData<EBCellFAB>&       a_ypc,
					LevelData<EBCellFAB>&       a_mean_energies,
					const LevelData<EBCellFAB>& a_E,
					const LevelData<EBCellFAB>& a_EdotJ,
					const int                   a_level,
					const Real                  a_dt);
      void advance_reaction_network_nwo(EBCellFAB&       a_particlesPerCell,
					EBCellFAB&       a_newPhotonsPerCell,
					EBCellFAB&       a_meanParticleEnergies,
					const EBCellFAB& a_E,
					const EBCellFAB& a_EdotJ,
					const int        a_level,
					const DataIndex  a_dit,
					const Box        a_box,
					const Real       a_dx,
					const Real       a_dt);

      // Routines for reconciling particles. 
      void reconcile_particles(const EBAMRCellData& a_newParticlesPerCell,
			       const EBAMRCellData& a_oldParticlesPerCell,
			       const EBAMRCellData& a_meanParticleEnergies,
			       const EBAMRCellData& a_newPhotonsPerCell);
      void reconcile_particles(const LevelData<EBCellFAB>& a_newParticlesPerCell,
			       const LevelData<EBCellFAB>& a_oldParticlesPerCell,
			       const LevelData<EBCellFAB>& a_meanParticleEnergies,
			       const LevelData<EBCellFAB>& a_newPhotonsPerCell,
			       const int                   a_level);
      void reconcile_particles(const EBCellFAB& a_newParticlesPerCell,
			       const EBCellFAB& a_oldParticlesPerCell,
			       const EBCellFAB& a_meanParticleEnergies,
			       const EBCellFAB& a_newPhotonsPerCell,
			       const int        a_level,
			       const DataIndex  a_dit,
			       const Box        a_box,
			       const Real       a_dx);


      // Old routines for chemistry advances. The physics module will switch between LFA and LEA.
      void advance_reaction_network(const Real a_dt);
      void advance_reaction_network(Vector<particle_container<ito_particle>* >& a_particles,
				    Vector<particle_container<photon>* >&       a_photons,
				    Vector<particle_container<photon>* >&       a_newPhotons,
				    const Vector<EBAMRCellData>&                a_sources,
				    const EBAMRCellData&                        a_E,
				    const Real                                  a_dt);
      void advance_reaction_network(Vector<AMRCellParticles<ito_particle>* >& a_particles,
				    Vector<AMRCellParticles<photon>* >&       a_photons,
				    Vector<AMRCellParticles<photon>* >&       a_newPhotons,
				    const Vector<EBAMRCellData>&              a_sources,
				    const EBAMRCellData&                      a_E,
				    const Real                                a_dt);
      void advance_reaction_network(Vector<LayoutData<BinFab<ito_particle> >* >& a_particles,
				    Vector<LayoutData<BinFab<photon> >* >&       a_photons,
				    Vector<LayoutData<BinFab<photon> >* >&       a_newPhotons,
				    const Vector<LevelData<EBCellFAB>* >&        a_sources,
				    const LevelData<EBCellFAB>&                  a_E,
				    const int                                    a_lvl,
				    const Real                                   a_dt);
      void advance_reaction_network(Vector<BinFab<ito_particle>* >& a_particles,
				    Vector<BinFab<photon>* >&       a_photons,
				    Vector<BinFab<photon>* >&       a_newPhotons,
				    const Vector<EBCellFAB*>&       a_sources,
				    const EBCellFAB&                a_E,
				    const int                       a_lvl,
				    const DataIndex                 a_dit,
				    const Box                       a_box,
				    const Real                      a_dx,
				    const Real                      a_dt);

      // Compute dt routines
      Real compute_physics_dt() const;
      Real compute_physics_dt(const EBAMRCellData& a_E, const Vector<EBAMRCellData*> a_densities) const;
      Real compute_physics_dt(const LevelData<EBCellFAB>& a_E, const Vector<LevelData<EBCellFAB> *> a_densities, const int a_level) const;
      Real compute_physics_dt(const EBCellFAB& a_E, const Vector<EBCellFAB*> a_densities, const int a_level, const DataIndex a_dit, const Box a_box) const;

      // Routines for advancing photons
      void advance_photons(const Real a_dt);

      // Routines for sorting particles. These routines are here because rte_layout is general for rte_solver, and I'm too
      // lazy to write an mc_photo_layout. So here you go. 
      void sort_photons_by_cell();
      void sort_photons_by_patch();
      void sort_source_photons_by_cell();
      void sort_source_photons_by_patch();
      void sort_bulk_photons_by_cell();
      void sort_bulk_photons_by_patch();
      void sort_eb_photons_by_cell();
      void sort_eb_photons_by_patch();
      void sort_domain_photons_by_cell();
      void sort_domain_photons_by_patch();

      // Routines for computing EdotJ per species, including the diffusive contribution. a_source is defined on the particle realm.
      void compute_EdotJ_source(const Real a_dt);
      void compute_EdotJ_source();
      void compute_EdotJ_source_nwo();
      void compute_EdotJ_source_nwo2(const Real a_dt);


    protected:

      std::string m_fluid_realm;     // Realm for Poisson and sigma solver. Defaults to realm::primal
      std::string m_particle_realm;  // Realm for particles. Defaults to realm::primal
      std::string m_name;
      phase::which_phase m_phase;

      BoxSorting m_boxSort;
      
      RefCountedPtr<ito_plasma_physics> m_physics;
      RefCountedPtr<ito_layout<ito_solver> > m_ito;
      RefCountedPtr<rte_layout<mc_photo> > m_rte;
      RefCountedPtr<FieldSolver> m_fieldSolver;
      RefCountedPtr<sigma_solver> m_sigma;

      // Filters that can be applied. 
      std::deque<std::tuple<Real, int, int> > m_filters;

      std::function<Real(const Real a_time)> m_potential;

      bool m_regrid_superparticles;
      bool m_load_balance;
      bool m_nwo_reactions; // For switching between the old and new reaction algorithms.
      bool m_filter_rho;
      
      int m_ppc;
      int m_halo_buffer;
      int m_pvr_buffer;
      int m_load_balance_idx;
      
      Real m_load_ppc;
      Real m_max_cells_hop;
      Real m_dt_relax;
      Real m_min_dt;
      Real m_max_dt;

      // Storage for current density. This is defined over the particle realm. 
      EBAMRCellData m_J;
      mutable EBAMRCellData m_fluid_scratch1;     // Scratch storage on the fluid realm with 1 component
      mutable EBAMRCellData m_fluid_scratchD;     // Scratch storage on the fluid realm with SpaceDim components

      mutable EBAMRCellData m_particle_scratch1;  // Scratch storage on the particle realm with 1 component
      mutable EBAMRCellData m_particle_scratchD;  // Scratch storage on the particle realm with SpaceDim components

      mutable Vector<EBAMRCellData> m_fscratch1; // Scratch storage with a single component for all species. On the fluid realm. 
      mutable Vector<EBAMRCellData> m_fscratch2; // Scratch storage with a single component for all species. On the fluid realm.

      EBAMRCellData m_particle_ppc; // For holding the number of particles per cell for all species. Contains compoents = number of species
      EBAMRCellData m_particle_ypc; // For holding the number of generated photons per cell on the particle realm.
      EBAMRCellData m_particle_eps; // For holding the mean energy on the particle realm
      EBAMRCellData m_particle_old; // For holding the mean energy on the particle realm
      
      EBAMRCellData m_fluid_ppc;    // For holding the number of particles per cell for all species. Contains compoents = number of species
      EBAMRCellData m_fluid_ypc;    // For holding the number of generated photons per cell on the fluid realm.
      EBAMRCellData m_fluid_eps;    // For holding the mean energy on the particle realm

      // Storage for the electric field on 
      EBAMRCellData m_fluid_E;
      EBAMRCellData m_particle_E;

      // Storage for energy source terms
      Vector<EBAMRCellData> m_energy_sources;
      EBAMRCellData m_EdotJ;

      // For parsing filters
      void parse_filters();

      // Set the particle buffers
      void set_particle_buffers();

      // Diagnostics tool
      void get_particle_statistics(Real& a_avg, Real& a_sigma, size_t& a_minParticles, size_t& a_maxParticles, int& a_minRank, int& a_maxRank);
      void print_timer_head();
      void print_timer_diagnostics(Real& a_timer, const std::string a_prefix); // Timer should be in percentage 
      void print_timer_tail();


      // Do load balancing of the particle realm. 
      void load_balance_particle_realm(Vector<Vector<int> >&            a_procs,
				       Vector<Vector<Box> >&            a_boxes,
				       const std::string                a_realm,
				       const Vector<DisjointBoxLayout>& a_grids,
				       const int                        a_lmin,
				       const int                        a_finest_level);
      
      Vector<RefCountedPtr<ito_solver> > get_lb_solvers() const;
    };
  }
};
#include "CD_NamespaceFooter.H"
#endif
