/*!
  @file   ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef _ITO_PLASMA_STEPPER_
#define _ITO_PLASMA_STEPPER_

#include "time_stepper.H"
#include "ito_plasma_physics.H"
#include "ito_layout.H"
#include "rte_layout.H"
#include "mc_photo.H"
#include "poisson_solver.H"
#include "sigma_solver.H"

namespace physics {
  namespace ito_plasma {

    class ito_plasma_stepper : public time_stepper {
    public:

      static Real s_constant_one(const RealVect a_pos);

      // Ctor and dtor
      ito_plasma_stepper();
      ito_plasma_stepper(RefCountedPtr<ito_plasma_physics>& a_physics);
      ~ito_plasma_stepper();

      // Setup routines
      virtual void parse_options() = 0;
      virtual void allocate_internals() = 0;
      virtual void setup_solvers() override;
      virtual void allocate() override;
      virtual void initial_data() override;
      virtual void initial_sigma();
      virtual void post_checkpoint_setup() override;
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
      int  get_num_plot_vars() const override;
  
      // Advance routines
      virtual Real advance(const Real a_dt) = 0;
      virtual void compute_dt(Real& a_dt, time_code& a_timecode) override;
      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override;

      // Registration routines
      void register_realms() override;
      void register_operators() override;
      
      // New regrid routines
      void deallocate() override;
      virtual void pre_regrid(const int a_lmin, const int a_old_finest_level) override;
      virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;
      void post_regrid() override;
      bool load_balance(Vector<Vector<int> >&            a_procs,
			Vector<Vector<Box> >&            a_boxes,
			const std::string                a_realm,
			const Vector<DisjointBoxLayout>& a_grids,
			const int                        a_lmin,
			const int                        a_finest_level) override;


      // Internal routines
      void set_verbosity(const int a_verbosity);

      // Set routines
      void set_ito(RefCountedPtr<ito_layout<ito_solver> >& a_cdr);
      void set_poisson(RefCountedPtr<poisson_solver>& a_poisson);
      void set_rte(RefCountedPtr<rte_layout<mc_photo> >& a_rte);
      void set_potential(Real (*a_potential)(const Real a_time));

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

      // Compute J routines
      void compute_J(EBAMRCellData& a_J, const Real a_dt);
      void compute_J(LevelData<EBCellFAB>& a_J, const int a_level, const Real a_dt);
      void compute_J(EBCellFAB& a_J, const int a_level, const DataIndex a_dit, const Box& a_box, const Real a_dt);

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

      // Deposition routines
      void deposit_particles();
      void deposit_mobile_particles();
      void deposit_diffusive_particles();
      void deposit_mobile_or_diffusive_particles();
      void deposit_stationary_particles();

      // Remap routines
      void remap_particles();
      void remap_mobile_particles();
      void remap_diffusive_particles();
      void remap_mobile_or_diffusive_particles();

      // This sets the velocity function for each species. This is sgn(charge)*E
      void set_ito_velocity_funcs();

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

      // Routines for computing velocities
      void compute_ito_velocities_lfa();
      void compute_ito_velocities_lea();


      // Routines for computing velocities.
      // This set of routines computes v = v(E) and then interpolates that to the particle positions. This
      void compute_ito_velocities();
      void compute_ito_velocities(Vector<EBAMRCellData*>&       a_velo_funcs,
				  const Vector<EBAMRCellData*>& a_densities,
				  const EBAMRCellData&          a_E,
				  const Real                    a_time);
      void compute_ito_velocities(Vector<LevelData<EBCellFAB>* >&       a_velo_funcs,
				  const Vector<LevelData<EBCellFAB>* >& a_densities,
				  const LevelData<EBCellFAB>&           a_E,
				  const int                             a_level,
				  const Real                            a_time);
      void compute_ito_velocities(Vector<EBCellFAB*>&       a_velo,
				  const Vector<EBCellFAB*>& a_densities,
				  const EBCellFAB&          a_E,
				  const int                 a_level,
				  const DataIndex           a_dit,
				  const Box                 a_box,
				  const Real                a_time);

      // Routines for computing diffusion coefficients.
      // These routines compute the diffusion coefficients by computing D(E)
      // as a mesh variable, and then interpolating that to the particle position. 
      void compute_ito_diffusion();
      void compute_ito_diffusion(Vector<EBAMRCellData*>&       a_diffco,
				 const Vector<EBAMRCellData*>& a_densities,
				 const EBAMRCellData&          a_E,
				 const Real                    a_time);
      void compute_ito_diffusion(Vector<LevelData<EBCellFAB>* >&       a_diffco,
				 const Vector<LevelData<EBCellFAB>* >& a_densities,
				 const LevelData<EBCellFAB>&           a_E,
				 const int                             a_level,
				 const Real                            a_time);
      void compute_ito_diffusion(Vector<EBCellFAB*>&       a_diffco,
				 const Vector<EBCellFAB*>& a_densities,
				 const EBCellFAB&          a_E,
				 const int                 a_level,
				 const DataIndex           a_dit,
				 const Box                 a_box,
				 const Real                a_time);

      // Routines for chemistry advances
      void advance_reaction_network(const Real a_dt);
      void advance_reaction_network(Vector<particle_container<ito_particle>* >& a_particles,
				    Vector<particle_container<photon>* >&       a_photons,
				    Vector<particle_container<photon>* >&       a_newPhotons,
				    const EBAMRCellData&                        a_E,
				    const Real                                  a_dt);
      void advance_reaction_network(Vector<AMRCellParticles<ito_particle>* >& a_particles,
				    Vector<AMRCellParticles<photon>* >&       a_photons,
				    Vector<AMRCellParticles<photon>* >&       a_newPhotons,
				    const EBAMRCellData&                      a_E,
				    const Real                                a_dt);
      void advance_reaction_network(Vector<LayoutData<BinFab<ito_particle> >* >& a_particles,
				    Vector<LayoutData<BinFab<photon> >* >&       a_photons,
				    Vector<LayoutData<BinFab<photon> >* >&       a_newPhotons,
				    const LevelData<EBCellFAB>&                  a_E,
				    const int                                    a_lvl,
				    const Real                                   a_dt);

      void advance_reaction_network(Vector<BinFab<ito_particle>* >& a_particles,
				    Vector<BinFab<photon>* >&       a_photons,
				    Vector<BinFab<photon>* >&       a_newPhotons,
				    const EBCellFAB&                a_E,
				    const int                       a_lvl,
				    const DataIndex                 a_dit,
				    const Box                       a_box,
				    const Real                      a_dx,
				    const Real                      a_dt);

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

    protected:

      std::string m_fluid_realm;     // Realm for Poisson and sigma solver. Defaults to realm::primal
      std::string m_particle_realm;  // Realm for particles. Defaults to realm::primal
      std::string m_name;
      phase::which_phase m_phase;
      

      RefCountedPtr<ito_plasma_physics> m_physics;
      RefCountedPtr<ito_layout<ito_solver> > m_ito;
      RefCountedPtr<rte_layout<mc_photo> > m_rte;
      RefCountedPtr<poisson_solver> m_poisson;
      RefCountedPtr<sigma_solver> m_sigma;

      Real (*m_potential)(const Real a_time);

      bool m_regrid_superparticles;
      bool m_load_balance;
      
      int m_ppc;
      Real m_max_cells_hop;
      Real m_dt_relax;
      Real m_min_dt;
      Real m_max_dt;

      // Storage for current density. This is defined over the particle realm. 
      EBAMRCellData m_J;
      EBAMRCellData m_fluid_scratch1;     // Scratch storage on the fluid realm with 1 component
      EBAMRCellData m_fluid_scratchD;     // Scratch storage on the fluid realm with SpaceDim components

      EBAMRCellData m_particle_scratch1;  // Scratch storage on the particle realm with 1 component
      EBAMRCellData m_particle_scratchD;  // Scratch storage on the particle realm with SpaceDim components

      // Storage for the electric field on 
      EBAMRCellData m_fluid_E;
      EBAMRCellData m_particle_E;


      bool load_balance_particle_realm(Vector<Vector<int> >&            a_procs,
				       Vector<Vector<Box> >&            a_boxes,
				       const std::string                a_realm,
				       const Vector<DisjointBoxLayout>& a_grids,
				       const int                        a_lmin,
				       const int                        a_finest_level);
    };
  }
};
#endif
