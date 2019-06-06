/*!
  @file   mc_photo.H
  @brief  PIC solver class. 
  @author Robert Marskar
  @date   June 2019
*/

#ifndef _PIC_SOLVER_
#define _PIC_SOLVER_

#include "type_definitions.H"
#include "computational_geometry.H"
#include "amr_mesh.H"
#include "physical_domain.H"

#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>
#include <MeshInterp.H>

/*!
  @brief Particle-In-Cell (PIC) solver
*/
class pic_solver {
public:
  
  pic_solver();
  virtual ~pic_solver();

  void advance(const Real a_dt);
  void deposit(EBAMRCellData& a_density) const;
  void interpolate_force(const EBAMRCellData& a_E);
  void set_deposition_type();
  void set_bisect_step();
  void insert_particles(EBAMRParticles& a_particles);
  void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);
  void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);
  void set_amr(const RefCountedPtr<amr_mesh>& a_amr);
  void set_time(const int a_step, const Real a_time, const Real a_dt);

  virtual void allocate_internals();
  virtual void cache_state();
  virtual void deallocate_internals();
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);
  virtual void write_plot_file() const; 

  Real get_time() const;
  
  int query_ghost() const;
  int get_step() const;

  EBAMRParticles& get_particles();        // Get all particles
  EBAMRParticles& get_eb_particles();     // Get particles that collided with an EB
  EBAMRParticles& get_domain_particles(); // Get particles that collided with the domain
  EBAMRParticles& get_coll_particles();   // Get particles for collisions

protected:

  std::string m_name; // Solver name

  RefCountedPtr<computational_geometry> m_compgeom; // Computational domain
  RefCountedPtr<physical_domain> m_physdom;         // Physical domain
  RefCountedPtr<amr_mesh> m_amr;                    // AMR grids
  
  Real m_bisect_step;  // Bisection step length
  Real m_time;         // Simulation time. For staggered methods this is the positional time
  Real m_dt;           // Time step size

  int m_verbosity;  // Solver verbosity
  int m_step;       // Time step #

  InterpType m_deposition; // Mass deposition type

  EBAMRPVR m_pvr;             // Particle valid region
  EBAMRParticles m_particles; // Particles
  EBAMRParticles m_eb_part;   // Particles that struck the EB in a time step. Passed externally for boundary treatment
  EBAMRParticles m_dom_part;  // Particles that struck a domain in a time step. Passed externally for boundary treatment

  void advance_leapfrog(const Real a_dt); // Leapfrog advance algorithm
  void advance_vervlet(const Real a_dt);  // Velocity vervlet advance algorithm
};

#endif
