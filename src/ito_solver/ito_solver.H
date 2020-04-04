/*!
  @file   ito_solver.H
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _ITO_SOLVER_
#define _ITO_SOLVER_

#include "amr_mesh.H"
#include "type_definitions.H"
#include "computational_geometry.H"

#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>
#include <MeshInterp.H>

#include <random>


/*!
  @brief Base class for Ito diffusion particle models
*/
class ito_solver {
public:

  /*!
    @brief Default ctor
  */
  ito_solver();

  /*!
    @Default dtor
  */
  virtual ~ito_solver();

  /*!
    @brief Get solver name
  */
  virtual std::string get_name();

  /*!
    @brief Parse options
  */
  virtual void parse_options();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> get_plotvar_names() const;

  /*!
    @brief Set computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

    /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set phase
  */
  virtual void set_phase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Write plot data
  */
  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);  

  
protected:

    /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Name
  */
  std::string m_name;

  /*!
    @brief Class name
  */
  std::string m_class_name;

  /*!
    @brief RNG engine
  */
  std::mt19937_64 m_rng; 

  /*!
    @brief Uniform real distribution on [0,1]
  */
  std::uniform_real_distribution<Real> m_udist01; 

  /*!
    @brief Uniform real distribution on [-1,1]
  */
  std::uniform_real_distribution<Real> m_udist11;

  /*!
    @brief Gaussian distribution centered on 0, std=1
  */
  std::normal_distribution<Real> m_gauss01;

  /*!
    @brief Initial seed for RNG
  */
  Real m_seed_rng;

  /*!
    @brief Bisection step for particle trajectories
  */
  Real m_bisect_step;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_step;

  /*!
    @brief Current time
  */
  Real m_time;

  /*!
    @brief Time step
  */
  Real m_dt;

  /*!
    @brief Solver is diffusive
  */
  bool m_diffusive;

  /*!
    @brief Also solve for advection/convection
  */
  bool m_mobile;

  /*!
    @brief Output state
  */
  bool m_plot_phi;

  /*!
    @brief Mass/charge deposition type
  */
  InterpType m_deposition; 

  /*!
    @brief Plot deposition type
  */
  InterpType m_plot_deposition;

  /*!
    @brief Deposited particles
  */
  EBAMRCellData m_state;

  /*!
    @brief Particles
  */
  EBAMRParticles m_particles;     

  /*!
    @brief Particle cache
  */
  EBAMRParticles m_particleCache; 

  /*!
    @brief Particle valid region
  */
  EBAMRPVR m_pvr;               
  
  /*!
    @brief Parse RNG options
  */
  void parse_rng();

  /*!
    @brief Parse deposition options
  */
  void parse_deposition();

  /*!
    @brief Parse bisection step
  */
  void parse_bisect_step();

  /*!
    @brief Parse PVR buffer
  */
  void parse_pvr_buffer();

  /*!
    @brief Parse plot variables
  */
  void parse_plot_vars();

};
#endif
