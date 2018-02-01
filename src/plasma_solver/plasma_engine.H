/*!
  @file plasma_engine.H
  @brief Declaration of the main class for time/space advancement of streamer equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PLASMA_ENGINE_
#define _PLASMA_ENGINE_

#include <RefCountedPtr.H>

#include "computational_geometry.H"
#include "plasma_kinetics.H"
#include "time_stepper.H"
#include "physical_domain.H"
#include "amr_mesh.H"
#include "cell_tagger.H"
#include "mfis.H"

/*!
  @brief Main class for time/space advancement of streamer equations
*/
class plasma_engine {
public:

  /*!
    @brief Full constructor. 
  */
  plasma_engine(const RefCountedPtr<physical_domain>&        a_physdom,
		const RefCountedPtr<computational_geometry>& a_compgeom,
		const RefCountedPtr<plasma_kinetics>&        a_plaskin,
		const RefCountedPtr<time_stepper>&           a_timestepper,
		const RefCountedPtr<amr_mesh>&               a_amr,
		const RefCountedPtr<cell_tagger>&            a_celltagger = RefCountedPtr<cell_tagger>(NULL));

  /*!
    @brief Constructor
  */
  ~plasma_engine();

  /*!
    @brief Setup function
  */
  virtual void setup_fresh();

  /*!
    @brief Set for restart
  */
  virtual void setup_for_restart(const std::string a_restart_file);

  /*!
    @brief Do initial regrids
  */
  virtual void initial_regrids(const int a_init_regrids);

  /*!
    @brief Verbosity
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set geometric refinement depth
  */
  virtual void set_geom_refinement_depth(const int a_depth = -1);

  /*!
    @brief Set geometrical refinement depths. This provides a little more flexibility than the basic call. 
    @param[in] a_depth1 Conductor tag depth
    @param[in] a_depth2 Dielectric tag depth
    @param[in] a_depth3 Conductor-gas interface tag depth
    @param[in] a_depth4 Dielectric-gas interface tag depth
    @param[in] a_depth5 Solid-gas interface tag depth
    @param[in] a_depth6 Solid-solid interface tag depth
  */
  virtual void set_geom_refinement_depth(const int a_depth1,
					 const int a_depth2,
					 const int a_depth3,
					 const int a_depth4,
					 const int a_depth5,
					 const int a_depth6);

  /*!
    @brief Set potential function
  */
  virtual void set_potential(Real (*a_potential)(const Real a_time, const RealVect a_pos));
					 

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();

protected:

  /*!
    @brief Constructor
  */
  plasma_engine();

  /*!
    @brief Index space
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief AMR
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Geometry
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief Plasma kinetics
  */
  RefCountedPtr<plasma_kinetics> m_plaskin;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<physical_domain> m_physdom;

  /*!
    @brief Time stpeper
  */
  RefCountedPtr<time_stepper> m_timestepper;

  /*!
    @brief Set the cell tagger
   */
  RefCountedPtr<cell_tagger> m_celltagger;

  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set the plasma kinetics
  */
  virtual void set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin);

  /*!
    @brief Set the time stepper
  */
  virtual void set_time_stepper(const RefCountedPtr<time_stepper>& a_timestepper);
  
    /*!
    @brief Set physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);



  /*!
    @brief Set the cell tagger
  */
  virtual void set_cell_tagger(const RefCountedPtr<cell_tagger>& a_celltagger);

protected:
  /*!
    @brief Write a plot file
  */
  virtual void write_plot_file();

  /*!
    @brief Write a checkpoint file
  */
  virtual void write_checkpoint_file();

  /*!
    @brief Read checkpoint file
  */
  virtual void read_checkpoint_file();

  /*!
    @brief Get geometric tags
  */
  virtual void get_geom_tags();

  /*!
    @brief Do a regrid step
  */
  virtual void regrid();

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
   */
  int m_time_step;

  /*!
    @brief Geometric tag depth
  */
  int m_geom_tag_depth;

  /*!
    @brief Geometric tag depth for dielectrics
  */
  int m_dielectric_tag_depth;

  /*!
    @brief Geometric tag depth for conductors
  */
  int m_conductor_tag_depth;

  /*!
    @brief Geometric tag depth for gas-solid interfaces
  */
  int m_gas_solid_interface_tag_depth;

  /*!
    @brief Geometry tag depth for solid-solid interfaces
  */
  int m_solid_solid_interface_tag_depth;

  /*!
    @brief Geometric tag depth for gas-dielectric interfaces
  */
  int m_gas_dielectric_interface_tag_depth;

  /*!
    @brief Geometric tag depth for gas-conductor interfaces
  */
  int m_gas_conductor_interface_tag_depth;

  /*!
    @brief Tags
  */
  Vector<IntVectSet> m_cell_tags;

  /*!
    @brief Tags
  */
  Vector<IntVectSet> m_geom_tags;

  /*!
    @brief Tag layouts
  */
  Vector<RefCountedPtr<LayoutData<IntVectSet> > > m_layout_tags;

  /*!
    @brief Potential
  */
  Real (*m_potential)(const Real a_time, const RealVect a_pos);

  /*!
    @brief Potential has been set
  */
  bool m_potential_set;
};
#endif
