/*!
  @file plasma_engine.H
  @brief Declaration of the main class for time/space advancement of streamer equations
  @author Robert Marskar
  @date Nov. 2017
  @todo Move the BC stuff to poisson_solver
*/

#ifndef _PLASMA_ENGINE_
#define _PLASMA_ENGINE_

#include <RefCountedPtr.H>

#include "computational_geometry.H"
#include "plasma_kinetics.H"
#include "time_stepper.H"
#include "wall_bc.H"
#include "physical_domain.H"
#include "amr_mesh.H"
#include "cell_tagger.H"
#include "mfis.H"

namespace Potential {
  enum GroundLive {
    Ground = false,
    Live   = true
  };
};

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
    @brief Set Dirichlet wall BC
  */
  virtual void set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const Potential::GroundLive a_live);

  /*!
    @brief Set Neumann wall BC
  */
  virtual void set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value);

  /*!
    @brief Get Poisson BC at specified domain wall
  */
  virtual wall_bc& get_wall_bc(const int a_dir, Side::LoHiSide a_side) const;

  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set geometric refinement depth
  */
  virtual void set_geom_refinement_depth(const int a_depth);

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
    @brief Wall boundary conditions for Poisson equation. Must be set by the plasma engine. 
  */
  Vector<RefCountedPtr<wall_bc> > m_wallbc;

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



  /*!
    @brief Allocate data holders for wall bcs so they can be set
  */
  virtual void allocate_wall_bc();

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
    @brief Internal bc map for Poisson stuff
  */
  virtual int map_bc(const int a_dir, const Side::LoHiSide a_side) const;

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
};
#endif
