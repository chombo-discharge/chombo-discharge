/*!
  @file cdr_solver.H
  @brief Declaration of an abstract class for evolving cdr equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_SOLVER_
#define _CDR_SOLVER_

#include "type_definitions.H"
#include "computational_geometry.H"

/*!
  @brief Abstract class for cdr-equations. This is (almost) a whole solver. 
*/
class cdr_solver {
public:
  
  /*!
    @brief Full constructor. Subsequently calls define. 
  */
  cdr_solver(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Constructor
  */
  virtual ~cdr_solver();

  /*!
    @brief Advance method
  */
  virtual void advance(const Real& a_dt);

  /*!
    @brief Advance state
  */
  virtual void advance(EBAMRCellData& a_state, const Real& a_dt) = 0;

  /*!
    @brief Compute right-hand side
  */
  virtual void compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt); 

  /*!
    @brief Compute the right-hand side using supplied solution
  */
  virtual void compute_rhs(EBAMRCellData& a_rhs, const EBAMRCellData& a_state, const Real& a_dt);

  /*!
    @brief Get the internal state
   */
  EBAMRCellData& get_state();

    /*!
    @brief Get the source term
   */
  EBAMRCellData& get_source();

  /*!
    @brief Get the face-centered velocities
   */
  EBAMRFluxData& get_velo_face();

  /*!
    @brief Get the eb-centered velocities
  */
  EBAMRIVData& get_velo_eb();

  /*!
    @brief Get the face-centered diffusion coefficient
  */
  EBAMRFluxData& get_diffco_face();

  /*!
    @brief Get the EB-centered diffusion coefficient
  */
  EBAMRIVData& get_diffco_eb();  


protected:

  /*!
    @brief Default constructor. Stay away.
  */
  cdr_solver();

  /*!
    @brief Set computational geometry
  */
  void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Allocate temporaries
  */
  virtual void allocate_temporaries();

  /*!
    @brief Delete temporaries
  */
  virtual void delete_temporaries();

  /*!
    @brief Extrapolate state to faces. 
  */
  virtual void extrapolate_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state) = 0;

  /*!
    @brief Compute the conservative divergence
  */
  virtual void conservative_divergence(EBAMRCellData&       a_cons_div,
				       const EBAMRFluxData& a_face_vel,
				       const EBAMRFluxData& a_face_state);
  /*!
    @brief Compute the non-conservative divergence
  */
  virtual void nonconservative_divergence(EBAMRCellData& a_noncons_div, const EBAMRCellData& a_cons_div);

  /*!
    @brief Compute the hybrid divergence
  */
  virtual void hybrid_divergence(EBAMRCellData&       a_hybrid_div,
				 EBAMRIVData&         a_mass_diff,
				 const EBAMRCellData& a_cons_div,
				 const EBAMRCellData& a_noncons_div);

  /*!
    @brief Increment flux registers. 
  */
  virtual void increment_flux_registers(const EBAMRCellData& a_face_state, const EBAMRCellData& a_velo_face);

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  virtual void coarse_fine_increment(const EBAMRIVData& m_mass_diff);

  /*!
    @brief Do mass redistribution
  */
  virtual void hyperbolic_redistribution(EBAMRCellData&       a_del_vel_rho,
					 const EBAMRIVData&   a_mass_diff,
					 const EBAMRCellData& a_redist_weights);

  /*!
    @brief Do the reflux-redist dance (only for EB x CF)
  */
  virtual void reflux_redist_dance(EBAMRCellData& a_del_vel_rho, const EBAMRCellData& a_redist_weights);

  /*!
    @brief Set weights for reflux-redist dance
  */
  virtual void set_weights(const EBAMRCellData& a_redist_weights);

  /*!
    @brief Increment redistribution flux registers
  */
  virtual void increment_redist_flux();

  /*!
    @brief Do coarse-fine redistribution
  */
  virtual void coarse_fine_redistribution(EBAMRCellData& a_state);

    /*!
    @brief Reflux state
  */
  virtual void reflux(EBAMRCellData& a_state);
  
  /*!
    @brief Computational geometry
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief State data (i.e. the ion density)
  */
  EBAMRCellData m_state;

  /*!
    @brief Source term 
  */
  EBAMRCellData m_source;

  /*!
    @brief Face-centered velocities (only normal components)
  */
  EBAMRFluxData m_velo_face;

  /*!
    @brief EB-centered velocities
  */
  EBAMRIVData m_velo_eb;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRFluxData m_diffco_face;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRIVData m_diffco_eb;
};

#endif
