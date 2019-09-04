/*!
  @file   cdr_fhd.H
  @brief  Flucuating hydrodynamics CDR solver based on implicit diffusion and Godunov
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _CDR_FHD_
#define _CDR_FHD_

#include "cdr_gdnv.H"

#include <random>
#include <time.h>

/*!
  @brief Godunov implementation for convection-diffusion-reaction equations. 

  @details Options that are available for this class
  
  Class options
  -------------

      cdr_fhd.divF_nc      = covered_face    # Valid options are "covered_face" and "conservative_average"
      cdr_fhd.limit_slopes = false           # Use van Leer slope limiter or not (default is false)

*/
class cdr_fhd : public cdr_gdnv {
public:

  /*!
    @brief Constructor
  */
  cdr_fhd();

  /*!
    @brief Destructor
  */
  virtual ~cdr_fhd();

  /*!
    @brief Parse class options
  */
  virtual void parse_options();

  /*!
    @brief Pure Euler advance without source terms
  */
  virtual void advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt);

  /*!
    @brief Pure Euler advance with source terms
  */
  virtual void advance_euler(EBAMRCellData&       a_new_state,
			     const EBAMRCellData& a_old_state,
			     const EBAMRCellData& a_source,
			     const Real           a_dt);

  /*!
    @brief Pure TGA advance without source terms
  */
  virtual void advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt);

  /*!
    @brief Pure Euler advance with source terms
  */
  virtual void advance_tga(EBAMRCellData&       a_new_state,
			   const EBAMRCellData& a_old_state,
			   const EBAMRCellData& a_source,
			   const Real           a_dt);

  virtual void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist);

  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);
protected:

  bool m_stochastic_diffusion;
  bool m_stochastic_advection;
  bool m_plot_numbers;

  void GWN_advection_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states);
  void GWN_diffusion_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states);
  void smooth_heaviside_faces(EBAMRFluxData& a_face_states, const EBAMRCellData& a_cell_states);
  void fill_GWN(EBAMRFluxData& a_noise, const Real a_mean);

  std::mt19937_64* m_rng;     // Engine
  int  m_seed;                // RNG seed

  // Various routines for setting class options. 
  virtual void parse_diffusion();    // Parses stochastic diffusion
  virtual void parse_rng_seed();     // Parses RNG seed
  virtual void parse_plotmode();     // Output with numbers rather 
};

#endif
