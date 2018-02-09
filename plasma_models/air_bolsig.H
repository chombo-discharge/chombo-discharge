/*!
  @file   air_bolsig.H
  @brief  Declaration of a BOLSIG-based transport model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _AIR_BOLSIG_
#define _AIR_BOLSIG_

/*!
  @brief BOLSIG-based kinetic model for air. Transport coefficients for electrons are taken by calls to BOLSIG+
*/
class air_bolsig : public plasma_kinetics {
public:

  /*!
    @brief Constructor
  */
  air_bolsig();

  /*!
    @brief Destructor
  */
  virtual ~air_bolsig();

protected:

  /*!
    @brief Minimum electric field to compute transport coefficients from
  */
  Real m_bolsig_min_townsend;
  
  /*!
    @brief Maximum electric field to compute transport coefficients from
  */
  Real m_bolsig_max_townsend;

  /*!
    @brief Number of grid points to use for BOLSIG computation
  */
  int m_bolsig_grid_points;

  /*!
    @brief Call BOLSIG in order to read in electron transport data
  */
  virtual void call_bolsig();

  /*!
    @brief Build the BOLSIG input script
  */
  virtual void build_script();

  /*!
    @brief Let BOLSIG build the input script
  */
  virtual void build_data();

  /*!
    @brief Extract necessary data from BOLSIG output file
  */
  virtual void extract_data();

  /*!
    @brief Delete the BOLSIG output file
  */
  virtual void delete_Data();

  /*!
    @brief Delete the BOLSIG script
  */
  virtual void delete_script();

};
#endif
