/*!
  @file species.H
  @brief Declaration of a species class that passes in names and initial data to the cdr_solver
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _SPECIES_
#define _SPECIES_

#include <RealVect.H>

/*!
  @brief Abstract class that describes an electron or species
*/
class species {
public:
  /*!
    @brief Weak constructor. 
  */
  species();

  /*!
    @brief Full constructor
  */
  species(const std::string a_name, const int a_charge, const bool a_diffusive);

  /*!
    @brief Destructor
  */
  ~species();

  /*!
    @brief Initial data. You really should overwrite this one. 
    @param[in] a_pos Positspecies
    @param[in] a_time Time
  */
  virtual const Real initial_data(const RealVect a_pos, const Real a_time) const;

  /*!
    @brief Return name
  */
  virtual const std::string get_name() const;

  /*!
    @brief Return charge
  */
  virtual const int get_charge() const;
    
  /*!
    @brief Return diffusive or not
  */
  virtual const bool is_diffusive() const;

protected:

  /*!
    @brief Species name
  */
  std::string m_name;
  
  /*!
    @brief Charge
  */
  int m_charge;

  /*!
    @brief Diffusive species or not
  */
  bool m_diffusive;
};

#endif
