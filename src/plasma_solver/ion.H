/*!
  @file ion.H
  @brief Declaration of an ion class needed for plasma_kinetics
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _ION_
#define _ION_

#include <RealVect.H>

/*!
  @brief Abstract class that describes an electron or ion
*/
class ion {
public:
  /*!
    @brief Weak construction with default parameters
  */
  ion();

  /*!
    @brief Destructor
  */
  ~ion();

  /*!
    @brief Initial data. Must be set. 
    @param[in] a_pos Position
    @param[in] a_time Time
  */
  virtual const Real initial_data(const RealVect a_pos, const Real a_time) const = 0;

  /*!
    @brief Return name
  */
  virtual const std::string get_name() const {
    return m_name;
  }

  /*!
    @brief Return charge
  */
  virtual const int get_charge() const {
    return m_charge;
  }

  /*!
    @brief Return diffusive or not
  */
  virtual const bool is_diffusive() const {
    return m_diffusive;
  }

protected:

  /*!
    @brief Ion name
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
