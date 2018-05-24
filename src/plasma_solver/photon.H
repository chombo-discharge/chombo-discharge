/*!
  @file   photon.H
  @brief  Declaration of a class that holds information about optical transitions
  @author Robert Marskar
  @date   June 2016
*/

#ifndef _PHOTON_
#define _PHOTON_

#include "RealVect.H"

/*!
  @brief Abstract class that describes an optical transition
*/
class photon {
public:
  /*!
    @brief Weak construction with default parameters
  */
  photon();

  /*!
    @brief Destructor
  */
  ~photon();

  /*!
    @brief Get absorption coefficient at position in space
  */
  virtual const Real get_absorption_coeff(const RealVect& a_pos) const = 0;

  /*!
    @brief Return name
  */
  virtual const std::string get_name() const {
    return m_name;
  }

protected:
    
  /*!
    @brief Ion name
  */
  std::string m_name;  
};

#endif
