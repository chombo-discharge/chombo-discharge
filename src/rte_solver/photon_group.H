/*!
  @file   photon_group.H
  @brief  Declaration of a photon-group class. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _PHOTON_GROUP_
#define _PHOTON_GROUP_

#include "RealVect.H"

/*!
  @brief Declaration of a class that defines a group within the multigroup approximation
  @details This is a lightweight class - it just holds a name for the group and the kappa for the group
*/
class photon_group {
public:
  /*!
    @brief Weak constructor, must subsequently call define
  */
  photon_group();

  /*!
    @brief Constructor, constant kappa
  */
  photon_group(const std::string a_name, const Real a_kappa);

  /*!
    @brief Constructor, spatially dependent kappa
  */
  photon_group(const std::string a_name, Real (*a_kappa)(const RealVect a_pos));

  /*!
    @brief Define function
  */
  virtual void define(const std::string a_name, const Real a_kappa);

  /*!
    @brief Define function
  */
  virtual void define(const std::string a_name, Real (*a_kappa)(const RealVect a_pos));

  /*!
    @brief Get group name
  */
  virtual std::string get_name() const;

  /*!
    @brief Get kappa at position. You really need to impement this, but it's not pure since this class may be useful for testing. 
  */
  virtual Real get_kappa(const RealVect a_pos) const;
	      
protected:

  /*!
    @brief Group name
  */
  std::string m_name;
  
  /*!
    @brief Constant kappa or not
  */
  bool m_constant;

  /*!
    @brief Constant kappa
  */
  Real m_kappa;

  /*!
    @brief Variable kappa
  */
  Real (*m_varkappa)(const RealVect a_pos);
};
#endif
  
