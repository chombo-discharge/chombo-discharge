/*!
  @file   rte_species.H
  @brief  Declaration of a photon-group class. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _RTE_SPECIES_
#define _RTE_SPECIES_

#include "RealVect.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Declaration of a class that defines a group within the multigroup approximation
  @details This is a lightweight class - it just holds a name for the group and the kappa for the group
*/
class rte_species {
public:
  /*!
    @brief Weak constructor, must subsequently call define
  */
  rte_species();

  /*!
    @brief Constructor, constant kappa
  */
  rte_species(const std::string a_name, const Real a_kappa);

  /*!
    @brief Constructor, spatially dependent kappa
  */
  rte_species(const std::string a_name, Real (*a_kappa)(const RealVect a_pos));

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
    @brief Constant kappa or not
  */
  virtual bool constant_kappa() const;

  /*!
    @brief Get kappa at position. You really need to impement this, but it's not pure since this class may be useful for testing. 
  */
  virtual Real get_kappa(const RealVect a_pos) const;

  /*!
    @brief Get scattering coefficient. 
  */
  virtual Real get_scatter(const RealVect a_pos) const;
	      
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
    @brief Scatter coefficient
  */
  Real m_scatter;

  /*!
    @brief Variable kappa
  */
  Real (*m_varkappa)(const RealVect a_pos);

  /*!
    @brief Variable kappa
  */
  Real (*m_varscatter)(const RealVect a_pos);
};
#include "CD_NamespaceFooter.H"
#endif
  
