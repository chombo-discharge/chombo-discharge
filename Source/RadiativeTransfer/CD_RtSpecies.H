/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtSpecies.H
  @brief  Declaration of a class which supplies a user interface to radiative transfer code. 
  @author Robert Marskar
*/

#ifndef CD_RtSpecies_H
#define CD_RtSpecies_H

// Std includes
#include <functional>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a class that defines an interface to radiative transfer solvers, i.e. RtSolver. 
  @details This is a lightweight class - it just holds a name for the species as well as a function signature for obtaining the absorption and 
  scattering coefficients. Note that this species supports spatially varying absorption and scattering coefficients. 
*/
class RtSpecies {
public:
  
  /*!
    @brief Weak constructor. Sets absorption length and scattering coefficients to one. 
  */
  RtSpecies();

  /*!
    @brief Weak constructor. Sets absorption length and scattering coefficients to one. 
  */
  virtual ~RtSpecies();  

  /*!
    @brief Get species name.
  */
  std::string getName() const;

  /*!
    @brief Get kappa (i.e. the inverse absorption length) at physical coordinates. 
    @param[in] a_pos Physical coordinates. 
  */
  virtual Real getAbsorptionCoefficient(const RealVect a_pos) const = 0;

  /*!
    @brief Get scattering coefficient (i.e. the inverse scattering length) at physical coordinates. 
    @param[in] a_pos Physical coordinates. 
    @note This is not used anywhere (yet).
  */
  virtual Real getScatteringCoefficient(const RealVect a_pos) const;
	      
protected:

  /*!
    @brief Group name
  */
  std::string m_name;
};

#include <CD_NamespaceFooter.H>

#endif
  
