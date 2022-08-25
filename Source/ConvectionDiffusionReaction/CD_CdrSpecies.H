/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSpecies.H
  @brief  Declaration of a class that passes information into CdrSolver (e.g., initial conditions). 
  @author Robert Marskar
*/

#ifndef CD_CdrSpecies_H
#define CD_CdrSpecies_H

// Chombo includes
#include <List.H>
#include <RealVect.H>

// Our includes
#include <CD_PointParticle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract class that describes the scalar quantity that is advected in CdrSolver. 
  @details Instantiations of CdrSpecies (i.e., derived classes thereof) can be passed into CdrSolver in order to provide run-time information
  to the class. The user will typically set a name for the advected quantity, the chargeNumber (if there is one), and whether or not advection 
  and diffusion is turned on or off. 
*/
class CdrSpecies
{
public:
  /*!
    @brief Weak constructor. 
  */
  CdrSpecies();

  /*!
    @brief Full constructor
    @param[in] a_name         Name of the species
    @param[in] a_chargeNumber Charge number
    @param[in] a_isMbile      Mobile or not
    @param[in] a_isDiffusive  Diffusive or not
  */
  CdrSpecies(const std::string a_name, const int a_chargeNumber, const bool a_isMobile, const bool a_isDiffusive);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~CdrSpecies();

  /*!
    @brief Initial data. 
    @param[in] a_pos  Position.
    @param[in] a_time Time
  */
  virtual Real
  initialData(const RealVect a_pos, const Real a_time) const = 0;

  /*!
    @brief Return name
  */
  virtual std::string
  getName() const;

  /*!
    @brief Return charge
  */
  virtual int
  getChargeNumber() const;

  /*!
    @brief Diffusive species or not
  */
  virtual bool
  isDiffusive() const;

  /*!
    @brief Mobile species or not
  */
  virtual bool
  isMobile() const;

  /*!
    @brief Get initial particles
  */
  virtual const List<PointParticle>&
  getInitialParticles() const;

  /*!
    @brief Get initial particles
  */
  List<PointParticle>&
  getInitialParticles();

protected:
  /*!
    @brief Cdr_Species name
  */
  std::string m_name;

  /*!
    @brief Charge
  */
  int m_chargeNumber;

  /*!
    @brief Diffusive CdrSpecies or not
  */
  bool m_isDiffusive;

  /*!
    @brief Mobile CdrSpecies or not
  */
  bool m_isMobile;

  /*!
    @brief Initial particles
  */
  List<PointParticle> m_initialParticles;
};

#include <CD_NamespaceFooter.H>

#endif
