/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoSpecies.H
  @brief  Declaration of an ItoSpecies class that passes in names and initial data to ItoSolvers
  @author Robert Marskar
*/

#ifndef CD_ItoSpecies_H
#define CD_ItoSpecies_H

// Chombo includes
#include <List.H>

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that acts as user interface for parsing initial data and mobility/diffusion kernels into ItoSolver. 
  @details This class is used for passing initial data and information into an ItoSolver. This class is required because the ItoSolver may not know about
  whether or not it will use advect/diffusion, and it is also does not know about initial conditions. All of that is supplied by this class. 
*/
class ItoSpecies
{
public:
  /*!
    @brief Default constructor, but user should set the name, charge/mobility/diffusivity etc. 
  */
  ItoSpecies();

  /*!
    @brief Full constructor
    @param[in] a_name         Species name
    @param[in] a_chargeNumber Charge number
    @param[in] a_mobile       Mobile species or not
    @param[in] a_diffusive    Diffusive species or not
  */
  ItoSpecies(const std::string a_name, const int a_chargeNumber, const bool a_mobile, const bool a_diffusive);

  /*!
    @brief Destructor (does nothing):
  */
  virtual ~ItoSpecies();

  /*!
    @brief Return name
    @return Returns m_name
  */
  std::string
  getName() const;

  /*!
    @brief Return charge
    @return Returns m_chargeNumber
  */
  int
  getChargeNumber() const;

  /*!
    @brief Return the initial density
    @return Returns m_initialDensity
  */
  const std::function<Real(const RealVect& x, const Real& t)>&
  getInitialDensity() const;

  /*!
    @brief Return diffusive or not
  */
  bool
  isDiffusive() const;

  /*!
    @brief Mobile ItoSpecies or not
  */
  bool
  isMobile() const;

  /*!
    @brief Compute the mobility as a function of energy.
    @details This is mostly used for the LEA-type models where the particle mobility has a parametric dependence on the energy. 
    @param[in] a_energy Particle energy
    @note Default implementation returns 0.0
  */
  virtual Real
  mobility(const Real a_energy) const;

  /*!
    @brief Compute the diffusion coefficient.
    @details This is mostly used for the LEA-type models where the particle diffusivity has a parametric dependence on the energy. 
    @param[in] a_energy Particle energy
    @note Default implementation returns 0.0
  */
  virtual Real
  diffusion(const Real a_energy) const;

  /*!
    @brief Set the initial species density
    @param[in] a_initialDensity Initial density.
  */
  virtual void
  setInitialDensity(const std::function<Real(const RealVect& x, const Real& t)>& a_initialDensity);

  /*!
    @brief Get initial particles -- this is called by ItoSolver when filling the solver with initial particles. 
    @return Returns m_initialParticles
  */
  List<ItoParticle>&
  getInitialParticles();

  /*!
    @brief Get initial particles -- this is called by ItoSolver when filling the solver with initial particles. 
    @return Returns m_initialParticles
  */
  const List<ItoParticle>&
  getInitialParticles() const;

protected:
  /*!
    @brief Ito_Species name
  */
  std::string m_name;

  /*!
    @brief Charge
  */
  int m_chargeNumber;

  /*!
    @brief Maximum number of initial numerical particles per cell
  */
  int m_maxInitialParticlesPerCell;

  /*!
    @brief Diffusive ItoSpecies or not
  */
  bool m_isDiffusive;

  /*!
    @brief Mobile ItoSpecies or not
  */
  bool m_isMobile;

  /*!
    @brief Initial particles
  */
  List<ItoParticle> m_initialParticles;

  /*!
    @brief Initial density, in case the user wants to generate particles from a density distribution
  */
  std::function<Real(const RealVect& x, const Real& t)> m_initialDensity;
};

#include <CD_NamespaceFooter.H>

#endif
