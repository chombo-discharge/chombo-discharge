/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ItoSpecies.H
  @brief  Declaration of an ItoSpecies class that passes in names and initial data to ito_solvers
  @author Robert Marskar
*/

#ifndef CD_ItoSpecies_H
#define CD_ItoSpecies_H

// Chombo includes
#include <RealVect.H>
#include <Particle.H>
#include <BinItem.H>
#include <MeshInterp.H>

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Metadata class that parses for ito_solver
  @details The ito_solver takes it's information from this class, e.g. mobile/diffusive, charge, name, initial particles etc. 
*/
class ItoSpecies {
public:

  /*!
    @brief Weak constructor. 
  */
  ItoSpecies();

  /*!
    @brief Full constructor
  */
  ItoSpecies(const std::string a_name, const int a_charge, const bool a_mobile, const bool a_diffusive);

  /*!
    @brief Destructor
  */
  ~ItoSpecies();

  /*!
    @brief Return name
  */
  std::string getName() const;

  /*!
    @brief Return charge
  */
  int getChargeNumber() const;
    
  /*!
    @brief Return diffusive or not
  */
  bool isDiffusive() const;

  /*!
    @brief Mobile ItoSpecies or not
  */
  bool isMobile() const;

  /*!
    @brief Compute the mobility. This is mostly used for the LEA. 
  */
  virtual Real mobility(const Real a_energy) const;

  /*!
    @brief Compute the mobility. This is mostly used for the LEA.  
  */
  virtual Real diffusion(const Real a_energy) const;

  /*!
    @brief Get initial particles
  */
  List<ItoParticle>& getInitialParticles();

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
};

#include <CD_NamespaceFooter.H>

#endif
