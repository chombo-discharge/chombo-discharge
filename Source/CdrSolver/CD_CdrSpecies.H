/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrSpecies.H
  @brief  Declaration of a CdrSpecies class that passes in names and initial data to the CdrSolver
  @author Robert Marskar
*/

#ifndef CD_CdrSpecies_H
#define CD_CdrSpecies_H

// Chombo includes
#include <RealVect.H>
#include <Particle.H>
#include <EBParticleInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Abstract class that describes an electron or CdrSpecies
*/
class CdrSpecies {
public:
  /*!
    @brief Weak constructor. 
  */
  CdrSpecies();

  /*!
    @brief Full constructor
    @param[in] a_name         Name of the species
    @param[in] a_chargeNumber Charge number
    @param[in] a_mobile       Mobile or not
    @param[in] a_diffusive    Diffusive or not
  */
  CdrSpecies(const std::string a_name, const int a_chargeNumber, const bool a_mobile, const bool a_diffusive);

  /*!
    @brief Destructor
  */
  virtual ~CdrSpecies();

  /*!
    @brief Initial data. 
    @param[in] a_pos  Position.
    @param[in] a_time Time
  */
  virtual Real initialData(const RealVect a_pos, const Real a_time) const = 0;

  /*!
    @brief Return name
  */
  virtual std::string getName() const;

  /*!
    @brief Return unit
  */
  virtual std::string get_unit() const;

  /*!
    @brief Return charge
  */
  virtual int getChargeNumber() const;
    
  /*!
    @brief Return diffusive or not
  */
  virtual bool isDiffusive() const;

  /*!
    @brief Mobile CdrSpecies or not
  */
  virtual bool isMobile() const;

  /*!
    @brief Flag CdrSpecies for output to the plasma solver
  */
  virtual bool forceOutput() const;

  /*!
    @brief Init with particles or not
  */
  virtual bool initializeWithParticles() const;

  /*!
    @brief Init with function or not
  */
  virtual bool initializeWithFunction() const;

  /*!
    @brief Get deposition
  */
  virtual DepositionType::Which getDeposition();

  /*!
    @brief Get initial particles
  */
  virtual List<Particle>& getInitialParticles();

protected:

  /*!
    @brief Cdr_Species name
  */
  std::string m_name;

  /*!
    @brief Unit
  */
  std::string m_unit;
  
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
    @brief Force output of CdrSpecies. This is relevant for plasma_engine only. 
  */
  bool m_forceOutput;

  /*!
    @brief Init with particles
  */
  bool m_initializeWithParticles;

  /*!
    @brief Init function values
  */
  bool m_initializeWithFunction;

  /*!
    @brief Initial particles
  */
  List<Particle> m_initialParticles;

  /*!
    @brief Deposition type for initial particles
  */
  DepositionType::Which m_deposition;
};

#include <CD_NamespaceFooter.H>

#endif
