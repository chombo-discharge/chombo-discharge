/*!
  @file   ito_species.H
  @brief  Declaration of an ito_species class that passes in names and initial data to ito_solvers
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _ITO_SPECIES_
#define _ITO_SPECIES_

#include <RealVect.H>
#include <Particle.H>
#include <BinItem.H>
#include <MeshInterp.H>

#include "ito_particle.H" 

/*!
  @brief Metadata class that parses for ito_solver
  @details The ito_solver takes it's information from this class, e.g. mobile/diffusive, charge, name, initial particles etc. 
*/
class ito_species {
public:


  /*!
    @brief Full constructor
  */
  ito_species(const std::string a_name, const int a_charge, const bool a_mobile, const bool a_diffusive);

  /*!
    @brief Destructor
  */
  ~ito_species();

  /*!
    @brief Return name
  */
  std::string get_name() const;

  /*!
    @brief Return charge
  */
  int get_charge() const;
    
  /*!
    @brief Return diffusive or not
  */
  bool is_diffusive() const;

  /*!
    @brief Mobile ito_species or not
  */
  bool is_mobile() const;

  /*!
    @brief Get initial particles
  */
  List<ito_particle>& get_initial_particles();

protected:

  /*!
    @brief Ito_Species name
  */
  std::string m_name;
  
  /*!
    @brief Charge
  */
  int m_charge;

  /*!
    @brief Diffusive ito_species or not
  */
  bool m_diffusive;

  /*!
    @brief Mobile ito_species or not
  */
  bool m_mobile;

  /*!
    @brief Initial particles
  */
  List<ito_particle> m_initial_particles;

  /*!
    @brief Weak constructor. 
  */
  ito_species();

};

#endif
