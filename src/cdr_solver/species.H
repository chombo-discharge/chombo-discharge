/*!
  @file species.H
  @brief Declaration of a species class that passes in names and initial data to the cdr_solver
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _SPECIES_
#define _SPECIES_

#include <RealVect.H>
#include <Particle.H>
#include <BinItem.H>
#include <MeshInterp.H>

/*!
  @brief Abstract class that describes an electron or species
*/
class species {
public:
  /*!
    @brief Weak constructor. 
  */
  species();

  /*!
    @brief Full constructor
  */
  species(const std::string a_name, const int a_charge, const bool a_mobile, const bool a_diffusive);

  /*!
    @brief Destructor
  */
  virtual ~species();

  /*!
    @brief Initial data. You really should overwrite this one. 
    @param[in] a_pos Positspecies
    @param[in] a_time Time
  */
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const = 0;

  /*!
    @brief Return name
  */
  virtual std::string get_name() const;

    /*!
    @brief Return unit
  */
  virtual std::string get_unit() const;

  /*!
    @brief Return charge
  */
  virtual int get_charge() const;
    
  /*!
    @brief Return diffusive or not
  */
  virtual bool is_diffusive() const;

  /*!
    @brief Mobile species or not
  */
  virtual bool is_mobile() const;

  /*!
    @brief Flag species for output to the plasma solver
  */
  virtual bool force_output() const;

  /*!
    @brief Init with particles or not
  */
  virtual bool init_with_particles() const;

  /*!
    @brief Init with function or not
  */
  virtual bool init_with_function() const;

  /*!
    @brief Get deposition
  */
  virtual InterpType& get_deposition();

  /*!
    @brief Get initial particles
  */
  virtual List<Particle>& get_initial_particles();

protected:

  /*!
    @brief Species name
  */
  std::string m_name;

  /*!
    @brief Unit
  */
  std::string m_unit;
  
  /*!
    @brief Charge
  */
  int m_charge;

  /*!
    @brief Diffusive species or not
  */
  bool m_diffusive;

  /*!
    @brief Mobile species or not
  */
  bool m_mobile;

  /*!
    @brief Force output of species. This is relevant for plasma_engine only. 
  */
  bool m_force_output;

  /*!
    @brief Init with particles
  */
  bool m_init_with_particles;

  /*!
    @brief Init function values
  */
  bool m_init_with_function;

  /*!
    @brief Initial particles
  */
  List<Particle> m_initial_particles;

  /*!
    @brief Deposition type for initial particles
  */
  InterpType m_deposition;
};

#endif
