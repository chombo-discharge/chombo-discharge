/*!
  @file   air3_species.H
  @brief  Declaration of species for the air3 model
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air3.H"

#ifndef _AIR3_SPECIES_
#define _AIR3_SPECIES_

/*!
  @brief Electron class
*/
class air3::electron : public species {
public:

  /*!
    @brief Constructor
  */
  electron();

  /*!
    @brief Destructor
  */
  virtual ~electron();

  /*!
    @brief Initial data
  */
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:

  /*!
    @brief Initial ionization level
  */
  Real m_initial_ionization;
  
};

/*!
  @brief O2 class
*/
class air3::positive_species : public species {
public:

  /*!
    @brief Constructor
  */
  positive_species();

  /*!
    @brief Destructor
  */
  virtual ~positive_species();

  /*!
    @brief Initial data
  */
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:
  
  Real m_initial_ionization;
};

/*!
  @brief Negative species class
*/
class air3::negative_species : public species {
public:

  /*!
    @brief Constructor
  */
  negative_species();

  /*!
    @brief Destructor
  */
  virtual ~negative_species();

  /*!
    @brief Initial data
  */
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
};

/*!
  @brief First Contribution to \gamma + O2 -> e + O2minus
*/
class air3::photon_one : public photon_group {
public:
  /*!
    @brief Constructor
  */
  photon_one();

  /*!
    @brief Destructor
  */
  ~photon_one();

  /*!
    @brief Get the absorption coefficient
  */
  virtual Real get_kappa(const RealVect a_pos) const;

  /*!
    @brief Get the lambda coefficient
  */
  virtual Real get_lambda() const;

  /*!
    @brief Get the A-coefficient
  */
  virtual Real get_A() const;

  /*!
    @brief Get the partial oxygen pressure
  */
  virtual Real get_pO2() const;

protected:
  /*!
    @brief Lambda-coefficient
  */
  Real m_lambda;

  /*!
    @brief A-coefficient
  */
  Real m_A;

  /*!
    @brief O2-fraction
  */
  Real m_pO2;
};

/*!
  @brief First Contribution to \gamma + O2 -> e + O2minus
*/
class air3::photon_two : public photon_group {
public:
  /*!
    @brief Constructor
  */
  photon_two();

  /*!
    @brief Destructor
  */
  ~photon_two();

  /*!
    @brief Get the absorption coefficient
  */
  virtual Real get_kappa(const RealVect a_pos) const;

  /*!
    @brief Get the lambda coefficient
  */
  virtual Real get_lambda() const;

  /*!
    @brief Get the A-coefficient
  */
  virtual Real get_A() const;

  /*!
    @brief Get the partial oxygen pressure
  */
  virtual Real get_pO2() const;

protected:
  /*!
    @brief Lambda-coefficient
  */
  Real m_lambda;

  /*!
    @brief A-coefficient
  */
  Real m_A;

  /*!
    @brief O2-fraction
  */
  Real m_pO2;
};

/*!
  @brief First Contribution to \gamma + O2 -> e + O2minus
*/
class air3::photon_three : public photon_group {
public:
  /*!
    @brief Constructor
  */
  photon_three();

  /*!
    @brief Destructor
  */
  ~photon_three();

  /*!
    @brief Get the absorption coefficient
  */
  virtual Real get_kappa(const RealVect a_pos) const;

  /*!
    @brief Get the lambda coefficient
  */
  virtual Real get_lambda() const;

  /*!
    @brief Get the A-coefficient
  */
  virtual Real get_A() const;

  /*!
    @brief Get the partial oxygen pressure
  */
  virtual Real get_pO2() const;

protected:
  /*!
    @brief Lambda-coefficient
  */
  Real m_lambda;

  /*!
    @brief A-coefficient
  */
  Real m_A;

  /*!
    @brief O2-fraction
  */
  Real m_pO2;
};

#endif
