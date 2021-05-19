/*!
  @file   air7.H
  @brief  Declaration of species for the BOLSIG-based 8-species air model
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air7.H"

#ifndef _AIR7_SPECIES_
#define _AIR7_SPECIES_

/*!
  @brief Electron class
*/
class air7::electron : public species {
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
  virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
protected:

  /*!
    @brief Initial ionization level
  */
  Real m_initial_ionization;
  Real m_seed_density;
  Real m_seed_radius;
  RealVect m_seed_position;
  
};

/*!
  @brief N2plus class
*/
class air7::N2plus : public species {
public:

  /*!
    @brief Constructor
  */
  N2plus();

  /*!
    @brief Destructor
  */
  virtual ~N2plus();

  /*!
    @brief Initial data
  */
  virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
protected:

  Real m_initial_ionization;
  Real m_seed_density;
  Real m_seed_radius;
  RealVect m_seed_position;
};

/*!
  @brief N4plus class
*/
class air7::N4plus : public species {
public:

  /*!
    @brief Constructor
  */
  N4plus();

  /*!
    @brief Destructor
  */
  virtual ~N4plus();

  /*!
    @brief Initial data
  */
  virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief O2plus class
*/
class air7::O2plus : public species {
public:

  /*!
    @brief Constructor
  */
  O2plus();

  /*!
    @brief Destructor
  */
  virtual ~O2plus();

  /*!
    @brief Initial data
  */
  virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
protected:

  Real m_initial_ionization;
  Real m_seed_density;
  Real m_seed_radius;
  RealVect m_seed_position;
};

/*!
  @brief O4plus class
*/
class air7::O4plus : public species {
public:

  /*!
    @brief Constructor
  */
  O4plus();

  /*!
    @brief Destructor
  */
  virtual ~O4plus();

  /*!
    @brief Initial data
  */
  virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief O2plusN2 class
*/
class air7::O2plusN2 : public species {
public:

  /*!
    @brief Constructor
  */
  O2plusN2();

  /*!
    @brief Destructor
  */
  virtual ~O2plusN2();

  /*!
    @brief Initial data
  */
  virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief O2minus class
*/
class air7::O2minus : public species {
public:

  /*!
    @brief Constructor
  */
  O2minus();

  /*!
    @brief Destructor
  */
  virtual ~O2minus();

  /*!
    @brief Initial data
  */
  virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief First Contribution to \gamma + O2 -> e + O2minus
*/
class air7::Photon_one : public Photon_group {
public:
  /*!
    @brief Constructor
  */
  Photon_one();

  /*!
    @brief Destructor
  */
  ~Photon_one();

  /*!
    @brief Get the absorption coefficient
  */
  virtual Real getKappa(const RealVect a_pos) const;

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
class air7::Photon_two : public Photon_group {
public:
  /*!
    @brief Constructor
  */
  Photon_two();

  /*!
    @brief Destructor
  */
  ~Photon_two();

  /*!
    @brief Get the absorption coefficient
  */
  virtual Real getKappa(const RealVect a_pos) const;

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
class air7::Photon_three : public Photon_group {
public:
  /*!
    @brief Constructor
  */
  Photon_three();

  /*!
    @brief Destructor
  */
  ~Photon_three();

  /*!
    @brief Get the absorption coefficient
  */
  virtual Real getKappa(const RealVect a_pos) const;

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
#include "CD_NamespaceFooter.H"

#endif
