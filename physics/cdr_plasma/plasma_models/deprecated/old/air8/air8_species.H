/*!
  @file   air8.H
  @brief  Declaration of species for the BOLSIG-based 8-species air model
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air8.H"

#ifndef _AIR8_SPECIES_
#define _AIR8_SPECIES_

/*!
  @brief Electron class
*/
class air8::electron : public species {
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
  @brief N2plus class
*/
class air8::N2plus : public species {
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
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:

  Real m_initial_ionization;
};

/*!
  @brief N4plus class
*/
class air8::N4plus : public species {
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
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief O2plus class
*/
class air8::O2plus : public species {
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
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:

  Real m_initial_ionization;
};

/*!
  @brief O4plus class
*/
class air8::O4plus : public species {
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
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief O2plusN2 class
*/
class air8::O2plusN2 : public species {
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
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief O2minus class
*/
class air8::O2minus : public species {
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
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:
};

/*!
  @brief Ominus class
*/
class air8::Ominus : public species {
public:

  /*!
    @brief Constructor
  */
  Ominus();

  /*!
    @brief Destructor
  */
  virtual ~Ominus();

  /*!
    @brief Initial data
  */
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const;
  
protected:
};


/*!
  @brief First Contribution to \gamma + O2 -> e + O2minus
*/
class air8::photon_one : public photon_group {
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
class air8::photon_two : public photon_group {
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
class air8::photon_three : public photon_group {
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
