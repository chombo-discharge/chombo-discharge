/*!
  @file   air9eed_bourdon.H
  @brief  Declaration of species for the BOLSIG-based 11-species air model
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air9eed_bourdon.H"

#ifndef _AIR9EED_BOURDON_SPECIES_
#define _AIR9EED_BOURDON_SPECIES_

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace cdr_plasma {


    /*!
      @brief Electron energy distribution
    */
    class air9eed_bourdon::eed : public CdrSpecies {
    public:

      /*!
	@brief Constructor
      */
      eed();

      /*!
	@brief Destructor
      */
      virtual ~eed();

      /*!
	@brief Initial data
      */
      virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
    protected:

      /*!
	@brief Gas temperature
      */
      Real m_init_energy;
    };

    /*!
      @brief Electron class
    */
    class air9eed_bourdon::electron : public CdrSpecies {
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
  
    };

    /*!
      @brief N2plus class
    */
    class air9eed_bourdon::N2plus : public CdrSpecies {
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
    };

    /*!
      @brief N4plus class
    */
    class air9eed_bourdon::N4plus : public CdrSpecies {
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
    class air9eed_bourdon::O2plus : public CdrSpecies {
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
    };

    /*!
      @brief O4plus class
    */
    class air9eed_bourdon::O4plus : public CdrSpecies {
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
    class air9eed_bourdon::O2plusN2 : public CdrSpecies {
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
    class air9eed_bourdon::O2minus : public CdrSpecies {
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
      @brief Ominus class
    */
    class air9eed_bourdon::Ominus : public CdrSpecies {
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
      virtual Real initialData(const RealVect a_pos, const Real a_time) const;
  
    protected:
    };


    /*!
      @brief First Contribution to \gamma + O2 -> e + O2minus
    */
    class air9eed_bourdon::photon_one : public rte_species {
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
    class air9eed_bourdon::photon_two : public rte_species {
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
    class air9eed_bourdon::photon_three : public rte_species {
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
  }
}
#include "CD_NamespaceFooter.H"

#endif
