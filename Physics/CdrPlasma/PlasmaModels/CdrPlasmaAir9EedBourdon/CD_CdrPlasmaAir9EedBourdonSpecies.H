/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir9EedBourdonSpecies.H
  @brief  Declaration of species for the CdrPlasmaAir9EedBourdon model
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaAir9EedBourdonSpecies_H
#define CD_CdrPlasmaAir9EedBourdonSpecies_H

// Our includes
#include <CD_CdrPlasmaAir9EedBourdon.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Electron energy distribution
    */
    class CdrPlasmaAir9EedBourdon::eed : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      /*!
	@brief Gas temperature
      */
      Real m_init_energy;
    };

    /*!
      @brief Electron class
    */
    class CdrPlasmaAir9EedBourdon::Electron : public CdrSpecies
    {
    public:
      /*!
	@brief Constructor
      */
      Electron();

      /*!
	@brief Destructor
      */
      virtual ~Electron();

      /*!
	@brief Initial data
      */
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      /*!
	@brief Initial ionization level
      */
      Real m_initial_ionization;
    };

    /*!
      @brief N2plus class
    */
    class CdrPlasmaAir9EedBourdon::N2plus : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      Real m_initial_ionization;
    };

    /*!
      @brief N4plus class
    */
    class CdrPlasmaAir9EedBourdon::N4plus : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
    };

    /*!
      @brief O2plus class
    */
    class CdrPlasmaAir9EedBourdon::O2plus : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      Real m_initial_ionization;
    };

    /*!
      @brief O4plus class
    */
    class CdrPlasmaAir9EedBourdon::O4plus : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
    };

    /*!
      @brief O2plusN2 class
    */
    class CdrPlasmaAir9EedBourdon::O2plusN2 : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
    };

    /*!
      @brief O2minus class
    */
    class CdrPlasmaAir9EedBourdon::O2minus : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
    };

    /*!
      @brief Ominus class
    */
    class CdrPlasmaAir9EedBourdon::Ominus : public CdrSpecies
    {
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
      virtual Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
    };

    /*!
      @brief First Contribution to \gamma + O2 -> e + O2minus
    */
    class CdrPlasmaAir9EedBourdon::PhotonOne : public RtSpecies
    {
    public:
      /*!
	@brief Constructor
      */
      PhotonOne();

      /*!
	@brief Destructor
      */
      ~PhotonOne();

      /*!
	@brief Get the absorption coefficient
      */
      virtual Real
      getAbsorptionCoefficient(const RealVect a_pos) const;

      /*!
	@brief Get the lambda coefficient
      */
      virtual Real
      getLambda() const;

      /*!
	@brief Get the A-coefficient
      */
      virtual Real
      getA() const;

      /*!
	@brief Get the partial oxygen pressure
      */
      virtual Real
      getPO2() const;

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
    class CdrPlasmaAir9EedBourdon::PhotonTwo : public RtSpecies
    {
    public:
      /*!
	@brief Constructor
      */
      PhotonTwo();

      /*!
	@brief Destructor
      */
      ~PhotonTwo();

      /*!
	@brief Get the absorption coefficient
      */
      virtual Real
      getAbsorptionCoefficient(const RealVect a_pos) const;

      /*!
	@brief Get the lambda coefficient
      */
      virtual Real
      getLambda() const;

      /*!
	@brief Get the A-coefficient
      */
      virtual Real
      getA() const;

      /*!
	@brief Get the partial oxygen pressure
      */
      virtual Real
      getPO2() const;

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
    class CdrPlasmaAir9EedBourdon::PhotonThree : public RtSpecies
    {
    public:
      /*!
	@brief Constructor
      */
      PhotonThree();

      /*!
	@brief Destructor
      */
      ~PhotonThree();

      /*!
	@brief Get the absorption coefficient
      */
      virtual Real
      getAbsorptionCoefficient(const RealVect a_pos) const;

      /*!
	@brief Get the lambda coefficient
      */
      virtual Real
      getLambda() const;

      /*!
	@brief Get the A-coefficient
      */
      virtual Real
      getA() const;

      /*!
	@brief Get the partial oxygen pressure
      */
      virtual Real
      getPO2() const;

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
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
