/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCSingleState.H
  @brief  Declaration of a simple state vector for running Kinetic Monte Carlo for plasma problems. 
  @author Robert Marskar
*/

#ifndef CD_KMCSingleState_H
#define CD_KMCSingleState_H

// Std includes
#include <vector>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a straightforward state for advancing with the Kinetic Monte Carlo module. 
  @details The state is simply a vector of integers (default long long type ints). 
  @note The template parameter indicates the integer type used for tracking the state. 
*/
template <typename T = long long>
class KMCSingleState
{
public:
  using State = std::vector<T>;

  /*!
    @brief Disallowed weak construction. 
  */
  inline KMCSingleState() = delete;

  /*!
    @brief Copy constructor.
  */
  inline KMCSingleState(KMCSingleState&) = default;

  /*!
    @brief Disallowed move constructor
  */
  inline KMCSingleState(KMCSingleState&&) = delete;

  /*!
    @brief Full constructor. 
    @param[in] a_numSpecies Number of reactive species in the allocated state vector.
  */
  inline KMCSingleState(const size_t a_numSpecies) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~KMCSingleState();

  /*!
    @brief Copy assignment constructor
  */
  inline KMCSingleState&
  operator=(const KMCSingleState&) = default;

  /*!
    @brief Disallowed move assignement
  */
  inline KMCSingleState&
  operator=(const KMCSingleState&&) = delete;

  /*!
    @brief Get the population of the input index
    @param[in] a_idx Index in state vector
  */
  inline T&
  operator[](const size_t a_idx) noexcept;

  /*!
    @brief Get the population of the input index
    @param[in] a_idx Index in state vector
  */
  inline const T&
  operator[](const size_t a_idx) const noexcept;

  /*!
    @brief Check if state is a valid state. 
    @return Returns false if any populations are negative. 
  */
  inline bool
  isValidState() const noexcept;

  /*!
    @brief Get modifiable state. 
    @return m_state
  */
  inline State&
  getState() noexcept;

  /*!
    @brief Get state. 
    @return m_state
  */
  inline const State&
  getState() const noexcept;

protected:
  /*!
    @brief State vector.
  */
  State m_state;
};

#include <CD_NamespaceFooter.H>

#include <CD_KMCSingleStateImplem.H>

#endif
