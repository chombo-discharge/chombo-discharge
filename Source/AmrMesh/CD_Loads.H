/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Loads.H
  @brief  Declaration of a simple class for holding computational loads
  @author Robert Marskar
*/

#ifndef CD_Loads_H
#define CD_Loads_H

// Std includes
#include <vector>
#include <map>

// Chombo includes
#include <Vector.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for holding computational loads. 
  @details This is typically used in the context of load-balancing with AMR. The class is, however, standalone.
*/
class Loads
{
public:
  /*!
    @brief Constructor. Assigns zero loads
  */
  Loads() noexcept;

  /*!
    @brief Destructor (does nothing).
  */
  virtual ~Loads() noexcept;

  /*!
    @brief Get the computational loads
    @return m_loads
  */
  virtual std::map<int, long int>&
  getLoads() noexcept;

  /*!
    @brief Get the computational loads
    @return m_loads
  */
  virtual const std::map<int, long int>&
  getLoads() const noexcept;

  /*!
    @brief Reset loads. Sets all loads to 0.
  */
  virtual void
  resetLoads() noexcept;

  /*!
    @brief Assign loads.
    @param[in] a_assignedLoads Loads assigned to ranks.
  */
  virtual void
  assignLoads(const std::map<int, long int>& a_assignedLoads) noexcept;

  /*!
    @brief Assign loads. In this version a_assignedLoads[0] is the load assigned to rank 0, and so on. 
    @param[in] a_assignedLoads Loads assigned to ranks.
  */
  virtual void
  assignLoads(const std::vector<long int>& a_assignedLoads) noexcept;

  /*!
    @brief Assign loads. In this version a_assignedLoads[0] is the load assigned to rank 0, and so on. 
    @param[in] a_assignedLoads Loads assigned to ranks.
  */
  virtual void
  assignLoads(const Vector<long int>& a_assignedLoads) noexcept;

  /*!
    @brief Increment loads. In this version a_increments[0] is the load increment for rank 0, and so on. 
    @param[in] a_increments Loads to increment by.
  */
  virtual void
  incrementLoads(const std::map<int, long int>& a_increments) noexcept;

  /*!
    @brief Increment loads. In this version a_increments[0] is the load increment for rank 0, and so on. 
    @param[in] a_increments Loads to increment by.
  */
  virtual void
  incrementLoads(const std::vector<long int>& a_increments) noexcept;

  /*!
    @brief Increment loads. In this version a_increments[0] is the load increment for rank 0, and so on. 
    @param[in] a_increments Loads to increment by.
  */
  virtual void
  incrementLoads(const Vector<long int>& a_increments) noexcept;

  /*!
    @brief Get sorted loads. 
    @details This routine sorts the computational loads from smallest to largest, and returns the result as a vector of ranks (and loads). 
  */
  virtual std::vector<std::pair<int, long int>>
  getSortedLoads() const noexcept;

protected:
  /*!
    @brief Actual computational loads per MPI rank. 
    @details The key is the rank, and the value is the actual load
  */
  std::map<int, long int> m_loads;
};

#include <CD_NamespaceFooter.H>

#endif
