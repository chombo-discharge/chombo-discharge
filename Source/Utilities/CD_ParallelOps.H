/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParallelOps.H
  @brief  Agglomeration of basic MPI reductions
  @author Robert Marskar
*/

#ifndef CD_PARALLELOPS_H
#define CD_PARALLELOPS_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief For putting a bunch of often-used MPI reductions in one place
*/
namespace ParallelOps {

  /**
    @brief MPI barrier
  */
  inline void
  barrier() noexcept;

  /**
    @brief Partition the input number such that each MPI rank gets an equal chunk. We assume
    that N is the length of some container whose indexing starts at 0.
    @details Returns start/stop indices along N. For example, if using two MPI ranks
    with N = 10 the first rank obtains pair(0,4) and the second rank get (5,9)
    @param[in] a_N Some container size.
    @return Return value
  */
  inline std::pair<int, int>
  partition(const int a_N) noexcept;

  /**
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
    @return Global minimum
  */
  inline Real
  min(const Real& a_input) noexcept;

  /**
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
    @return Global minimum
  */
  inline int
  min(const int& a_input) noexcept;

  /**
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
    @return Global minimum
  */
  inline long long int
  min(const long long int& a_input) noexcept;

  /**
    @brief Get the minimum value and the position corresponding to the minimum value of a_val
    @details This returns the value a_val and position a_pos such that they correspond to the 
    input variables for the processor which held the smallest value a_val. 
    @return Global minimum
    @param[in,out] a_val Description
    @param[in,out] a_pos Description
  */
  inline std::pair<Real, RealVect>
  min(const Real& a_val, const RealVect& a_pos) noexcept;

  /**
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
    @return Global maximum
  */
  inline Real
  max(const Real& a_input) noexcept;

  /**
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
    @return Global maximum
  */
  inline int
  max(const int& a_input) noexcept;

  /**
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
    @return Global maximum
  */
  inline long long
  max(const long long& a_input) noexcept;

  /**
    @brief Get the maximum value and the position corresponding to the maximum value of a_val
    @details This returns the value a_val and position a_pos such that they correspond to the 
    input variables for the processor which held the largest value a_val. 
    @return Global maximum
    @param[in,out] a_val Description
    @param[in,out] a_pos Description
  */
  inline std::pair<Real, RealVect>
  max(const Real& a_val, const RealVect& a_pos) noexcept;

  /**
    @brief Get the maximum value and the rank having the maximum value. 
    @param[in] a_val Value to be compared. 
    @return Return value
  */
  inline std::pair<Real, int>
  maxRank(const Real& a_val) noexcept;

  /**
    @brief Get the minimum value and the rank having the minimum value. 
    @param[in] a_val Value to be compared. 
    @return Return value
  */
  inline std::pair<Real, int>
  minRank(const Real& a_val) noexcept;

  /**
    @brief Compute the average (across MPI ranks) of the input value
    @param[in] a_val Value to take the average of
    @return Return value
  */
  inline Real
  average(const Real& a_val) noexcept;

  /**
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
    @return Global sum
  */
  inline Real
  sum(const Real& a_value) noexcept;

  /**
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
    @return Global sum
  */
  inline int
  sum(const int& a_value) noexcept;

  /**
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
    @return Global sum
  */
  inline long
  sum(const long& a_value) noexcept;

  /**
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
    @return Global sum
  */
  inline long long
  sum(const long long& a_value) noexcept;

  /**
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
    @return Global sum
  */
  inline unsigned long long
  sum(const unsigned long long& a_value) noexcept;

  /**
    @brief Compute the standard deviation of the input value.
    @param[in] a_value Value to compute the standard deviation of.
    @return Return value
  */
  inline Real
  standardDeviation(const Real& a_value) noexcept;

  /**
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[in,out] a_data On 
  */
  inline void
  sum(Vector<Real>& a_data) noexcept;

  /**
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[in,out] a_data On 
  */
  inline void
  sum(Vector<int>& a_data) noexcept;

  /**
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[in,out] a_data On 
  */
  inline void
  sum(Vector<long int>& a_data) noexcept;

  /**
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[in,out] a_data On 
  */
  inline void
  sum(Vector<long long int>& a_data) noexcept;

  /**
    @brief Gather local values -- return a vector of local values.
    @param[in] a_localValue Local value per MPI rank
    @return Return value
  */
  inline Vector<int>
  gather(const int& a_localValue) noexcept;

  /**
    @brief Gather local values -- return a vector of local values.
    @param[in] a_localValue Local value per MPI rank
    @return Return value
  */
  inline Vector<Real>
  gather(const Real& a_localValue) noexcept;
} // namespace ParallelOps

#include <CD_NamespaceFooter.H>

#include <CD_ParallelOpsImplem.H>

#endif
