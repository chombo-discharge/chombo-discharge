/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParallelOps.H
  @brief  Agglomeration of basic MPI reductions
  @author Robert Marskar
*/

#ifndef CD_ParallelOps_H
#define CD_ParallelOps_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief For putting a bunch of often-used MPI reductions in one place
*/
namespace ParallelOps {

  /*!
    @brief MPI barrier
  */
  inline void
  barrier() noexcept;

  /*!
    @brief Partition the input number such that each MPI rank gets an equal chunk. We assume
    that N is the length of some container whose indexing starts at 0.
    @details Returns start/stop indices along N. For example, if using two MPI ranks
    with N = 10 the first rank obtains pair(0,4) and the second rank get (5,9)
    @param[in] a_N Some container size.
  */
  inline std::pair<int, int>
  partition(const int a_N) noexcept;

  /*!
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
  */
  inline Real
  min(const Real& a_input) noexcept;

  /*!
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
  */
  inline int
  min(const int& a_input) noexcept;

  /*!
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
  */
  inline long long int
  min(const long long int& a_input) noexcept;

  /*!
    @brief Get the minimum value and the position corresponding to the minimum value of a_val
    @details This returns the value a_val and position a_pos such that they correspond to the 
    input variables for the processer which held the smallest value a_val. 
  */
  inline std::pair<Real, RealVect>
  min(const Real& a_val, const RealVect& a_pos) noexcept;

  /*!
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
  */
  inline Real
  max(const Real& a_input) noexcept;

  /*!
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
  */
  inline int
  max(const int& a_input) noexcept;

  /*!
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
  */
  inline long long
  max(const long long& a_input) noexcept;

  /*!
    @brief Get the maximum value and the position corresponding to the maximum value of a_val
    @details This returns the value a_val and position a_pos such that they correspond to the 
    input variables for the processer which held the largest value a_val. 
  */
  inline std::pair<Real, RealVect>
  max(const Real& a_val, const RealVect& a_pos) noexcept;

  /*!
    @brief Get the maximum value and the rank having the maximum value. 
    @param[in] a_val Value to be compared. 
  */
  inline std::pair<Real, int>
  maxRank(const Real& a_val) noexcept;

  /*!
    @brief Get the minimum value and the rank having the minimum value. 
    @param[in] a_val Value to be compared. 
  */
  inline std::pair<Real, int>
  minRank(const Real& a_val) noexcept;

  /*!
    @brief Compute the average (across MPI ranks) of the input value
    @param[in] a_val Value to take the average of
  */
  inline Real
  average(const Real& a_val) noexcept;

  /*!
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
  */
  inline Real
  sum(const Real& a_value) noexcept;

  /*!
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
  */
  inline int
  sum(const int& a_value) noexcept;

  /*!
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
  */
  inline long
  sum(const long& a_value) noexcept;

  /*!
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
  */
  inline long long
  sum(const long long& a_value) noexcept;

  /*!
    @brief Compute the sum across all MPI ranks
    @param[in] a_value Value to the take the sum of. 
  */
  inline unsigned long long
  sum(const unsigned long long& a_value) noexcept;

  /*!
    @brief Compute the standard deviation of the input value.
    @param[in] a_value Value to compute the standard deviation of.
  */
  inline Real
  standardDeviation(const Real& a_value) noexcept;

  /*!
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[inout] a_data On 
  */
  inline void
  vectorSum(Vector<Real>& a_data) noexcept;
  /*!
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[inout] a_data On 
  */
  inline void
  vectorSum(Vector<int>& a_data) noexcept;

  /*!
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[inout] a_data On 
  */
  inline void
  vectorSum(Vector<long int>& a_data) noexcept;

  /*!
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[inout] a_data On 
  */
  inline void
  vectorSum(Vector<long long int>& a_data) noexcept;
} // namespace ParallelOps

#include <CD_NamespaceFooter.H>

#include <CD_ParallelOpsImplem.H>

#endif
