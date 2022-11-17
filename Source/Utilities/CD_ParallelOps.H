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
  barrier();

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
  min(const Real& a_input);

  /*!
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
  */
  inline int
  min(const int& a_input);

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
  max(const Real& a_input);

  /*!
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
  */
  inline int
  max(const int& a_input);

  /*!
    @brief Get the maximum value and the position corresponding to the maximum value of a_val
    @details This returns the value a_val and position a_pos such that they correspond to the 
    input variables for the processer which held the largest value a_val. 
  */
  inline std::pair<Real, RealVect>
  max(const Real& a_val, const RealVect& a_pos) noexcept;

  /*!
    @brief Compute the sum across all MPI ranks
  */
  inline Real
  sum(const Real& a_value);

  /*!
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[inout] a_data On 
  */
  inline void
  vectorSum(Vector<Real>& a_data);

  /*!
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[inout] a_data On 
  */
  inline void
  vectorSum(Vector<long int>& a_data);
} // namespace ParallelOps

#include <CD_NamespaceFooter.H>

#include <CD_ParallelOpsImplem.H>

#endif
