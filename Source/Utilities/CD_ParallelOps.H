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
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
  */
  inline
  Real Min(const Real& a_input);

  /*!
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
  */
  inline
  Real Max(const Real& a_input);

  /*!
    @brief Get the minimum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. Returns the smallest. 
  */
  inline
  int Min(const int& a_input);

  /*!
    @brief Get the maximum of the input, reduced over MPI ranks (in the Chombo communicator)
    @param[in] a_input Data input. 
  */
  inline
  int Max(const int& a_input);

  /*!
    @brief Compute the sum across all MPI ranks
  */
  inline
  Real sum(const Real& a_value);

  /*!
    @brief Perform a summation of all the MPI ranks's input data. 
    @details This performs a rank-wise summation. If rank 1 has data (1,2,3) and rank 2 has data (3,4,5), the output data
    on both ranks is (4,6,8).
    @param[inout] a_data On 
  */
  inline
  void VectorSum(Vector<long int>& a_data);
}

#include <CD_NamespaceFooter.H>

#include <CD_ParallelOpsImplem.H>

#endif
