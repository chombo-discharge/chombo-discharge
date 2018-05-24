/*!
  @file load_balance.H
  @brief Declaration of various load-balancing routines
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _LOAD_BALANCE_
#define _LOAD_BALANCE_

#include "mfis.H"

/*!
  @brief Static class that does various kinds of load-balancing
*/
class load_balance {
public:

  /*!
    @brief Load-balancing types
  */
  enum which_balance {
    knapsack,
    elliptic,
    multifluid
  };

  /*!
    @brief Conventional load balancing
  */
  static void balance_knapsack(Vector<int>& a_procs, const Vector<Box>& a_boxes);
  
  /*!
    @brief Do elliptic load-balancing using EBAMRPoissonOp
    @param[out] a_procs   Processor ownership of boxes
    @param[in]  a_boxes   Boxes (Remember to order them by using space filling curves first!)
    @param[in]  a_ebis    Index space
    @param[in]  a_domain  Problem domain
    @param[in]  a_verbose Verbose
  */
  static void balance_elliptic(Vector<int>&                       a_procs,
			       const Vector<Box>&                 a_boxes,
			       const RefCountedPtr<EBIndexSpace>& a_ebis,
			       const ProblemDomain&               a_domain,
			       const bool                         a_verbose);

  /*!
    @brief Do elliptic load-balancing using EBAMRPoissonOp
    @param[out] a_procs   Processor ownership of boxes
    @param[in]  a_boxes   Boxes (Remember to order them by using space filling curves first!)
    @param[in]  a_mfis    EBIndex-spaces
    @param[in]  a_domain  Problem domain
    @param[in]  a_verbose Verbose
  */
  static void balance_multifluid(Vector<int>&               a_procs,
				 const Vector<Box>&         a_boxes,
				 const RefCountedPtr<mfis>& a_mfis,
				 const ProblemDomain&       a_domain,
				 const bool                 a_verbose);

protected:

  static void get_mfpoisson_loads(Vector<unsigned long long>& a_loads,
				  Vector<Box>&                a_boxes,
				  RefCountedPtr<mfis>&        a_mfis,
				  const DisjointBoxLayout&    a_dbl,
				  const ProblemDomain&        a_domain);
};
#endif
