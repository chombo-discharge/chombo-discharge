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
    volume,
    elliptic,
    multifluid
  };

  /*!
    @brief Load balancing. Loads are box volumes. 
  */
  static void make_balance(Vector<int>& a_levelRanks, const Vector<Box>& a_levelBoxes);

  /*!
    @brief Load balancing
  */
  template <class T>
  static void make_balance(Vector<int>& a_ranks, const Vector<T>& a_loads, const Vector<Box>& a_boxes);

  /*!
    @brief Load balancing function which load balances level by level. 
  */
  template <class T>
  static void level_by_level(Vector<Vector<int> >&       a_procs,
			     const Vector<Vector<T> >&   a_loads,
			     const Vector<Vector<Box> >& a_boxes);

  /*!
    @brief Load balancing function which load balances all boxes at once, not caring about actual levels. 
  */
  template <class T>
  static void hierarchy(Vector<Vector<int> >&       a_procs,
			const Vector<Vector<T> >&   a_loads,
			const Vector<Vector<Box> >& a_boxes);

  /*!
    @brief Gather boxes
  */
  static void gather_boxes(Vector<Box>& a_boxes);

  /*!
    @brief Gather compute loads
  */
  static void gather_loads(Vector<Real>& a_loads);

  /*!
    @brief Gather compute loads
  */
  static void gather_loads(Vector<int>& a_loads);

  /*!
    @brief Gather boxes and loads
  */
  static void gather_boxes_and_loads(Vector<Box>& a_boxes, Vector<int>& a_loads);

protected:
};

#include "load_balanceI.H"

#endif
