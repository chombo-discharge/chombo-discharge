/*!
  @file wall_bc.H
  @brief Indexing class for BCs on domain edges
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _WALL_BC_
#define _WALL_BC_

#include "type_definitions.H"

#include <LoHiSide.H>
#include <RealVect.H>

/*!
  @brief Types of wall boundary conditions for various solvers. I hate this structure. 
*/
namespace wallbc {
  enum which_bc {
    dirichlet     = 0,
    dirichlet_map = 1,
    neumann       = 2,
    neumann_map   = 3,
    robin         = 4,
    robin_map     = 5,
    outflow       = 6,
    symmetry      = 7,
    wall          = 8,
  };
};

/*!
  @brief Indexing class of BCs on domain edges. Very lightweight stuff.
*/
class wall_bc {
public:

  /*!
    @brief Create BC of a type
  */
  wall_bc(const int a_dir, const Side::LoHiSide a_side, wallbc::which_bc a_which);

  /*!
    @brief Destructor
  */
  ~wall_bc();

  /*!
    @brief Mapping of bcs
  */
  static int map_bc(const int a_dir, const Side::LoHiSide a_side);

  /*!
    @brief Direction
  */
  virtual int direction(){
    return m_dir;
  }

  /*!
    @brief Side
  */
  virtual Side::LoHiSide side(){
    return m_side;
  }

  /*!
    @brief Set value
  */
  virtual void set_value(Real a_value);

  /*!
    @brief Set live or not
  */
  virtual void set_live(bool a_live);

  /*!
    @brief Set function
  */
  virtual void set_function(Real (*m_func)(const RealVect a_pos));

  /*!
    @brief Check if domain wall is live or not
  */
  virtual bool is_live();

  /*!
    @brief Get value. 
  */
  virtual Real get_value();

  /*!
    @brief Get BC type
  */
  virtual wallbc::which_bc which_bc();

protected:

  /*!
    @brief Coordinate direction
  */
  int m_dir;

  /*!
    @brief Side (lo/hi)
  */
  Side::LoHiSide m_side;

  /*!
    @brief Value. This is for Neumann
  */
  Real m_value;

  /*!
    @brief Live or not. This is for Dirichlet.
  */
  bool m_live;

  /*!
    @brief Function based. This is for both Dirichlet and Neumann
  */
  Real (*m_func)(const RealVect a_pos);

  /*!
    @brief BC type
  */
  wallbc::which_bc m_which;
};

#endif
