/*!
  @file wall_bc.H
  @brief Indexing class for BCs on domain edges
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _WALL_BC_
#define _WALL_BC_

#include <LoHiSide.H>
#include <RealVect.H>

namespace WallBC {
  enum WhichBC {
    Dirichlet = 0,
    Neumann = 1,
    Robin = 2,
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
  wall_bc(const int a_dir, const Side::LoHiSide a_side, WallBC::WhichBC a_which);

  /*!
    @brief Destructor
  */
  ~wall_bc();

  /*!
    @brief Mapping of bcs
  */
  static int map_bc(const int a_dir, const Side::LoHiSide a_side);

  /*!
    @brief Set value
  */
  virtual void set_value(Real a_value);

  /*!
    @brief Set live or not
  */
  virtual void set_live(bool a_live);

  /*!
    @brief Get BC type
  */
  virtual WallBC::WhichBC which_bc();

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
    @brief BC type
  */
  WallBC::WhichBC m_which;
};

#endif
