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
    @brief BC type
  */
  wallbc::which_bc m_which;
};

#endif
