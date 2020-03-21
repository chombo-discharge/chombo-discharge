/*!
  @file   geo_coarsener.H
  @brief  Declaration of a utility class for removing geometrically tagged cells
  @author Robert Marskar
  @date   June 2018
*/

#ifndef _GEO_COARSENER_
#define _GEO_COARSENER_

#include "real_box.H"

#include <Vector.H>
#include <IntVectSet.H>

/*!
  @brief Utility class for removing geometrically tagged cells
  @details This class is an optimization class for coarsening of EBs. The usage of this class is simple; the user defines m_coarsen_boxes and a level for each coarsening box. These are passed to plasma_engine which will remove geometric tags in those regions (above the level specified)
*/
class geo_coarsener {
public:

  /*!
    @brief Constructor
  */
  geo_coarsener();

  /*!
    @brief Destructor
  */
  virtual ~geo_coarsener();

  /*!
    @brief Coarsen tags
  */
  virtual void coarsen_tags(Vector<IntVectSet>& a_tags, const Vector<Real>& a_dx, const RealVect& a_origin) const;

  /*!
    @brief Get coarsening boxes
  */
  virtual Vector<real_box> get_coarsen_boxes();

  /*!
    @brief Get coarsening levels
  */
  virtual Vector<int> get_coarsen_levels();
  
protected:

  /*!
    @brief Coarsening boxes
  */
  Vector<real_box> m_coarsen_boxes;

  /*!
    @brief Coarsening levels. Tags are removed above this level
  */
  Vector<int> m_coarsen_levels;

  /*!
    @brief Inverse box or not
  */
  Vector<int> m_inverse;
};
#endif
