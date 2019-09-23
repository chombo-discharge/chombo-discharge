/*!
  @brief  fast_gshop.H
  @brief  Version of GeometryShop that also may use a bounding box for location EB volumes
  @author Robert Marskar
  @date   2019
*/

#ifndef _FAST_GSHOP_
#define _FAST_GSHOP_

#include <GeometryShop.H>
#include <GeometryService.H>

#include "real_box.H"

/*!
  @brief Derived version of GeometryShop that can also use a bounding box for EB volumes
*/
class fast_gshop : public GeometryShop {
public:

  fast_gshop(const BaseIF&       a_localGeom,
	     const int           a_verbosity,
	     const Real          a_dx,
	     const RealVect      a_origin,
	     const ProblemDomain a_finest_domain,
	     const Real          a_thrshdVoF = 1.0e-16);
  ~fast_gshop();
    
  void makeGrids(const ProblemDomain&      a_domain,
		 DisjointBoxLayout&        a_grids,
		 const int&                a_maxGridSize,
		 const int&                a_maxIrregGridSize);

  void set_bboxes(const Vector<real_box> a_bboxes);

  void set_regular_boxes(const Vector<real_box> a_regboxes);

  void set_covered_boxes(const Vector<real_box> a_regboxes);

  void add_bbox(const real_box a_rbox);

  void add_regular_box(const real_box a_rbox);

  void add_covered_box(const real_box a_rbox);

  GeometryService::InOut InsideOutside(const Box&           a_region,
				       const ProblemDomain& a_domain,
				       const RealVect&      a_origin,
				       const Real&          a_dx) const;




protected:

  RealVect              m_origin;  // Origin
  Vector<Real>          m_dx;      // Index 0 is the finest EBIS level
  Vector<ProblemDomain> m_domains; // Index 0 is the finest EBIS level

  Vector<real_box> m_bboxes;         // Bou
  Vector<real_box> m_covered_boxes; // Boxes that are known to be fully regular
  Vector<real_box> m_regular_boxes; // Boxes that are known to be fully covered

  /*!
    @brief Check if a box is regular, i.e. lies COMPLETELY within one of m_regular_boxes
  */
  bool is_regular_box(const real_box& a_box) const;

  /*!
    @brief Check if a box is fully covered, i.e. lies COMPLETELY within one of m_covered_boxes
  */
  bool is_covered_box(const real_box& a_box) const;

  /*!
    @brief Check if a box intersects a bounding box
  */
  bool is_bounding_box(const real_box& a_box) const;
};

#endif
