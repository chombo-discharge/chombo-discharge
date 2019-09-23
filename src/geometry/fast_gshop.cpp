/*!
  @brief  fast_gshop.cpp
  @brief  Implementation of fast_gshop
  @author Robert Marskar
  @date   2019
*/

#include "fast_gshop.H"

#include <BRMeshRefine.H>
#include <LoadBalance.H>

#include <ParmParse.H>

fast_gshop::fast_gshop(const BaseIF&       a_localGeom,
		       const int           a_verbosity,
		       const Real          a_dx,
		       const RealVect      a_origin,
		       const ProblemDomain a_finest_domain,
		       const Real          a_thrshdVoF)
  : GeometryShop(a_localGeom, a_verbosity, a_dx*RealVect::Unit, a_thrshdVoF) {

  // EBISLevel doesn't give resolution, origin, and problem domains through makeGrids, so we
  // need to construct these here, and then extract the proper resolution when we actually do makeGrids
  
  m_origin = a_origin;
  
  m_dx.resize(0);
  m_domains.resize(0);
  
  m_domains.push_back(a_finest_domain);
  m_dx.push_back(a_dx);
  
  const int ref = 2;
  for (int lvl=0; ; lvl++){
    Real  dx             = m_dx[lvl];
    ProblemDomain domain = m_domains[lvl];
    
    if(domain.domainBox().coarsenable(ref)){
      domain.coarsen(ref);
      dx *= 1./ref;
      
      m_dx.push_back(dx);
      m_domains.push_back(domain);
    }
    else{
      break;
    }
  }


  m_regular_boxes.resize(0);
  m_covered_boxes.resize(0);
  
#if 1 // Debug code
  ParmParse pp("devel");
  int use_bbox = 0;
  pp.get("use_bbox", use_bbox);
  if(use_bbox != 0){



    // Evertying to regular
    m_regular_boxes.push_back(real_box(-RealVect::Unit, RealVect::Unit));

    // Except this boxc

    // Everything below needle
    Real radius = 1.E-3;
    RealVect lo = RealVect::Zero;
    RealVect hi = RealVect(0, 0, 1);

    lo = lo - radius*1.1*RealVect::Unit;
    hi = hi + radius*1.1*RealVect::Unit;

    m_bboxes.push_back(real_box(lo, hi));
  }

#endif
}

fast_gshop::~fast_gshop(){

}


void fast_gshop::set_bboxes(const Vector<real_box> a_bboxes){
  m_bboxes = a_bboxes;
}

void fast_gshop::set_regular_boxes(const Vector<real_box> a_regular_boxes){
  m_regular_boxes = a_regular_boxes;
}

void fast_gshop::set_covered_boxes(const Vector<real_box> a_covered_boxes){
  m_covered_boxes = a_covered_boxes;
}

void fast_gshop::add_bbox(const real_box a_rbox){
  m_bboxes.push_back(a_rbox);
}

void fast_gshop::add_regular_box(const real_box a_rbox){
  m_regular_boxes.push_back(a_rbox);
}

void fast_gshop::add_covered_box(const real_box a_rbox){
  m_covered_boxes.push_back(a_rbox);
}

void fast_gshop::makeGrids(const ProblemDomain&      a_domain,
			   DisjointBoxLayout&        a_grids,
			   const int&                a_maxGridSize,
			   const int&                a_maxIrregGridSize){


    
  // 1. Find the resolution that corresponds to a_domain
  Real dx;
  bool found_dx;
  for (int lvl = 0 ; lvl < m_domains.size(); lvl++){
    if(a_domain == m_domains[lvl]){
      dx = m_dx[lvl];
      found_dx = true;
      break;
    }
  }

  // 1.1 Debug
  if(!found_dx){
    MayDay::Abort("fast_gshop::makeGrids - logic bust");
  }

  
  // 2. Do a domain split of a_domain and get all boxes
  Vector<Box> boxes;
  domainSplit(a_domain, boxes, a_maxGridSize, 1);

  // 3. Find out which boxes intersect with the bounding volume boxes
  Vector<Box> regular_boxes;
  Vector<Box> irregular_boxes;
  
  for (int ibox = 0; ibox < boxes.size(); ibox++){
    const Box cur_box = boxes[ibox];
    
    const real_box rbox = real_box(cur_box, m_origin, dx);
    if(is_regular_box(rbox) || is_covered_box(rbox)){
      regular_boxes.push_back(cur_box);
    }
    else{
      irregular_boxes.push_back(cur_box);
    }
  }

  // 4. Load balance the boxes that are outside & inside the bounding volumes independently
  Vector<int> regular_procs;
  Vector<int> irregular_procs;

  mortonOrdering(regular_boxes);
  mortonOrdering(irregular_boxes);

  LoadBalance(regular_procs, regular_boxes);
  LoadBalance(irregular_procs, irregular_boxes);

  Vector<int> procs;
  boxes.resize(0);
  
  procs.append(regular_procs);
  procs.append(irregular_procs);

  boxes.append(regular_boxes);
  boxes.append(irregular_boxes);

  a_grids.define(boxes, procs, a_domain);
  
#if 0 // Debug
  if(procID() == 0 && dx == m_dx[0]){
    for (int i = 0; i < irregular_boxes.size(); i++){
      std::cout << i << "\t" << irregular_procs[i] << std::endl;
    }
  }
#endif
}

bool fast_gshop::is_covered_box(const real_box& a_rbox) const {

  bool is_covered = false;

  for(int ibox = 0; ibox < m_covered_boxes.size(); ibox++){
    const real_box& covered_box = m_covered_boxes[ibox];

    if(covered_box.is_box_inside(a_rbox)){
      is_covered = true;
      break;
    }
  }

  return is_covered;
}

bool fast_gshop::is_regular_box(const real_box& a_rbox) const {

  bool is_regular = false;

  for(int ibox = 0; ibox < m_regular_boxes.size(); ibox++){
    const real_box& regular_box = m_regular_boxes[ibox];

    if(regular_box.is_box_inside(a_rbox)){
      is_regular = true;
      break;
    }
  }

  return is_regular;
}

bool fast_gshop::is_bounding_box(const real_box& a_rbox) const {

  bool is_bbox = false;

  for(int ibox = 0; ibox < m_bboxes.size(); ibox++){
    const real_box& bbox = m_bboxes[ibox];

    if(bbox.intersect(a_rbox)){
      is_bbox = true;
      break;
    }
  }

  return is_bbox;
}



GeometryService::InOut fast_gshop::InsideOutside(const Box&           a_region,
						 const ProblemDomain& a_domain,
						 const RealVect&      a_origin,
						 const Real&          a_dx) const {

  // Check if a_region falls outside of any of the bounding boxes. If it does, this is a regular box. If we
  // don't have any bounding boxes, we call GeometryShop parent function which does a point-wise iteration
  const real_box rbox(a_region, a_origin, a_dx);

  const bool regular  = is_regular_box(rbox);
  const bool covered  = is_covered_box(rbox);
  const bool bounding = is_bounding_box(rbox);

  if(regular && !bounding){
    return GeometryService::Regular;
  }
  else if(covered && !bounding){
    return GeometryService::Covered;
  }
  else{
    return GeometryShop::InsideOutside(a_region, a_domain, a_origin, a_dx);
  }
}
