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

bool fast_gshop::s_recursive = true;
int fast_gshop::s_grow = 4;

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
      dx *= ref;
      
      m_dx.push_back(dx);
      m_domains.push_back(domain);
    }
    else{
      break;
    }
  }


  m_bounded_boxes.resize(0);
  m_regular_boxes.resize(0);
  m_covered_boxes.resize(0);
  
#if 0 // Debug code
  ParmParse pp("devel");
  int use_bbox = 0;
  pp.get("use_bbox", use_bbox);
  if(use_bbox != 0){



    // Evertying to regular
    m_regular_boxes.push_back(real_box(-RealVect::Unit, RealVect::Unit));

    // Except this boxc

    // Everything below needle
#if 0
    Real radius = 1.E-3;
    RealVect lo = RealVect::Zero;
    RealVect hi = RealVect(0, 0, 1);

    lo = lo - radius*1.001*RealVect::Unit;
    hi = hi + radius*1.001*RealVect::Unit;

    //    m_bounded_boxes.push_back(real_box(lo, hi));
#endif
  }

#endif
}

fast_gshop::~fast_gshop(){

}

void fast_gshop::set_bounded_boxes(const Vector<real_box> a_bounded_boxes){
  m_bounded_boxes = a_bounded_boxes;
}

void fast_gshop::set_regular_boxes(const Vector<real_box> a_regular_boxes){
  m_regular_boxes = a_regular_boxes;

  
#if 0
  std::cout << m_regular_boxes.size() << std::endl;
#endif
}

void fast_gshop::set_covered_boxes(const Vector<real_box> a_covered_boxes){

  m_covered_boxes = a_covered_boxes;
#if 0
    std::cout << m_covered_boxes.size() << std::endl;
#endif
}

void fast_gshop::add_bounded_box(const real_box a_rbox){
  m_bounded_boxes.push_back(a_rbox);
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

  if(fast_gshop::s_recursive){
    makeGrids_recursive(a_domain, a_grids, a_maxGridSize, a_maxIrregGridSize);
  }
  else{
    makeGrids_domainSplit(a_domain, a_grids, a_maxGridSize, a_maxIrregGridSize);
  }
}

void fast_gshop::makeGrids_recursive(const ProblemDomain&      a_domain,
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
    MayDay::Abort("fast_gshop::makeGrids_recursive - logic bust");
  }

  Vector<Box> regular_boxes;
  Vector<Box> irregular_boxes;

  Vector<int> regular_procs;
  Vector<int> irregular_procs;

  // Make boxes recursively
#if 0
  pout() << "starting make boxes" << endl;
  pout() << "dx = " << dx << endl;
  pout() << "domain = " << a_domain.domainBox() << endl;
#endif
  makeBoxes(regular_boxes, irregular_boxes, a_domain.domainBox(), a_domain, dx, a_maxGridSize);
#if 0
  pout() << "end make boxes" << endl;
  pout() << "num boxes = " << regular_boxes.size() + irregular_boxes.size() << endl;
#endif


  mortonOrdering(regular_boxes);
  mortonOrdering(irregular_boxes);

  LoadBalance(regular_procs, regular_boxes);
  LoadBalance(irregular_procs, irregular_boxes);

  Vector<int> procs;
  Vector<Box> boxes;
  
  procs.append(regular_procs);
  procs.append(irregular_procs);

  boxes.append(regular_boxes);
  boxes.append(irregular_boxes);

  pout() << "fast_gshop::makeGrids_recursive - domain = " << a_domain << "\t num boxes = " << boxes.size() << endl;

  a_grids.define(boxes, procs, a_domain);

}

void fast_gshop::makeBoxes(Vector<Box>&         a_reg_boxes,
			   Vector<Box>&         a_irreg_boxes,
			   const Box&           a_region,
			   const ProblemDomain& a_domain,
			   const Real           a_dx, 
			   const int            a_maxGridSize){

  int longdir;
  int length = a_region.longside(longdir); // Not necessarily multiple of two, but always multiple of a_maxGridSize

  // std::cout << length << "\t" << a_maxGridSize << std::endl;
#if 0
  if(length % a_maxGridSize != 0){
    MayDay::Abort("fast_gshop::makeBoxes - shouldn't happen!");
  }
#endif

  // Split domain along the longest direction. We happen to know that a_region[longdir] is a integer multiple of maxgridsize
  if(length > a_maxGridSize){

    Box box = a_region;
    box.grow(s_grow);
    const real_box rbox(box, m_origin, a_dx);
    const bool regular  = is_regular_box(rbox);
    const bool covered  = is_covered_box(rbox);
    const bool bounded  = is_bounded_box(rbox);

    if((regular || covered) && !bounded && length < 1024){
#if 0
      pout() << "pushing regular box = " << rbox.get_lo() << "\t" << rbox.get_hi() << endl;
      pout() << "dx = " << a_dx << endl;
      pout() << "domain box = " << a_region << endl;
      pout() << "regular = " << regular << endl;
      pout() << "covered = " << covered << endl;
      pout() << "bounded = " << bounded << endl;
#endif
      a_reg_boxes.push_back(a_region);
    }
    else{
      //      pout() << "splitting boxes" << endl;
      const int npatches = length/a_maxGridSize;
      const int n_left   = npatches/2;
      const int n_right  = npatches - n_left;


      // Divide the input box into two halves, left and right with n_left and n_right patches each
      Box lo(a_region), hi;
      hi = lo.chop(longdir, a_region.smallEnd(longdir) + n_left*a_maxGridSize);
      
      makeBoxes(a_reg_boxes, a_irreg_boxes, lo, a_domain, a_dx, a_maxGridSize);
      makeBoxes(a_reg_boxes, a_irreg_boxes, hi, a_domain, a_dx, a_maxGridSize);
    }
  }
  else{
    // Just put them somewhere, not important right now
    //    pout() << "pushing irregular box" << endl;
    a_irreg_boxes.push_back(a_region);
  }
}


void fast_gshop::makeGrids_domainSplit(const ProblemDomain&      a_domain,
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
    MayDay::Abort("fast_gshop::makeGrids_split - logic bust");
  }

  
  // 2. Do a domain split of a_domain and get all boxes
  Vector<Box> boxes;
  domainSplit(a_domain, boxes, a_maxGridSize, 1);

  // 3. Find out which boxes intersect with the bounding volume boxes
  Vector<Box> regular_boxes;
  Vector<Box> irregular_boxes;
  
  for (int ibox = 0; ibox < boxes.size(); ibox++){
    const Box cur_box = boxes[ibox];
    Box grown_box = cur_box;
    grown_box.grow(s_grow);
    
    
    const real_box rbox = real_box(grown_box, m_origin, dx);

    const bool regular = is_regular_box(rbox);
    const bool covered = is_covered_box(rbox);
    const bool bounded = is_bounded_box(rbox);
    
    if((regular || covered) && !bounded){
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

    const bool inside_box = regular_box.is_box_inside(a_rbox);
    if(inside_box){
      is_regular = true;
      break;
    }
  }
  
  //  std::cout << is_regular << "\t" << m_regular_boxes[0].get_lo() << "\t" << m_regular_boxes[0].get_hi() << std::endl;

  return is_regular;
}

bool fast_gshop::is_bounded_box(const real_box& a_rbox) const {

  bool is_bbox = false;

  for(int ibox = 0; ibox < m_bounded_boxes.size(); ibox++){
    const real_box& bbox = m_bounded_boxes[ibox];

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
  Box box = a_region;
  box.grow(s_grow);
  const real_box rbox(box, a_origin, a_dx);

  const bool regular = is_regular_box(rbox);
  const bool covered = is_covered_box(rbox);
  const bool bounded = is_bounded_box(rbox);

#if 1 // Original code
  if(regular && !bounded){
    return GeometryService::Regular;
  }
  else if(covered && !bounded){
    return GeometryService::Covered;
  }
  else{
    return GeometryShop::InsideOutside(a_region, a_domain, a_origin, a_dx);
  }
#else // This code defaults to GeometryShops insideoutside stuff
  return GeometryShop::InsideOutside(a_region, a_domain, a_origin, a_dx);
#endif
}
