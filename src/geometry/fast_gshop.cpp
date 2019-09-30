/*!
  @brief  fast_gshop.cpp
  @brief  Implementation of fast_gshop
  @author Robert Marskar
  @date   2019
*/

#include "fast_gshop.H"

#include <BRMeshRefine.H>
#include <LoadBalance.H>
#include <EBLevelDataOps.H>

#include <ParmParse.H>


bool fast_gshop::s_irregular_balance = true;
bool fast_gshop::s_recursive         = true;
int fast_gshop::s_grow               = 4;

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


  m_bounded_voxels.resize(0);
  m_regular_voxels.resize(0);
  m_covered_voxels.resize(0);
}

fast_gshop::~fast_gshop(){

}

void fast_gshop::set_bounded_voxels(const Vector<real_box> a_bounded_voxels){
  m_bounded_voxels = a_bounded_voxels;
}

void fast_gshop::set_regular_voxels(const Vector<real_box> a_regular_voxels){
  m_regular_voxels = a_regular_voxels;
}

void fast_gshop::set_covered_voxels(const Vector<real_box> a_covered_voxels){
  m_covered_voxels = a_covered_voxels;
}

void fast_gshop::add_bounded_voxel(const real_box a_rbox){
  m_bounded_voxels.push_back(a_rbox);
}

void fast_gshop::add_regular_voxel(const real_box a_rbox){
  m_regular_voxels.push_back(a_rbox);
}

void fast_gshop::add_covered_voxel(const real_box a_rbox){
  m_covered_voxels.push_back(a_rbox);
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


  // Find the resolution that corresponds to a_domain
  Real dx;
  bool found_dx;
  for (int lvl = 0 ; lvl < m_domains.size(); lvl++){
    if(a_domain == m_domains[lvl]){
      dx = m_dx[lvl];
      found_dx = true;
      break;
    }
  }

  if(!found_dx) MayDay::Abort("fast_gshop::makeGrids_recursive - logic bust");

  Vector<Box> regular_boxes;
  Vector<Box> irregular_boxes;
  Vector<Box> boxes;

  Vector<int> regular_procs;
  Vector<int> irregular_procs;
  Vector<int> procs;

  // Make boxes recursively and load balance them
  makeBoxes(regular_boxes, irregular_boxes, a_domain.domainBox(), a_domain, dx, a_maxGridSize);

  mortonOrdering(regular_boxes);
  mortonOrdering(irregular_boxes);

  LoadBalance(regular_procs, regular_boxes);
  LoadBalance(irregular_procs, irregular_boxes);

  // Define grids based on bounding voxels 
  procs.append(regular_procs);
  procs.append(irregular_procs);

  boxes.append(regular_boxes);
  boxes.append(irregular_boxes);

  a_grids.define(boxes, procs, a_domain);

  //  pout() << "fast_gshop::makeGrids_recursive - domain = " << a_domain << "\t num boxes = " << boxes.size() << endl;
  
  // 2. If we're load balancing with the irregular boxes, each grid now looks through his current boxes and reassign these
  //    by calling the GeometryShop::InsideOutside function for boxes that intersect with one of the bounded voxels
  if(s_irregular_balance){ 

    Vector<Box> reg_boxes; // Currently local to each rank
    Vector<Box> irr_boxes; // Currently local to each rank

    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit){
      const Box box = a_grids.get(dit());
      GeometryService::InOut inout = this->InsideOutside(box, a_domain, m_origin, dx);
      if(inout == GeometryService::Irregular){
	irr_boxes.push_back(box);
      }
      else{
	reg_boxes.push_back(box);
      }
    }
    
    // Gather the irregular and regular boxes
    gather_boxes_parallel(irr_boxes);
    gather_boxes_parallel(reg_boxes);

    // Load balance again
    Vector<int> irr_procs;
    Vector<int> reg_procs;

    LoadBalance(reg_procs, reg_boxes);
    LoadBalance(irr_procs, irr_boxes);

    boxes.resize(0);
    procs.resize(0);

    boxes.append(reg_boxes);
    boxes.append(irr_boxes);

    procs.append(reg_procs);
    procs.append(irr_procs);

    a_grids = DisjointBoxLayout();
    a_grids.define(boxes, procs, a_domain);
  }
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

    if((regular || covered) && !bounded && length){
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
}

bool fast_gshop::is_covered_box(const real_box& a_rbox) const {

  bool is_covered = false;

  for(int ibox = 0; ibox < m_covered_voxels.size(); ibox++){
    const real_box& covered_box = m_covered_voxels[ibox];

    if(covered_box.is_box_inside(a_rbox)){
      is_covered = true;
      break;
    }
  }

  return is_covered;
}

bool fast_gshop::is_regular_box(const real_box& a_rbox) const {

  bool is_regular = false;
  
  for(int ibox = 0; ibox < m_regular_voxels.size(); ibox++){
    const real_box& regular_box = m_regular_voxels[ibox];

    const bool inside_box = regular_box.is_box_inside(a_rbox);
    if(inside_box){
      is_regular = true;
      break;
    }
  }

  return is_regular;
}

bool fast_gshop::is_bounded_box(const real_box& a_rbox) const {

  bool is_bbox = false;

  for(int ibox = 0; ibox < m_bounded_voxels.size(); ibox++){
    const real_box& bbox = m_bounded_voxels[ibox];

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

  if(regular && !bounded){
    return GeometryService::Regular;
  }
  else if(covered && !bounded){
    return GeometryService::Covered;
  }
  else{
    return GeometryShop::InsideOutside(a_region, a_domain, a_origin, a_dx);
  }
  // This code defaults to GeometryShops insideoutside stuff
  // return GeometryShop::InsideOutside(a_region, a_domain, a_origin, a_dx);
}


void fast_gshop::gather_boxes_parallel(Vector<Box>& a_boxes){

  // 1. Linearize local boxes
  int send_size    = 2*CH_SPACEDIM;                 // Message size for one box
  int send_count   = a_boxes.size()*send_size;      // Number of elements sent from this rank
  int* send_buffer = new int[send_count*send_size]; // Send buffer for this rank
  int* send_buf2   = send_buffer;                   // Backup address. Going to monkey with pointer increments on send buffer

  // Linearize a_boxes onto send_buffer
  for (int i = 0; i < a_boxes.size(); i++, send_buffer+=send_size){
    const Box& b = a_boxes[i];
    D_TERM6(send_buffer[0] =b.smallEnd(0); send_buffer[1] =b.bigEnd(0);,
	    send_buffer[2] =b.smallEnd(1); send_buffer[3] =b.bigEnd(1);,
	    send_buffer[4] =b.smallEnd(2); send_buffer[5] =b.bigEnd(2);,
	    send_buffer[6] =b.smallEnd(3); send_buffer[7] =b.bigEnd(3);,
	    send_buffer[8] =b.smallEnd(4); send_buffer[9] =b.bigEnd(4);,
	    send_buffer[10]=b.smallEnd(5); send_buffer[11]=b.bigEnd(5););
  }
  send_buffer = send_buf2; // Revert point to start of array


  // 2. Get the number of elements sent from each rank
  int* send_counts = new int[numProc()];
  MPI_Allgather(&send_count, 1, MPI_INT, send_counts, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0] = 0;
  for (int i = 0; i < numProc()-1; i++){
    offsets[i+1] = offsets[i] + send_counts[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++){
    total_count += send_counts[i];
  }
  int* recv_buffer = new int[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer, send_count, MPI_INT, recv_buffer, send_counts, offsets, MPI_INT, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_boxes.resize(0);
  int* recv_buf2 = recv_buffer; // Going to monkey with pointer increments again
  for (int i = 0; i < total_count/send_size; i++, recv_buffer+=send_size){
    IntVect lo, hi;
    D_TERM6(lo[0] = recv_buffer[0];  hi[0] = recv_buffer[1];,
	    lo[1] = recv_buffer[2];  hi[1] = recv_buffer[3];,
	    lo[2] = recv_buffer[4];  hi[2] = recv_buffer[5];,
	    lo[3] = recv_buffer[6];  hi[3] = recv_buffer[7];,
	    lo[4] = recv_buffer[8];  hi[4] = recv_buffer[9];,
	    lo[5] = recv_buffer[10]; hi[5] = recv_buffer[11];);

    a_boxes.push_back(Box(lo, hi));
  }
  recv_buffer = recv_buf2;

  delete recv_buffer;
  delete send_buffer;
}
