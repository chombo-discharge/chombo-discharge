/*!
  @brief  ScanShop.cpp
  @brief  Implementation of ScanShop
  @author Robert Marskar
  @date   2019
*/

#include "ScanShop.H"

#include <BRMeshRefine.H>
#include <LoadBalance.H>
#include <EBLevelDataOps.H>

#include <ParmParse.H>


bool ScanShop::s_irregularBalance = true;
bool ScanShop::s_recursive        = true;
int ScanShop::s_grow              = 8;

ScanShop::ScanShop(const BaseIF&       a_localGeom,
		   const int           a_verbosity,
		   const Real          a_dx,
		   const RealVect      a_origin,
		   const ProblemDomain a_finestDomain,
		   const ProblemDomain a_scanLevel,
		   const Real          a_thrshdVoF)
  : GeometryShop(a_localGeom, a_verbosity, a_dx*RealVect::Unit, a_thrshdVoF) {


  m_baseif = &a_localGeom;
  m_hasScanLevel = false;

    // EBISLevel doesn't give resolution, origin, and problem domains through makeGrids, so we
  // need to construct these here, and then extract the proper resolution when we actually do makeGrids
  ScanShop::makeDomains(a_dx, a_origin, a_finestDomain, a_scanLevel);
}

ScanShop::~ScanShop(){

}

void ScanShop::makeDomains(const Real          a_dx,
			   const RealVect      a_origin,
			   const ProblemDomain a_finestDomain,
			   const ProblemDomain a_scanLevel){
  CH_TIME("ScanShop::makeDomains");

  m_origin = a_origin;
  
  m_dx.resize(0);
  m_domains.resize(0);
  
  m_domains.push_back(a_finestDomain);
  m_dx.push_back(a_dx);
  
  const int ref = 2;
  for (int lvl=0; ; lvl++){
    Real  dx             = m_dx[lvl];
    ProblemDomain domain = m_domains[lvl];

    if(a_scanLevel.domainBox() == domain.domainBox()){
      m_scanLevel = lvl;
    }
    
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

  // These will be built when they are needed. 
  m_grids.resize(m_domains.size());
  m_boxMaps.resize(m_domains.size());
  m_hasThisLevel.resize(m_domains.size(), 0);
}

void ScanShop::makeGrids(const ProblemDomain& a_domain,
			 DisjointBoxLayout&   a_grids,
			 const int&           a_maxGridSize,
			 const int&           a_maxIrregGridSize){


  // Build the scan level first
  if(!m_hasScanLevel){
    for (int lvl = m_domains.size()-1; lvl >= m_scanLevel; lvl--){
      ScanShop::buildCoarseLevel(lvl, a_maxGridSize); // Coarser levels built in the same way as the scan level
    }
    ScanShop::buildFinerLevels(m_scanLevel, a_maxGridSize);   // Traverse towards finer levels
  }

  // Development code. Break up a_domnain in a_maxGridSize chunks, load balance trivially and return the dbl
  Vector<Box> boxes;
  Vector<int> procs;
  
  domainSplit(a_domain, boxes, a_maxGridSize, 1);
  mortonOrdering(boxes);
  LoadBalance(procs, boxes);

  a_grids.define(boxes, procs, a_domain);
}

void ScanShop::buildCoarseLevel(const int a_level, const int a_maxGridSize){
  CH_TIME("ScanShop::buildCoarseLevel");

  // This function does the following:
  // 1. Break up the domain into chunks of max size a_maxGridSize and load balance them trivially
  // 2. Search through all the boxes and label them as regular/covered/irregular. Store the result in a LevelData map
  // 3. Gather cut-cell and regular/covered boxes separately and load balance them separately
  // 4. Create a new DBL with the newly load-balanced boxes; there should be approximately the same amount
  //    of cut-cell boxes for each rank
  // 5. Copy the map created over the initial DisjointBoxLayout onto the final grid
  
  // 1. 
  Vector<Box> boxes;
  Vector<int> procs;
  domainSplit(m_domains[a_level], boxes, a_maxGridSize, 1);
  LoadBalance(procs, boxes);

  // 2. 
  DisjointBoxLayout dbl(boxes, procs, m_domains[a_level]);
  LevelData<BoxType> map(dbl, 1, IntVect::Zero, BoxTypeFactory());
  Vector<Box> CutCellBoxes;
  Vector<Box> ReguCovBoxes;
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());
    
    Box grownBox = box;
    grownBox.grow(ScanShop::s_grow);
    
    bool foundNegative = false;
    bool foundPositive = false;
    for (BoxIterator bit(grownBox); bit.ok(); ++bit){
      const IntVect iv   = bit();
      const RealVect pos = m_origin + RealVect(iv)*m_dx[a_level];

      const Real value = m_baseif->value(pos);

      if(value < 0.0){
	foundNegative = true;
      }
      else if(value > 0.0){
	foundPositive = true;
      }
    }

    if(foundNegative && !foundPositive){ // Box must be covered
      map[dit()].setCovered();
      ReguCovBoxes.push_back(box);
    }
    else if(!foundNegative && foundPositive){ // Box is all regular
      map[dit()].setRegular();
      ReguCovBoxes.push_back(box);
    }
    else{ // Box must be a cut-cell box
      map[dit()].setCutCell();
      CutCellBoxes.push_back(box);
    }
  }

  // 3. Gather the uncut boxes and the cut boxes separately
  ScanShop::gatherBoxesParallel(CutCellBoxes);
  ScanShop::gatherBoxesParallel(ReguCovBoxes);

  Vector<int> CutCellProcs;
  Vector<int> ReguCovProcs;

  LoadBalance(CutCellProcs, CutCellBoxes);
  LoadBalance(ReguCovProcs, ReguCovBoxes);

  Vector<Box> allBoxes;
  Vector<int> allProcs;

  allBoxes.append(CutCellBoxes);
  allBoxes.append(ReguCovBoxes);

  allProcs.append(CutCellProcs);
  allProcs.append(ReguCovProcs);

  // 4. Define the grids on this level
  m_grids[a_level]   = DisjointBoxLayout(allBoxes, allProcs, m_domains[a_level]);
  m_boxMaps[a_level] = RefCountedPtr<LevelData<BoxType> >
    (new LevelData<BoxType>(m_grids[a_level], 1, IntVect::Zero, BoxTypeFactory()));

  // 5. Finally, copy the map
  map.copyTo(*m_boxMaps[a_level]);

  // We're done!
  m_hasThisLevel[a_level] = 1;
}

void ScanShop::buildFinerLevels(const int a_coarserLevel, const int a_maxGridSize){
  CH_TIME("ScanShop::buildFinerLevels");


  if(a_coarserLevel > 0){
    const int coarLvl = a_coarserLevel;
    const int fineLvl = coarLvl - 1;


    m_hasThisLevel[fineLvl] = 1;

    // This does the next level, we stop when level 0 has been built
    ScanShop::buildFinerLevels(fineLvl-1, a_maxGridSize);
  }
}

GeometryService::InOut ScanShop::InsideOutside(const Box&           a_region,
					       const ProblemDomain& a_domain,
					       const RealVect&      a_origin,
					       const Real&          a_dx,
					       const DataIndex&     a_dit) const{
  CH_TIME("ScanShop::InsideOutSide");

  // During development stage, return the standard stuff
  return GeometryService::InsideOutside(a_region, a_domain, a_origin, a_dx, a_dit);
}

void ScanShop::gatherBoxesParallel(Vector<Box>& a_boxes){

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
