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

#define DEBUG 0

#include "CD_NamespaceHeader.H"

bool ScanShop::s_irregularBalance = true;
bool ScanShop::s_recursive        = true;
int ScanShop::s_grow              = 4;

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

#if DEBUG
    for (int lvl = 0; lvl < m_domains.size(); lvl++){
      ScanShop::printNumBoxesLevel(lvl);
    }
#endif

    m_hasScanLevel = true;
  }

  // Find the level corresponding to a_domain
  int whichLevel;
  for (int lvl = 0; lvl < m_domains.size(); lvl++){
    if(m_domains[lvl].domainBox() == a_domain.domainBox()){
      whichLevel = lvl;
      break;
    }
  }

  if(m_hasThisLevel[whichLevel] != 0){
    a_grids = m_grids[whichLevel];
  }
  else{

    // Development code. Break up a_domnain in a_maxGridSize chunks, load balance trivially and return the dbl
    Vector<Box> boxes;
    Vector<int> procs;
  
    domainSplit(a_domain, boxes, a_maxGridSize, 1);
    mortonOrdering(boxes);
    LoadBalance(procs, boxes);

    a_grids.define(boxes, procs, a_domain);
  }
}

bool ScanShop::isRegular(const Box a_box, const RealVect a_origin, const Real a_dx){
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    const RealVect a_point = a_origin + a_dx*(0.5*RealVect::Unit + RealVect(iv));
    if(m_baseif->value(a_point) > -0.5*a_dx){
      return false;
    }
  }

  return true;
}

bool ScanShop::isCovered(const Box a_box, const RealVect a_origin, const Real a_dx){
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    const RealVect a_point = a_origin + a_dx*(0.5*RealVect::Unit + RealVect(iv));
    if(m_baseif->value(a_point) < 0.5*a_dx) return false;
  }

  return true;
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
    grownBox &= m_domains[a_level];

    const bool isRegular = ScanShop::isRegular(grownBox, m_origin, m_dx[a_level]);
    const bool isCovered = ScanShop::isCovered(grownBox, m_origin, m_dx[a_level]);

    if(isRegular && !isCovered){
      map[dit()].setRegular();
      ReguCovBoxes.push_back(box);
    }
    else if(isCovered && !isRegular){
      map[dit()].setCovered();
      ReguCovBoxes.push_back(box);
    }
    else if(!isRegular && !isCovered){
      map[dit()].setCutCell();
      CutCellBoxes.push_back(box);
    }
    else{
      MayDay::Abort("ScanShop::buildCoarseLevel - logic bust");
    }
  }

#if DEBUG // Debug first map
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    if(map[dit()].m_which == -1) MayDay::Abort("ScanShop::buildCoarseLevel - initial map shouldn't get -1");
  }
#endif

  // 3. Gather the uncut boxes and the cut boxes separately
  ScanShop::gatherBoxesParallel(CutCellBoxes);
  ScanShop::gatherBoxesParallel(ReguCovBoxes);

  mortonOrdering(CutCellBoxes);
  mortonOrdering(ReguCovBoxes);

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
  m_grids[a_level] = DisjointBoxLayout(allBoxes, allProcs, m_domains[a_level]);
  m_boxMaps[a_level] = RefCountedPtr<LevelData<BoxType> >
    (new LevelData<BoxType>(m_grids[a_level], 1, IntVect::Zero, BoxTypeFactory()));

  // 5. Finally, copy the map
  map.copyTo(*m_boxMaps[a_level]);

#if DEBUG // Debug first map
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    if((*m_boxMaps[a_level])[dit()].m_which == -1) MayDay::Abort("ScanShop::buildCoarseLevel - final map shouldn't get -1");
  }
#endif

  // We're done!
  m_hasThisLevel[a_level] = 1;
}

void ScanShop::buildFinerLevels(const int a_coarserLevel, const int a_maxGridSize){
  CH_TIME("ScanShop::buildFinerLevels");


  if(a_coarserLevel > 0){

    const int coarLvl = a_coarserLevel;
    const int fineLvl = coarLvl - 1;

    // Coar stuff
    const DisjointBoxLayout& dblCoar  = m_grids[coarLvl];
    const LevelData<BoxType>& mapCoar = *m_boxMaps[coarLvl];


    Vector<Box> CutCellBoxes;
    Vector<Box> ReguCovBoxes;

    Vector<int> CutCellProcs;
    Vector<int> ReguCovProcs;

    // Refined coarse stuff.
    DisjointBoxLayout dblFineCoar;
    refine(dblFineCoar, dblCoar, 2);
    LevelData<BoxType> mapFineCoar(dblFineCoar, 1, IntVect::Zero, BoxTypeFactory());

    // Break up cut cell boxes
    for (DataIterator dit = dblFineCoar.dataIterator(); dit.ok(); ++dit){
      mapFineCoar[dit()] = mapCoar[dit()];
      const Box box = dblFineCoar.get(dit());
      if(mapFineCoar[dit()].isCutCell()){
	Vector<Box> boxes;
	domainSplit(box, boxes, a_maxGridSize, 1);
	CutCellBoxes.append(boxes);
      }
      else{
	ReguCovBoxes.push_back(box);
      }
    }

    // Make a DBL out of the cut-cell boxes and again check if they actually contain cut cells
    gatherBoxesParallel(CutCellBoxes);
    mortonOrdering(CutCellBoxes);
    LoadBalance(CutCellProcs, CutCellBoxes);
    DisjointBoxLayout  CutCellDBL(CutCellBoxes, CutCellProcs, m_domains[fineLvl]);
    LevelData<BoxType> CutCellMap(CutCellDBL, 1, IntVect::Zero, BoxTypeFactory());

    // Redo the cut cell boxes
    CutCellBoxes.resize(0);
    for (DataIterator dit = CutCellDBL.dataIterator(); dit.ok(); ++dit){
      const Box box = CutCellDBL.get(dit());
      
      Box grownBox = box;
      grownBox.grow(ScanShop::s_grow);
      grownBox &= m_domains[fineLvl];

      const bool isRegular = ScanShop::isRegular(grownBox, m_origin, m_dx[fineLvl]);
      const bool isCovered = ScanShop::isCovered(grownBox, m_origin, m_dx[fineLvl]);

      if(isRegular){
	CutCellMap[dit()].setRegular();
	ReguCovBoxes.push_back(box);
      }
      else if(isCovered){
	CutCellMap[dit()].setCovered();
	ReguCovBoxes.push_back(box);
      }
      else if(!isRegular && !isCovered){
	CutCellMap[dit()].setCutCell();
	CutCellBoxes.push_back(box);
      }
      else{
	MayDay::Abort("ScanShop::buildCoarseLevel - logic bust");
      }
    }
    
    // Gather boxes with cut cells and the ones don't contain cut cells
    ScanShop::gatherBoxesParallel(CutCellBoxes);
    ScanShop::gatherBoxesParallel(ReguCovBoxes);

    mortonOrdering(CutCellBoxes);
    mortonOrdering(ReguCovBoxes);
    
    LoadBalance(CutCellProcs, CutCellBoxes);
    LoadBalance(ReguCovProcs, ReguCovBoxes);
    
    Vector<Box> allBoxes;
    Vector<int> allProcs;
    
    allBoxes.append(CutCellBoxes);
    allBoxes.append(ReguCovBoxes);
    
    allProcs.append(CutCellProcs);
    allProcs.append(ReguCovProcs);


    // Define the grids and copy the maps over. First copy the refine coarse map, then update with the
    // cut cell map
    m_grids[fineLvl] = DisjointBoxLayout(allBoxes, allProcs, m_domains[fineLvl]);
    m_boxMaps[fineLvl] = RefCountedPtr<LevelData<BoxType> >
      (new LevelData<BoxType>(m_grids[fineLvl], 1, IntVect::Zero, BoxTypeFactory()));

    mapFineCoar.copyTo(*m_boxMaps[fineLvl]);
    CutCellMap.copyTo(*m_boxMaps[fineLvl]);

#if DEBUG // Dummy check the box maps
    for (DataIterator dit = CutCellDBL.dataIterator(); dit.ok(); ++dit){
      if(CutCellMap[dit()].m_which == -1) MayDay::Abort("ScanShop::buildFinerLevels - CutCellMap shouldn't get -1");
    }
    for (DataIterator dit = dblFineCoar.dataIterator(); dit.ok(); ++dit){
      if(mapFineCoar[dit()].m_which == -1) MayDay::Abort("ScanShop::buildFinerLevels - mapFineCoar shouldn't get -1");
    }
    for (DataIterator dit = m_grids[fineLvl].dataIterator(); dit.ok(); ++dit){
      if((*m_boxMaps[fineLvl])[dit()].m_which == -1){ MayDay::Abort("ScanShop::buildFinerLevels - final map shouldn't get -1");}
    }
#endif
    
    m_hasThisLevel[fineLvl] = 1;

    // This does the next level, we stop when level 0 has been built
    ScanShop::buildFinerLevels(fineLvl, a_maxGridSize);
  }
}

GeometryService::InOut ScanShop::InsideOutside(const Box&           a_region,
					       const ProblemDomain& a_domain,
					       const RealVect&      a_origin,
					       const Real&          a_dx,
					       const DataIndex&     a_dit) const{
  CH_TIME("ScanShop::InsideOutSide");

  // Find the level corresponding to a_domain
  int whichLevel = -1;
  bool foundLevel = false;
  for (int lvl = 0; lvl < m_domains.size(); lvl++){
    if(m_domains[lvl].domainBox() == a_domain.domainBox()){
      whichLevel = lvl;
      foundLevel = true;
      break;
    }
  }

  // A strang but true thing. This function is used in EBISLevel::simplifyGraphFromGeo and that function can send in a_domain
  // and a_dx on different levels....
  ProblemDomain domain;
  if(a_dx < m_dx[whichLevel] && whichLevel > 0){
    domain = m_domains[whichLevel-1];

    MayDay::Abort("ScanShop::InsideOutside - shouldn't happen...");
  }

  GeometryService::InOut ret;

  if(foundLevel && m_hasThisLevel[whichLevel]){
    const LevelData<BoxType>& map = (*m_boxMaps[whichLevel]);
    const BoxType& boxType        = map[a_dit];
    
    if(boxType.isRegular()){
      ret = GeometryService::Regular;
    }
    else if(boxType.isCovered()){
      ret = GeometryService::Covered;
    }
    else{
      ret= GeometryService::Irregular;
      //      return GeometryService::InsideOutside(a_region, domain, a_origin, a_dx, a_dit);
    }
  }
  else{
    ret = GeometryService::InsideOutside(a_region, domain, a_origin, a_dx, a_dit);
  }

  return ret;
}

void ScanShop::printNumBoxesLevel(const int a_level){
  CH_TIME("ScanShop::printNumBoxesLevel");

  const DisjointBoxLayout&     dbl = m_grids[a_level];
  const LevelData<BoxType>& boxMap = *m_boxMaps[a_level];

  int myNumRegular = 0;
  int myNumCovered = 0;
  int myNumCutCell = 0;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const BoxType& bType = boxMap[dit()];

    if(bType.isRegular()){
      myNumRegular++;
    }
    else if(bType.isCovered()){
      myNumCovered++;
    }
    else if(bType.isCutCell()){
      myNumCutCell++;
    }
  }

  const int numRegular = EBLevelDataOps::parallelSum(myNumRegular);
  const int numCovered = EBLevelDataOps::parallelSum(myNumCovered);
  const int numCutCell = EBLevelDataOps::parallelSum(myNumCutCell);

  if(procID() == 0){
    std::cout << "ScanShop::printNumBoxesLevel on domain = " << m_domains[a_level] << "\n"
	      << "\t Regular boxes = " << numRegular << "\n"
	      << "\t Covered boxes = " << numCovered << "\n"
	      << "\t CutCell boxes = " << numCutCell << "\n"
	      << std::endl;
  }
    
}

void ScanShop::gatherBoxesParallel(Vector<Box>& a_boxes){

  // 1. Linearize local boxes
  int send_size    = 2*CH_SPACEDIM;                 // Message size for one box
  int send_count   = a_boxes.size()*send_size;      // Number of elements sent from this rank
#if 0 // Why did I multiply by send_size TWICE...?
  int* send_buffer = new int[send_count*send_size]; // Send buffer for this rank. 
#else
  int* send_buffer = new int[send_count]; // Send buffer for this rank
#endif
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
#include "CD_NamespaceFooter.H"
