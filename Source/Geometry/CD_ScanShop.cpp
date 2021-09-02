/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @brief  CD_ScanShop.cpp
  @brief  Implementation of CD_ScanShop.H
  @author Robert Marskar
*/

// Chombo includes
#include <BRMeshRefine.H>
#include <LoadBalance.H>
#include <EBLevelDataOps.H>
#include <ParmParse.H>

// Our includes
#include <CD_ScanShop.H>
#include <CD_NamespaceHeader.H>

#define DEBUG_ScanShop 0

ScanShop::ScanShop(const BaseIF&       a_localGeom,
		   const int           a_verbosity,
		   const Real          a_dx,
		   const RealVect      a_probLo,
		   const ProblemDomain a_finestDomain,
		   const ProblemDomain a_scanLevel,
		   const int           a_maxGhostEB,
		   const Real          a_thrshdVoF)
  : GeometryShop(a_localGeom, a_verbosity, a_dx*RealVect::Unit, a_thrshdVoF) {

  CH_TIME("ScanShop::ScanShop(BaseIF, int, Real, RealVect, ProblemDomain, ProblemDomain, int, Real)");

  m_baseIF       = &a_localGeom;
  m_hasScanLevel = false;
  m_maxGhostEB   = a_maxGhostEB;

  // EBISLevel doesn't give resolution, origin, and problem domains through makeGrids, so we
  // need to construct these here, and then extract the proper resolution when we actually call makeGrids
  ScanShop::makeDomains(a_dx, a_probLo, a_finestDomain, a_scanLevel);
}

ScanShop::~ScanShop(){

}

void ScanShop::makeDomains(const Real          a_dx,
			   const RealVect      a_probLo,
			   const ProblemDomain a_finestDomain,
			   const ProblemDomain a_scanLevel){
  CH_TIME("ScanShop::makeDomains(Real, RealVect, ProblemDomain, ProblemDomain)");

  m_probLo = a_probLo;
  
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
  CH_TIME("ScanShop::makeDomains(ProblemDomain, DisjointBoxLayout, int, int)");  

  // Build the scan level first
  if(!m_hasScanLevel){
    for (int lvl = m_domains.size()-1; lvl >= m_scanLevel; lvl--){
      ScanShop::buildCoarseLevel(lvl, a_maxGridSize); // Coarser levels built in the same way as the scan level
    }
    ScanShop::buildFinerLevels(m_scanLevel, a_maxGridSize);   // Traverse towards finer levels

#if DEBUG_ScanShop
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

bool ScanShop::isRegular(const Box a_box, const RealVect a_probLo, const Real a_dx) const {
  CH_TIME("ScanShop::isRegular(Box, RealVect, Real)");

  bool ret = true;
  
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    const RealVect a_point = a_probLo + a_dx*(0.5*RealVect::Unit + RealVect(iv));
    if(m_baseIF->value(a_point) > -0.5*a_dx){
      ret = false;
      break;
    }
  }

  return ret;
}

bool ScanShop::isCovered(const Box a_box, const RealVect a_probLo, const Real a_dx) const {
  CH_TIME("ScanShop::isCovered(Box, RealVect, Real)");

  bool ret = true;
  
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    const RealVect a_point = a_probLo + a_dx*(0.5*RealVect::Unit + RealVect(iv));
    if(m_baseIF->value(a_point) < 0.5*a_dx) {
      ret = false;
      break;
    }
  }

  return ret;
}

void ScanShop::buildCoarseLevel(const int a_level, const int a_maxGridSize){
  CH_TIME("ScanShop::buildCoarseLevel(int, int)");

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
    grownBox.grow(m_maxGhostEB);
    grownBox &= m_domains[a_level];

    const bool isRegular = ScanShop::isRegular(grownBox, m_probLo, m_dx[a_level]);
    const bool isCovered = ScanShop::isCovered(grownBox, m_probLo, m_dx[a_level]);

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

#if DEBUG_ScanShop // Debug first map
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

#if DEBUG_ScanShop // Debug first map
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    if((*m_boxMaps[a_level])[dit()].m_which == -1) MayDay::Abort("ScanShop::buildCoarseLevel - final map shouldn't get -1");
  }
#endif

  // We're done!
  m_hasThisLevel[a_level] = 1;
}

void ScanShop::buildFinerLevels(const int a_coarserLevel, const int a_maxGridSize){
  CH_TIME("ScanShop::buildFinerLevels(int, int)");

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
      grownBox.grow(m_maxGhostEB);
      grownBox &= m_domains[fineLvl];

      const bool isRegular = ScanShop::isRegular(grownBox, m_probLo, m_dx[fineLvl]);
      const bool isCovered = ScanShop::isCovered(grownBox, m_probLo, m_dx[fineLvl]);

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
	MayDay::Abort("ScanShop::buildFinerLevels - logic bust!");
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

#if DEBUG_ScanShop // Dummy check the box maps
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
					       const RealVect&      a_probLo,
					       const Real&          a_dx,
					       const DataIndex&     a_dit) const{
  CH_TIME("ScanShop::InsideOutSide(Box, ProblemDomain, RealVect, Real, DataIndex)");

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
      //      return GeometryService::InsideOutside(a_region, domain, a_probLo, a_dx, a_dit);
    }
  }
  else{
    ret = GeometryService::InsideOutside(a_region, domain, a_probLo, a_dx, a_dit);
  }

  return ret;
}

void ScanShop::printNumBoxesLevel(const int a_level) const {
  CH_TIME("ScanShop::printNumBoxesLevel(int)");

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

void ScanShop::gatherBoxesParallel(Vector<Box>& a_boxes) const {
  CH_TIME("ScanShop::gatherBoxesParallel(Vector<Box>)");
  
#ifdef CH_MPI
  // 1. Linearize local boxes
  int sendSize     = 2*CH_SPACEDIM;                 // Message size for one box
  int sendCount    = a_boxes.size()*sendSize;      // Number of elements sent from this rank
  int* sendBuffer  = new int[sendCount]; // Send buffer for this rank
  int* sendBuffer2 = sendBuffer;                   // Backup address. Going to monkey with pointer increments on send buffer

  // Linearize a_boxes onto sendBuffer
  for (int i = 0; i < a_boxes.size(); i++, sendBuffer+=sendSize){
    const Box& b = a_boxes[i];
    D_TERM6(sendBuffer[0] =b.smallEnd(0); sendBuffer[1] =b.bigEnd(0);,
	    sendBuffer[2] =b.smallEnd(1); sendBuffer[3] =b.bigEnd(1);,
	    sendBuffer[4] =b.smallEnd(2); sendBuffer[5] =b.bigEnd(2);,
	    sendBuffer[6] =b.smallEnd(3); sendBuffer[7] =b.bigEnd(3);,
	    sendBuffer[8] =b.smallEnd(4); sendBuffer[9] =b.bigEnd(4);,
	    sendBuffer[10]=b.smallEnd(5); sendBuffer[11]=b.bigEnd(5););
  }
  sendBuffer = sendBuffer2; // Revert point to start of array


  // 2. Get the number of elements sent from each rank
  int* sendCounts = new int[numProc()];
  MPI_Allgather(&sendCount, 1, MPI_INT, sendCounts, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0] = 0;
  for (int i = 0; i < numProc()-1; i++){
    offsets[i+1] = offsets[i] + sendCounts[i];
  }

  // 4. Allocate storage for total buffer size
  int totalCount = 0;
  for (int i = 0; i < numProc(); i++){
    totalCount += sendCounts[i];
  }
  int* receiveBuffer = new int[totalCount];

  // 5. MPI send
  MPI_Allgatherv(sendBuffer, sendCount, MPI_INT, receiveBuffer, sendCounts, offsets, MPI_INT, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_boxes.resize(0);
  int* recv_buf2 = receiveBuffer; // Going to monkey with pointer increments again
  for (int i = 0; i < totalCount/sendSize; i++, receiveBuffer+=sendSize){
    IntVect lo, hi;
    D_TERM6(lo[0] = receiveBuffer[0];  hi[0] = receiveBuffer[1];,
	    lo[1] = receiveBuffer[2];  hi[1] = receiveBuffer[3];,
	    lo[2] = receiveBuffer[4];  hi[2] = receiveBuffer[5];,
	    lo[3] = receiveBuffer[6];  hi[3] = receiveBuffer[7];,
	    lo[4] = receiveBuffer[8];  hi[4] = receiveBuffer[9];,
	    lo[5] = receiveBuffer[10]; hi[5] = receiveBuffer[11];);

    a_boxes.push_back(Box(lo, hi));
  }
  receiveBuffer = recv_buf2;

  delete[] sendCounts;
  delete[] offsets;
  delete[] receiveBuffer;
  delete[] sendBuffer;
#endif
}

#include <CD_NamespaceFooter.H>
