/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @brief  CD_ScanShop.cpp
  @brief  Implementation of CD_ScanShop.H
  @author Robert Marskar
*/

// Std includes
#include <chrono>

// Chombo includes
#include <BRMeshRefine.H>
#include <LoadBalance.H>
#include <EBLevelDataOps.H>
#include <ParmParse.H>

// Our includes
#include <CD_ScanShop.H>
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

#define ScanShop_Debug 0

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
  m_profile      = false;
  m_loadBalance  = false;
  m_maxGhostEB   = 1;//a_maxGhostEB;
  m_fileName     = "ScanShopReport.dat";
  m_boxSorting   = BoxSorting::Morton;

  // EBISLevel doesn't give resolution, origin, and problem domains through makeGrids, so we
  // need to construct these here, and then extract the proper resolution when we actually call makeGrids
  ScanShop::makeDomains(a_dx, a_probLo, a_finestDomain, a_scanLevel);

  // These are settings that I want to hide from the user.
  ParmParse pp("ScanShop");

  std::string str;
  pp.query("profile",      m_profile);
  pp.query("box_sorting",  str);

  if(str == "none"){
    m_boxSorting = BoxSorting::None;
  }
  else if(str == "std"){
    m_boxSorting = BoxSorting::Std;
  }
  else if(str == "morton"){
    m_boxSorting = BoxSorting::Morton;
  }
  else if(str == "shuffle"){
    m_boxSorting = BoxSorting::Shuffle;
  }

  m_timer = Timer("ScanShop");
}

ScanShop::~ScanShop(){
  CH_TIME("ScanShop::~ScanShop()");

  if(m_profile){
    m_timer.writeReportToFile(m_fileName);
    m_timer.eventReport(pout(), false);
  }
}

void ScanShop::setProfileFileName(const std::string a_fileName){
  m_fileName = a_fileName;

  m_timer = Timer(m_fileName);
}

void ScanShop::makeDomains(const Real          a_dx,
			   const RealVect      a_probLo,
			   const ProblemDomain a_finestDomain,
			   const ProblemDomain a_scanLevel){
  CH_TIME("ScanShop::makeDomains(Real, RealVect, ProblemDomain, ProblemDomain)");

  m_probLo = a_probLo;
  
  m_dx.     resize(0);
  m_domains.resize(0);

  m_dx.     push_back(a_dx);  
  m_domains.push_back(a_finestDomain);

  
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
      
      m_dx.     push_back(dx);
      m_domains.push_back(domain);
    }
    else{
      break;
    }
  }

  // These will be built when they are needed. 
  m_grids.       resize(m_domains.size());
  m_boxMap.      resize(m_domains.size());  
  m_hasThisLevel.resize(m_domains.size(), 0);
}

void ScanShop::makeGrids(const ProblemDomain& a_domain,
			 DisjointBoxLayout&   a_grids,
			 const int&           a_maxGridSize,
			 const int&           a_maxIrregGridSize){
  CH_TIME("ScanShop::makeGrids(ProblemDomain, DisjointBoxLayout, int, int)");  

  m_timer.startEvent("Make grids");

  // Build the scan level first
  if(!m_hasScanLevel){
    for (int lvl = m_domains.size()-1; lvl >= m_scanLevel; lvl--){
      ScanShop::buildCoarseLevel(lvl, a_maxGridSize); // Coarser levels built in the same way as the scan level
    }
    ScanShop::buildFinerLevels(m_scanLevel, a_maxGridSize);   // Traverse towards finer levels

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
    
    MayDay::Warning("ScanShop::makeGrids -- decomposing by breaking up the domain into maxGridSize chunks. This should not happen!");

    Vector<Box> boxes;
    Vector<int> procs;
  
    domainSplit(a_domain, boxes, a_maxGridSize);
    LoadBalancing::sort(boxes, m_boxSorting);
    LoadBalancing::makeBalance(procs, boxes);

    a_grids.define(boxes, procs, a_domain);
  }

  m_timer.stopEvent("Make grids");
}

bool ScanShop::isRegular(const Box a_box, const RealVect a_probLo, const Real a_dx) const {
  CH_TIME("ScanShop::isRegular(Box, RealVect, Real)");

  bool ret = true;
  
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    const RealVect a_point = a_probLo + a_dx*(0.5*RealVect::Unit + RealVect(iv));
    if(m_baseIF->value(a_point) >= -0.5*a_dx){
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
    if(m_baseIF->value(a_point) <= 0.5*a_dx) {
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
  // 2. Search through all the boxes and label them as regular/covered/irregular. 
  // 3. Gather cut-cell and regular/covered boxes separately and load balance them separately
  // 4. Create a new DBL with the newly load-balanced boxes; there should be approximately the same amount
  //    of cut-cell boxes for each rank
  // 5. Copy the map created over the initial DisjointBoxLayout onto the final grid
  
  // 1.
  Vector<Box> boxes;
  Vector<int> procs;
  domainSplit(m_domains[a_level], boxes, a_maxGridSize);
  LoadBalancing::makeBalance(procs, boxes);

  // 2. 
  DisjointBoxLayout dbl(boxes, procs, m_domains[a_level]);

  Vector<Box> coveredBoxes;
  Vector<Box> cutCellBoxes;
  Vector<Box> regularBoxes;  

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box = dbl[dit()];
    
    Box grownBox = box;
    grownBox.grow(m_maxGhostEB);
    grownBox &= m_domains[a_level];

    const bool isRegular = ScanShop::isRegular(grownBox, m_probLo, m_dx[a_level]);
    const bool isCovered = ScanShop::isCovered(grownBox, m_probLo, m_dx[a_level]);

    if(isCovered && !isRegular){
      coveredBoxes.push_back(box);
    }    
    else if(isRegular && !isCovered){
      regularBoxes.push_back(box);
    }
    else if(!isRegular && !isCovered){
      cutCellBoxes.push_back(box);
    }
    else{
      MayDay::Error("ScanShop::buildCoarseLevel - logic bust");
    }
  }

  this->defineLevel(coveredBoxes, regularBoxes, cutCellBoxes, a_level);

  m_hasThisLevel[a_level] = 1;
}

void ScanShop::buildFinerLevels(const int a_coarserLevel, const int a_maxGridSize){
  CH_TIME("ScanShop::buildFinerLevels(int, int)");

  if(a_coarserLevel > 0){

    const int coarLvl = a_coarserLevel;
    const int fineLvl = coarLvl - 1;

    // Coar stuff
    Vector<Box> coveredBoxes;
    Vector<Box> regularBoxes;    
    Vector<Box> cutCellBoxes;

    // Find out which boxes were covered/regular/cut on the coarse level. Every box that had cut-cells
    // is split up into new boxes. We will redo these boxes later. 
    const DisjointBoxLayout& dblCoar =  m_grids[coarLvl];        
    for (DataIterator dit(dblCoar); dit.ok(); ++dit){
      const Box coarBox = dblCoar[dit()];
      const Box fineBox = refine(coarBox, 2);

      const GeometryService::InOut& boxType = (*m_boxMap[coarLvl])[dit()];

      if(boxType == GeometryService::Covered){
	coveredBoxes.push_back(fineBox);
      }            
      else if(boxType == GeometryService::Regular){
	regularBoxes.push_back(fineBox);
      }
      else if(boxType == GeometryService::Irregular){
	Vector<Box> boxes;
	domainSplit(fineBox, boxes, a_maxGridSize);
	cutCellBoxes.append(boxes);
      }
    }

    // Make a DBL out of the cut-cell boxes. We distribute these boxes equally among the available ranks so we can efficiently iterate through them again.
    Vector<int> procs;
    LoadBalancing::gatherBoxes(cutCellBoxes);
    LoadBalancing::sort(cutCellBoxes, m_boxSorting);
    LoadBalancing::makeBalance(procs, cutCellBoxes);
    DisjointBoxLayout CutCellDBL(cutCellBoxes, procs, m_domains[fineLvl]);

    // Redo the cut cell boxes
    cutCellBoxes.resize(0);
    for (DataIterator dit(CutCellDBL); dit.ok(); ++dit){
      const Box box = CutCellDBL[dit()];
      
      Box grownBox = box;
      grownBox.grow(m_maxGhostEB);
      grownBox &= m_domains[fineLvl];

      const bool isRegular = ScanShop::isRegular(grownBox, m_probLo, m_dx[fineLvl]);
      const bool isCovered = ScanShop::isCovered(grownBox, m_probLo, m_dx[fineLvl]);

      // Designate boxes and regular/covered or irregular. For regular covered boxes we use a constant load = 1 and for the cut-cell boxes
      // we compute the load introspectively. We do this because patches can contain a different amount of cut-cells and thus there's inherent
      // load imbalance. Worse, the cost of the implicit function can vary in space and thus we time these things introspectively.
      if(isCovered){
	coveredBoxes.push_back(box);
      }      
      else if(isRegular){
	regularBoxes.push_back(box);
      }
      else if(!isRegular && !isCovered){
	cutCellBoxes.push_back(box);
      }
      else{
	MayDay::Error("ScanShop::buildFinerLevels - logic bust!");
      }
    }

    this->defineLevel(coveredBoxes, regularBoxes, cutCellBoxes, fineLvl);
    
    m_hasThisLevel[fineLvl] = 1;

    // Now recurse into finer levels. 
    ScanShop::buildFinerLevels(fineLvl, a_maxGridSize);
  }
}

void ScanShop::defineLevel(Vector<Box>& a_coveredBoxes,
			   Vector<Box>& a_regularBoxes,
			   Vector<Box>& a_cutCellBoxes,
			   const int    a_level){
  CH_TIME("ScanShop::defineLevel");

  // Gather boxes and loads for regular/covered and cut-cell regions.
  LoadBalancing::gatherBoxes(a_coveredBoxes);        
  LoadBalancing::gatherBoxes(a_regularBoxes);    
  LoadBalancing::gatherBoxes(a_cutCellBoxes);

  LoadBalancing::sort(a_coveredBoxes, m_boxSorting);
  LoadBalancing::sort(a_regularBoxes, m_boxSorting);
  LoadBalancing::sort(a_cutCellBoxes, m_boxSorting);

  const Vector<int> coveredTypes(a_coveredBoxes.size(), 0);
  const Vector<int> regularTypes(a_regularBoxes.size(), 1);
  const Vector<int> cutCellTypes(a_cutCellBoxes.size(), 2);    

  // Load balance the boxes. 
  const Vector<long> coveredLoads(a_coveredBoxes.size(), 1L);
  const Vector<long> regularLoads(a_regularBoxes.size(), 1L);
  const Vector<long> cutCellLoads(a_cutCellBoxes.size(), 1L);

  Vector<int> coveredProcs;
  Vector<int> regularProcs;    
  Vector<int> cutCellProcs;    
    
  LoadBalancing::makeBalance(coveredProcs, coveredLoads, a_coveredBoxes);
  LoadBalancing::makeBalance(regularProcs, regularLoads, a_regularBoxes);
  LoadBalancing::makeBalance(cutCellProcs, cutCellLoads, a_cutCellBoxes);    

  // We load balanced the regular/covered and cut-cell regions independently, but now we need to create a box-to-rank map
  // that is usable by Chombo's DisjointBoxLayout.
  Vector<Box> allBoxes;
  Vector<int> allProcs;
  Vector<int> allTypes;    
    
  allBoxes.append(a_coveredBoxes);
  allBoxes.append(a_regularBoxes);
  allBoxes.append(a_cutCellBoxes);    
    
  allProcs.append(coveredProcs);
  allProcs.append(regularProcs);
  allProcs.append(cutCellProcs);

  allTypes.append(coveredTypes);
  allTypes.append(regularTypes);
  allTypes.append(cutCellTypes);

  // This is something that I fucking HATE, but DisjointBoxLayout sorts the boxes using lexicographical sorting, and we must do the same if
  // we want to be able to globally index correctly into the box types. So, create a view of the boxes and box types that is consistent with
  // what we will have in the DataIterator. This means that we must lexicographically sort the boxes. 
  const std::vector<std::pair<Box, int> > sortedBoxesAndTypes = this->getSortedBoxesAndTypes(allBoxes, allTypes);

  // Define shit. 
  m_grids [a_level] = DisjointBoxLayout(allBoxes, allProcs, m_domains[a_level]);
  m_boxMap[a_level] = RefCountedPtr<LayoutData<GeometryService::InOut> > (new LayoutData<GeometryService::InOut>(m_grids[a_level]));    

  for (DataIterator dit(m_grids[a_level]); dit.ok(); ++dit){
    const Box  box  = sortedBoxesAndTypes[dit().intCode()].first;    
    const long type = sortedBoxesAndTypes[dit().intCode()].second;

    // This is an error.
    if(box != m_grids[a_level][dit()]){
      MayDay::Error("ScanShop::defineLevel -- logic bust, boxes should be the same!");
    }
    
    // Otherwise we are fine, set the map to what it should be. 
    if(type == 0L){
      (*m_boxMap[a_level])[dit()] = GeometryService::Covered;
    }
    else if(type == 1L){
      (*m_boxMap[a_level])[dit()] = GeometryService::Regular;
    }
    else if(type == 2L){
      (*m_boxMap[a_level])[dit()] = GeometryService::Irregular;
    }
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

    MayDay::Error("ScanShop::InsideOutside - logic bust 1");
  }

  GeometryService::InOut ret;

  if(foundLevel && m_hasThisLevel[whichLevel]){
    ret = (*m_boxMap[whichLevel])[a_dit];
  }
  else{
    MayDay::Error("ScanShop::InsideOutSide -- logic bust 2");
    
    ret = GeometryService::InsideOutside(a_region, domain, a_probLo, a_dx, a_dit);
  }

  return ret;
}

void ScanShop::fillGraph(BaseFab<int>&        a_regIrregCovered,
                         Vector<IrregNode>&   a_nodes,
                         const Box&           a_validRegion,
                         const Box&           a_ghostRegion,
                         const ProblemDomain& a_domain,
                         const RealVect&      a_probLo,
                         const Real&          a_dx,
                         const DataIndex&     a_di) const {
  CH_TIME("ScanShop::fillGraph");

  if(m_profile){
    m_timer.startEvent("Fill graph");
  }

  GeometryShop::fillGraph(a_regIrregCovered, a_nodes, a_validRegion, a_ghostRegion, a_domain, a_probLo, a_dx, a_di);

  if(m_profile){
    m_timer.stopEvent("Fill graph");
  }
}

std::vector<std::pair<Box, int> > ScanShop::getSortedBoxesAndTypes(const Vector<Box>& a_boxes, const Vector<int>& a_types) const {

  std::vector<std::pair<Box, int> > sortedBoxesAndTypes;
  
  for (int i = 0; i < a_boxes.size(); i++){
    sortedBoxesAndTypes.emplace_back(std::make_pair(a_boxes[i], a_types[i]));
  }


  auto comparator = [](const std::pair<Box, int>& a, const std::pair<Box, int>& b) -> bool {
    return a.first < b.first;
  };
  
  std::sort(sortedBoxesAndTypes.begin(), sortedBoxesAndTypes.end(), comparator);

  return sortedBoxesAndTypes;
}

#include <CD_NamespaceFooter.H>
