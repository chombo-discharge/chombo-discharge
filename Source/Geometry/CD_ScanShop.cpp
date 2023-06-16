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

ScanShop::ScanShop(const BaseIF&       a_localGeom,
                   const int           a_verbosity,
                   const Real          a_dx,
                   const RealVect      a_probLo,
                   const ProblemDomain a_finestDomain,
                   const ProblemDomain a_scanLevel,
                   const int           a_ebGhost,
                   const Real          a_thrshdVoF)
  : GeometryShop(a_localGeom, a_verbosity, a_dx * RealVect::Unit, a_thrshdVoF)
{

  CH_TIME("ScanShop::ScanShop(BaseIF, int, Real, RealVect, ProblemDomain, ProblemDomain, int, Real)");

  m_baseIF       = &a_localGeom;
  m_hasScanLevel = false;
  m_profile      = false;
  m_ebGhost      = a_ebGhost;
  m_fileName     = "ScanShopReport.dat";
  m_boxSorting   = BoxSorting::Morton;
  m_threshold    = a_thrshdVoF;

  // EBISLevel doesn't give resolution, origin, and problem domains through makeGrids, so we
  // need to construct these here, and then extract the proper resolution when we actually call makeGrids
  ScanShop::makeDomains(a_dx, a_probLo, a_finestDomain, a_scanLevel);

  // These are settings that I want to hide from the user.
  ParmParse pp("ScanShop");

  std::string str;
  pp.query("profile", m_profile);
  pp.query("box_sorting", str);

  if (str == "none") {
    m_boxSorting = BoxSorting::None;
  }
  else if (str == "std") {
    m_boxSorting = BoxSorting::Std;
  }
  else if (str == "morton") {
    m_boxSorting = BoxSorting::Morton;
  }
  else if (str == "shuffle") {
    m_boxSorting = BoxSorting::Shuffle;
  }

  m_timer = Timer("ScanShop");
}

ScanShop::~ScanShop()
{
  CH_TIME("ScanShop::~ScanShop()");

  if (m_profile) {
    m_timer.writeReportToFile(m_fileName);
    m_timer.eventReport(pout(), false);
  }
}

void
ScanShop::setProfileFileName(const std::string a_fileName)
{
  m_fileName = a_fileName;

  m_timer = Timer(m_fileName);
}

void
ScanShop::makeDomains(const Real          a_dx,
                      const RealVect      a_probLo,
                      const ProblemDomain a_finestDomain,
                      const ProblemDomain a_scanLevel)
{
  CH_TIME("ScanShop::makeDomains(Real, RealVect, ProblemDomain, ProblemDomain)");

  m_probLo = a_probLo;

  m_dx.resize(0);
  m_domains.resize(0);

  m_dx.push_back(a_dx);
  m_domains.push_back(a_finestDomain);

  const int ref = 2;
  for (int lvl = 0;; lvl++) {
    Real          dx     = m_dx[lvl];
    ProblemDomain domain = m_domains[lvl];

    if (a_scanLevel.domainBox() == domain.domainBox()) {
      m_scanLevel = lvl;
    }

    if (domain.domainBox().coarsenable(ref)) {
      domain.coarsen(ref);
      dx *= ref;

      m_dx.push_back(dx);
      m_domains.push_back(domain);
    }
    else {
      break;
    }
  }

  // These will be built when they are needed.
  m_grids.resize(m_domains.size());
  m_boxMap.resize(m_domains.size());
  m_hasThisLevel.resize(m_domains.size(), false);
}

void
ScanShop::makeGrids(const ProblemDomain& a_domain,
                    DisjointBoxLayout&   a_grids,
                    const int&           a_maxGridSize,
                    const int&           a_maxIrregGridSize)
{
  CH_TIME("ScanShop::makeGrids(ProblemDomain, DisjointBoxLayout, int, int)");
  m_timer.startEvent("Make grids");

  // Build the scan level first
  if (!m_hasScanLevel) {
    m_timer.startEvent("Build coarse");
    for (int lvl = m_domains.size() - 1; lvl >= m_scanLevel; lvl--) {
      ScanShop::buildCoarseLevel(lvl, a_maxGridSize); // Coarser levels built in the same way as the scan level
    }
    m_timer.stopEvent("Build coarse");
    m_timer.startEvent("Build fine");
    ScanShop::buildFinerLevels(m_scanLevel, a_maxGridSize); // Traverse towards finer levels
    m_timer.stopEvent("Build fine");

    m_hasScanLevel = true;
  }

  // Find the level corresponding to a_domain
  int whichLevel;
  for (int lvl = 0; lvl < m_domains.size(); lvl++) {
    if (m_domains[lvl].domainBox() == a_domain.domainBox()) {
      whichLevel = lvl;
      break;
    }
  }

  if (m_hasThisLevel[whichLevel]) {
    a_grids = m_grids[whichLevel];
  }
  else {
    // Development code. Break up a_domnain in a_maxGridSize chunks, load balance trivially and return the dbl

    MayDay::Warning(
      "ScanShop::makeGrids -- decomposing by breaking up the domain into maxGridSize chunks. This should not happen!");

    Vector<Box> boxes;
    Vector<int> procs;

    domainSplit(a_domain, boxes, a_maxGridSize);
    LoadBalancing::sort(boxes, m_boxSorting);
    LoadBalancing::makeBalance(procs, boxes);

    a_grids.define(boxes, procs, a_domain);
  }

  m_timer.stopEvent("Make grids");
}

void
ScanShop::buildCoarseLevel(const int a_level, const int a_maxGridSize)
{
  CH_TIME("ScanShop::buildCoarseLevel(int, int)");

  // This function does the following:
  // 1. Break up the domain into chunks of max size a_maxGridSize and load balance them based on their volumes.
  // 2. Search through all the boxes and label them as regular/covered/irregular.
  // 3. Call defineLevel which creates a "map" of the grid boxes.

  Vector<Box> boxes;
  Vector<int> procs;
  domainSplit(m_domains[a_level], boxes, a_maxGridSize);
  LoadBalancing::makeBalance(procs, boxes);

  // 2.
  DisjointBoxLayout dbl(boxes, procs, m_domains[a_level]);

  Vector<Box> coveredBoxes;
  Vector<Box> cutCellBoxes;
  Vector<Box> regularBoxes;

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box box      = dbl[dit()];
    const Box grownBox = grow(box, m_ebGhost) & m_domains[a_level];

    const bool isRegular = ScanShop::isRegular(grownBox, m_probLo, m_dx[a_level]);
    const bool isCovered = ScanShop::isCovered(grownBox, m_probLo, m_dx[a_level]);

    if (isCovered && !isRegular) {
      coveredBoxes.push_back(box);
    }
    else if (isRegular && !isCovered) {
      regularBoxes.push_back(box);
    }
    else if (!isRegular && !isCovered) {
      cutCellBoxes.push_back(box);
    }
    else {
      MayDay::Error("ScanShop::buildCoarseLevel - logic bust");
    }
  }

  this->defineLevel(coveredBoxes, regularBoxes, cutCellBoxes, a_level);

  m_hasThisLevel[a_level] = true;
}

void
ScanShop::buildFinerLevels(const int a_coarserLevel, const int a_maxGridSize)
{
  CH_TIME("ScanShop::buildFinerLevels(int, int)");

  if (a_coarserLevel > 0) {

    m_timer.startEvent("Fine from coar");
    const int coarLvl = a_coarserLevel;
    const int fineLvl = coarLvl - 1;

    // Coar stuff
    Vector<Box> coveredBoxes;
    Vector<Box> regularBoxes;
    Vector<Box> cutCellBoxes;

    // Find out which boxes were covered/regular/cut on the coarse level. Every box that had cut-cells
    // is split up into new boxes. We will redo these boxes later.
    const DisjointBoxLayout& dblCoar = m_grids[coarLvl];
    for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
      const Box coarBox = dblCoar[dit()];
      const Box fineBox = refine(coarBox, 2);

      const GeometryService::InOut& boxType = (*m_boxMap[coarLvl])[dit()];

      const int  yiv       = coarBox.bigEnd()[1];
      const Real y         = m_probLo[1] + yiv * m_dx[coarLvl];
      const bool doThisBox = y < 0.5;
      if (doThisBox) {
        if (boxType == GeometryService::Covered) {
          coveredBoxes.push_back(fineBox);
        }
        else if (boxType == GeometryService::Regular) {
          regularBoxes.push_back(fineBox);
        }
        else if (boxType == GeometryService::Irregular) {
          Vector<Box> boxes;
          domainSplit(fineBox, boxes, a_maxGridSize, a_maxGridSize);

          for (const auto& box : boxes.stdVector()) {
            const Box grownBox = grow(box, m_ebGhost) & m_domains[fineLvl];

            const bool isRegular = ScanShop::isRegular(grownBox, m_probLo, m_dx[fineLvl]);
            const bool isCovered = ScanShop::isCovered(grownBox, m_probLo, m_dx[fineLvl]);

            if (isCovered) {
              coveredBoxes.push_back(box);
            }
            else if (isRegular) {
              regularBoxes.push_back(box);
            }
            else if (!isRegular && !isCovered) {
              cutCellBoxes.push_back(box);
            }
            else {
              MayDay::Error("ScanShop::buildFinerLevels - logic bust!");
            }
          }
        }
      }
    }
    m_timer.stopEvent("Fine from coar");

    m_timer.startEvent("Define level");
    this->defineLevel(coveredBoxes, regularBoxes, cutCellBoxes, fineLvl);
    m_timer.stopEvent("Define level");

    m_hasThisLevel[fineLvl] = true;

    // Now recurse into finer levels.
    ScanShop::buildFinerLevels(fineLvl, a_maxGridSize);
  }
}

void
ScanShop::defineLevel(Vector<Box>& a_coveredBoxes,
                      Vector<Box>& a_regularBoxes,
                      Vector<Box>& a_cutCellBoxes,
                      const int    a_level)
{
  CH_TIME("ScanShop::defineLevel");

  // Gather boxes and loads for regular/covered and cut-cell regions.
  m_timer.startEvent("Gather boxes");
  LoadBalancing::gatherBoxes(a_coveredBoxes);
  LoadBalancing::gatherBoxes(a_regularBoxes);
  LoadBalancing::gatherBoxes(a_cutCellBoxes);
  m_timer.stopEvent("Gather boxes");

  // This is needed because when we join the regular and covered boxes onto the same "load" when we sort them, but we need to have
  // the indexing correctly. For the cut-cell boxes we load balance them independently.
  m_timer.startEvent("Sort boxes");
  const Vector<int> coveredTypes(a_coveredBoxes.size(), 0);
  const Vector<int> regularTypes(a_regularBoxes.size(), 1);
  const Vector<int> cutCellTypes(a_cutCellBoxes.size(), 2);

  Vector<Box> reguCovBoxes;
  Vector<int> reguCovTypes;

  reguCovBoxes.append(a_coveredBoxes);
  reguCovBoxes.append(a_regularBoxes);

  reguCovTypes.append(coveredTypes);
  reguCovTypes.append(regularTypes);

  //   LoadBalancing::sort(  reguCovBoxes, reguCovTypes, m_boxSorting); Don't need to sort these, I think.
  // We don't need to track the "type" for cut-cell boxes because this call only sorts one type of box.
  LoadBalancing::sort(a_cutCellBoxes, m_boxSorting);
  m_timer.stopEvent("Sort boxes");

  // Load balance the boxes.
  m_timer.startEvent("Make balance");
  const Vector<long> reguCovLoads(reguCovBoxes.size(), 1L);
  const Vector<long> cutCellLoads(a_cutCellBoxes.size(), 1L);

  Vector<int> reguCovProcs;
  Vector<int> cutCellProcs;

  LoadBalancing::makeBalance(reguCovProcs, reguCovLoads, reguCovBoxes);
  LoadBalancing::makeBalance(cutCellProcs, cutCellLoads, a_cutCellBoxes);
  m_timer.stopEvent("Make balance");

  // We load balanced the regular/covered and cut-cell regions independently, but now we need to create a box-to-rank map
  // that is usable by Chombo's DisjointBoxLayout.
  m_timer.startEvent("Vector append");
  Vector<Box> allBoxes;
  Vector<int> allProcs;
  Vector<int> allTypes;

  allBoxes.append(reguCovBoxes);
  allBoxes.append(a_cutCellBoxes);

  allProcs.append(reguCovProcs);
  allProcs.append(cutCellProcs);

  allTypes.append(reguCovTypes);
  allTypes.append(cutCellTypes);
  m_timer.stopEvent("Vector append");

  // DisjointBoxLayout sorts the boxes using lexicographical sorting, and we must do the same if we want to be able
  // to globally index correctly into the box types. So, create a view of the boxes and box types that is consistent with
  // what we will have in the DataIterator. This means that we must lexicographically sort the boxes.
  m_timer.startEvent("Lexi-sort");
  const std::vector<std::pair<Box, int>> sortedBoxesAndTypes = this->getSortedBoxesAndTypes(allBoxes, allTypes);
  m_timer.stopEvent("Lexi-sort");

  // Define shit.
  m_timer.startEvent("Define grids");
  m_grids[a_level] = DisjointBoxLayout(allBoxes, allProcs, m_domains[a_level]);
  m_timer.stopEvent("Define grids");
  m_timer.startEvent("Define map");
  m_boxMap[a_level] = RefCountedPtr<LayoutData<GeometryService::InOut>>(
    new LayoutData<GeometryService::InOut>(m_grids[a_level]));
  m_timer.stopEvent("Define map");

  m_timer.startEvent("Set box types");
  for (DataIterator dit(m_grids[a_level]); dit.ok(); ++dit) {
    const Box box  = sortedBoxesAndTypes[dit().intCode()].first;
    const int type = sortedBoxesAndTypes[dit().intCode()].second;

    // This is an error.
    if (box != m_grids[a_level][dit()]) {
      MayDay::Error("ScanShop::defineLevel -- logic bust, boxes should be the same!");
    }

    // Otherwise we are fine, set the map to what it should be.
    if (type == 0) {
      (*m_boxMap[a_level])[dit()] = GeometryService::Covered;
    }
    else if (type == 1) {
      (*m_boxMap[a_level])[dit()] = GeometryService::Regular;
    }
    else if (type == 2) {
      (*m_boxMap[a_level])[dit()] = GeometryService::Irregular;
    }
    else {
      MayDay::Error("ScanShop::defineLevel - logic bust, did not find type!");
    }
  }

  if (m_profile) {
    pout() << "ScanShop::defineLevel  domain = " << m_domains[a_level] << ":" << endl
           << "\t Covered  boxes = " << a_coveredBoxes.size() << endl
           << "\t Regular  boxes = " << a_regularBoxes.size() << endl
           << "\t Cut-cell boxes = " << a_cutCellBoxes.size() << endl
           << endl;
  }

  m_timer.stopEvent("Set box types");
}

GeometryService::InOut
ScanShop::InsideOutside(const Box&           a_region,
                        const ProblemDomain& a_domain,
                        const RealVect&      a_probLo,
                        const Real&          a_dx,
                        const DataIndex&     a_dit) const
{
  CH_TIME("ScanShop::InsideOutSide(Box, ProblemDomain, RealVect, Real, DataIndex)");

  // Find the level corresponding to a_domain
  int  whichLevel = -1;
  bool foundLevel = false;
  for (int lvl = 0; lvl < m_domains.size(); lvl++) {
    if (m_domains[lvl].domainBox() == a_domain.domainBox()) {
      whichLevel = lvl;
      foundLevel = true;
      break;
    }
  }

  // A strang but true thing. This function is used in EBISLevel::simplifyGraphFromGeo and that function can send in a_domain
  // and a_dx on different levels....
  ProblemDomain domain;
  if (a_dx < m_dx[whichLevel] && whichLevel > 0) {
    domain = m_domains[whichLevel - 1];

    MayDay::Error("ScanShop::InsideOutside - logic bust 1");
  }

  GeometryService::InOut ret;

  if (foundLevel && m_hasThisLevel[whichLevel]) {
    ret = (*m_boxMap[whichLevel])[a_dit];
  }
  else {
    MayDay::Error("ScanShop::InsideOutSide -- logic bust 2");

    ret = GeometryService::InsideOutside(a_region, domain, a_probLo, a_dx, a_dit);
  }

  return ret;
}

void
ScanShop::fillGraph(BaseFab<int>&        a_regIrregCovered,
                    Vector<IrregNode>&   a_nodes,
                    const Box&           a_validRegion,
                    const Box&           a_ghostRegion,
                    const ProblemDomain& a_domain,
                    const RealVect&      a_probLo,
                    const Real&          a_dx,
                    const DataIndex&     a_dit) const
{
  CH_TIMERS("ScanShop::fillGraph");
  CH_TIMER("ScanShop::fillGraph::identify_cells", t1);
  CH_TIMER("ScanShop::fillGraph::make_ndoes", t2);
  CH_TIMER("ScanShop::fillGraph::remove_threshold_vofs", t3);

  CH_assert(a_domain.contains(a_ghostRegion));
  CH_assert(a_ghostRegion.contains(a_validRegion));

  if (m_profile) {
    m_timer.startEvent("Fill graph");
  }

  const RealVect vectDx = a_dx * RealVect::Unit;

  IntVectSet ivsIrreg = IntVectSet(DenseIntVectSet(a_ghostRegion, false));
  IntVectSet ivsDrop  = IntVectSet(DenseIntVectSet(a_ghostRegion, false)); // CP

  CH_START(t1);
  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    const GeometryService::InOut inout = this->insideOutside(iv, a_domain, a_probLo, vectDx);

    if (inout == GeometryService::Covered) {
      a_regIrregCovered(iv, 0) = -1;
    }
    else if (inout == GeometryService::Regular) {
      a_regIrregCovered(iv, 0) = 1;
    }
    else {
      a_regIrregCovered(iv, 0) = 0;
      if (a_validRegion.contains(iv)) {
        ivsIrreg |= iv;
      }
    }
  }

  //if a regular is next to a  covered, change to irregular with correct arcs and so on.
  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (a_regIrregCovered(iv, 0) == -1) {
      fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);
    }
  }
  CH_STOP(t1);

  // now loop through irregular cells and make nodes for each  one.
  CH_START(t2);

  for (IVSIterator ivsit(ivsIrreg); ivsit.ok(); ++ivsit) {
    const IntVect iv = ivsit();
    VolIndex      vof(iv, 0);

    Real             volFrac, bndryArea;
    RealVect         normal, volCentroid, bndryCentroid;
    Vector<int>      loArc[SpaceDim];
    Vector<int>      hiArc[SpaceDim];
    Vector<Real>     loAreaFrac[SpaceDim];
    Vector<Real>     hiAreaFrac[SpaceDim];
    Vector<RealVect> loFaceCentroid[SpaceDim];
    Vector<RealVect> hiFaceCentroid[SpaceDim];
    computeVoFInternals(volFrac,
                        loArc,
                        hiArc,
                        loAreaFrac,
                        hiAreaFrac,
                        bndryArea,
                        normal,
                        volCentroid,
                        bndryCentroid,
                        loFaceCentroid,
                        hiFaceCentroid,
                        a_regIrregCovered,
                        ivsIrreg,
                        vof,
                        a_domain,
                        a_probLo,
                        a_dx,
                        vectDx,
                        ivsit());

    if (m_threshold > 0. && volFrac < m_threshold) {
      ivsDrop |= ivsit();
      a_regIrregCovered(ivsit(), 0) = -1;
      pout() << "Removing vof " << vof << " with volFrac " << volFrac << endl;
    }
    else {
      IrregNode newNode;
      newNode.m_cell          = ivsit();
      newNode.m_volFrac       = volFrac;
      newNode.m_cellIndex     = 0;
      newNode.m_volCentroid   = volCentroid;
      newNode.m_bndryCentroid = bndryCentroid;

      for (int faceDir = 0; faceDir < SpaceDim; faceDir++) {
        int loNodeInd                     = newNode.index(faceDir, Side::Lo);
        int hiNodeInd                     = newNode.index(faceDir, Side::Hi);
        newNode.m_arc[loNodeInd]          = loArc[faceDir];
        newNode.m_arc[hiNodeInd]          = hiArc[faceDir];
        newNode.m_areaFrac[loNodeInd]     = loAreaFrac[faceDir];
        newNode.m_areaFrac[hiNodeInd]     = hiAreaFrac[faceDir];
        newNode.m_faceCentroid[loNodeInd] = loFaceCentroid[faceDir];
        newNode.m_faceCentroid[hiNodeInd] = hiFaceCentroid[faceDir];
      }
      a_nodes.push_back(newNode);
    }
  }
  CH_STOP(t2);

  CH_START(t3);
  // Sweep that removes cells with volFrac less than a certain threshold
  for (IVSIterator ivsit(ivsDrop); ivsit.ok(); ++ivsit) {
    VolIndex vof(ivsit(), 0);
    IntVect  iv = vof.gridIndex();
    for (int faceDir = 0; faceDir < SpaceDim; faceDir++) {
      for (SideIterator sit; sit.ok(); ++sit) {
        int     isign   = sign(sit());
        IntVect otherIV = iv + isign * BASISV(faceDir);
        if (a_validRegion.contains(otherIV)) {
          if (a_regIrregCovered(otherIV, 0) == 0) {
            // i am in the case where the other cell
            // is also irregular.   I just made a previously
            // irregular cell covered so I have to check to
            // see if it had any faces pointed this way.
            int  inode = -1;
            bool found = false;
            for (int ivec = 0; ivec < a_nodes.size() && !found; ivec++) {
              if (a_nodes[ivec].m_cell == otherIV) {
                inode = ivec;
                found = true;
              }
            }
            if (!found && a_validRegion.contains(otherIV)) {
              MayDay::Error("something wrong in our logic");
            }
            if (found) {
              int arcindex = a_nodes[inode].index(faceDir, flip(sit()));
              a_nodes[inode].m_arc[arcindex].resize(0);
              a_nodes[inode].m_areaFrac[arcindex].resize(0);
              a_nodes[inode].m_faceCentroid[arcindex].resize(0);
            }
          }
        }
      }
    }

    fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);
  }
  CH_STOP(t3);

  if (m_profile) {
    m_timer.stopEvent("Fill graph");
  }
}

GeometryService::InOut
ScanShop::insideOutside(const IntVect        a_iv,
                        const ProblemDomain& a_domain,
                        const RealVect       a_probLo,
                        const RealVect       a_dx) const noexcept
{
  CH_TIME("ScanShop::insideOutside");

  // TLDR: To classify a cell as completely inside/outside we assume that the value of the implicit function at all edges must be the same.
  GeometryService::InOut inOut = GeometryService::Irregular;

  // Compute implicit function value at all corners
  const Box  corners    = surroundingNodes(Box(a_iv, a_iv));
  const Real firstValue = m_baseIF->value(a_probLo + corners.smallEnd() * a_dx);

  if (firstValue < 0.0) {
    inOut = GeometryService::Regular;
  }
  else if (firstValue > 0.0) {
    inOut = GeometryService::Covered;
  }

  if (inOut != GeometryService::Irregular) {
    for (BoxIterator bit(corners); bit.ok(); ++bit) {
      const Real value = m_baseIF->value(a_probLo + bit() * a_dx);

      if (firstValue * value <= 0.0) {
        inOut = GeometryService::Irregular;

        break;
      }
    }
  }

  return inOut;
}

#include <CD_NamespaceFooter.H>
