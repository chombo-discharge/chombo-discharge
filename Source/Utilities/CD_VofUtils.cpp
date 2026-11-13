/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_VofUtils.cpp
 * @brief  Implementation of CD_VofUtils.H
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

// Chombo includes
#include <BoxIterator.H>

// Our includes
#include <CD_VofUtils.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

Vector<VolIndex>
VofUtils::getVofsInRadius(const VolIndex&    a_startVof,
                          const EBISBox&     a_ebisbox,
                          const int          a_radius,
                          const Connectivity a_connectivity,
                          const bool         a_addStartVof)
{
  Vector<VolIndex> vofs;
  switch (a_connectivity) {
  case VofUtils::Connectivity::MonotonePath: {
    vofs = VofUtils::getVofsInMonotonePath(a_startVof, a_ebisbox, a_radius, a_addStartVof);

    break;
  }
  case VofUtils::Connectivity::SimplyConnected: {
    vofs = VofUtils::getConnectedVofsInRadius(a_startVof, a_ebisbox, a_radius, a_addStartVof);

    break;
  }
  case VofUtils::Connectivity::All: {
    vofs = VofUtils::getAllVofsInRadius(a_startVof, a_ebisbox, a_radius, a_addStartVof);

    break;
  }
  }

  return vofs;
}

Vector<VolIndex>
VofUtils::getVofsInQuadrant(const VolIndex&    a_startVof,
                            const EBISBox&     a_ebisbox,
                            const RealVect&    a_normal,
                            const int          a_radius,
                            const Connectivity a_connectivity,
                            const bool         a_addStartVof)
{

  Vector<VolIndex> vofs;

  if (a_normal != RealVect::Zero) {
    vofs = VofUtils::getVofsInRadius(a_startVof, a_ebisbox, a_radius, a_connectivity, false);
    vofs = VofUtils::restrictVofsToQuadrant(vofs, a_startVof, a_ebisbox, a_normal, a_radius);

    if (a_addStartVof) {
      vofs.push_back(a_startVof);
    }
  }

  return vofs;
}

Vector<VolIndex>
VofUtils::restrictVofsToQuadrant(const Vector<VolIndex>& /*a_vofs*/,
                                 const VolIndex& a_startVof,
                                 const EBISBox&  a_ebisbox,
                                 const RealVect& a_normal,
                                 const int       a_radius)
{
  Vector<VolIndex> vofs;

  // Find the quadrant or "half-plane" where the quadrant vofs live, and restrict to that quadrant
  Box quadBox;
  if (VofUtils::isQuadrantWellDefined(a_normal)) {
    quadBox = VofUtils::getQuadrant(a_normal, a_startVof, a_ebisbox, a_radius);
  }
  else {
    const std::pair<int, Side::LoHiSide> cardinal = VofUtils::getCardinalDirection(a_normal);
    quadBox = VofUtils::getSymmetricQuadrant(cardinal, a_startVof, a_ebisbox, a_radius);
  }

  VofUtils::includeCells(vofs, quadBox);

  return vofs;
}

Vector<VolIndex>
VofUtils::getVofsInSemiCircle(const VolIndex&      a_startVof,
                              const EBISBox&       a_ebisbox,
                              const RealVect&      a_normal,
                              const int            a_radius,
                              const Real           a_deltaThresh,
                              const Connectivity   a_connectivity,
                              const Location::Cell a_vofLocation,
                              const Location::Cell a_cellLocation,
                              const bool           a_addStartVof)
{
  Vector<VolIndex> vofs;

  if (a_normal != RealVect::Zero) {
    vofs = VofUtils::getVofsInRadius(a_startVof, a_ebisbox, a_radius, a_connectivity, false);
    vofs = VofUtils::restrictVofsToSemiCircle(vofs,
                                              a_startVof,
                                              a_ebisbox,
                                              a_normal,
                                              a_deltaThresh,
                                              a_vofLocation,
                                              a_cellLocation);
    if (a_addStartVof) {
      vofs.push_back(a_startVof);
    }
  }

  return vofs;
}

Vector<VolIndex>
VofUtils::restrictVofsToSemiCircle(const Vector<VolIndex>& a_vofs,
                                   const VolIndex&         a_startVof,
                                   const EBISBox&          a_ebisbox,
                                   const RealVect&         a_normal,
                                   const Real              a_deltaThresh,
                                   const Location::Cell    a_vofLocation,
                                   const Location::Cell    a_cellLocation)
{
  Vector<VolIndex> vofs;

  if (a_normal != RealVect::Zero) {
    const RealVect x0 = Location::position(a_vofLocation, a_startVof, a_ebisbox, 1.0);

    for (int i = 0; i < a_vofs.size(); i++) {
      const VolIndex& curVoF = a_vofs[i];
      const RealVect  delta  = Location::position(a_cellLocation, curVoF, a_ebisbox, 1.0) - x0;

      if (delta.dotProduct(a_normal) > a_deltaThresh) {
        vofs.push_back(curVoF);
      }
    }
  }

  return vofs;
}

Vector<VolIndex>
VofUtils::getVofsInMonotonePath(const VolIndex& a_startVoF,
                                const EBISBox&  a_ebisbox,
                                const int       a_radius,
                                const bool      a_addStartVof)
{
  IntVect timesMoved = IntVect::Zero;
  IntVect pathSign   = IntVect::Zero;

  const ProblemDomain& domain   = a_ebisbox.getDomain();
  const Box&           region   = a_ebisbox.getRegion();
  const Box            validBox = domain & region;

#if 0 // Old way
  Vector<VolIndex> vofList;
  VofUtils::getVofsInMonotonePath(vofList, a_startVoF, a_ebisbox, a_radius, timesMoved, pathSign);

  Vector<VolIndex> ret;  
  if(!a_addStartVof){
    for (int i = 0; i < vofList.size(); i++){
      if(vofList[i] != a_startVoF) ret.push_back(vofList[i]);
    }
  }
  else{
    ret = vofList;
  }

  return ret;
#else // New way
  // Get the Vofs in the monotone path
  std::set<VolIndex> vofSet;
  VofUtils::getVofsInMonotonePath(vofSet, a_startVoF, a_ebisbox, validBox, a_radius, timesMoved, pathSign);

  // Build the return vector.
  Vector<VolIndex> ret;
  for (const auto& vof : vofSet) {
    if (vof != a_startVoF) {
      ret.push_back(vof);
    }
  }
  if (a_addStartVof) {
    ret.push_back(a_startVoF);
  }

  return ret;

#endif
}

Vector<VolIndex>
VofUtils::getConnectedVofsInRadius(const VolIndex& a_startVof,
                                   const EBISBox&  a_ebisbox,
                                   const int       a_radius,
                                   const bool      a_addStartVof)
{

  // Simple: First get all the vofs in a radius, and then only get the ones that are connected to a_startVof
  Vector<VolIndex> allVofs       = VofUtils::getAllVofsInRadius(a_startVof, a_ebisbox, a_radius, true);
  Vector<VolIndex> connectedVofs = VofUtils::connectedVofsOnly(a_startVof, allVofs, a_ebisbox);

  if (a_addStartVof) {
    connectedVofs.push_back(a_startVof);
  }

  return connectedVofs;
}

Vector<VolIndex>
VofUtils::getAllVofsInRadius(const VolIndex& a_startVof,
                             const EBISBox&  a_ebisbox,
                             const int       a_radius,
                             const bool      a_addStartVof)
{
  Vector<VolIndex> allVofs;
  Box              bx(a_startVof.gridIndex(), a_startVof.gridIndex());
  bx.grow(a_radius);
  bx &= a_ebisbox.getDomain();
  bx &= a_ebisbox.getRegion();

  // Not auto-vectorizable: this gathers cut-cell VoFs through the out-of-line ebisbox.getVoFs(iv)
  // (which allocates a Vector per cell) and appends them to a growing Vector. The access is
  // irregular, allocating, and order-dependent, so the loop is left as a plain cell traversal.
  auto kernel = [&](const IntVect& iv) -> void {
    const std::vector<VolIndex> vofsInThisCell = a_ebisbox.getVoFs(iv).stdVector();

    for (const auto& ivof : vofsInThisCell) {
      if (ivof != a_startVof) {
        allVofs.push_back(ivof);
      }
    }
  };

  BoxLoops::loop<D_DECL(1, 1, 1)>(bx, kernel);

  if (a_addStartVof) {
    allVofs.push_back(a_startVof);
  }

  return allVofs;
}

void
VofUtils::excludeCells(Vector<VolIndex>& a_vofs, const Box& a_excludeBox)
{

  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()) {
    if (!a_excludeBox.contains(ivof.gridIndex())) {
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void
VofUtils::excludeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_excludeIVS)
{

  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()) {
    if (!a_excludeIVS.contains(ivof.gridIndex())) {
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void
VofUtils::includeCells(Vector<VolIndex>& a_vofs, const Box& a_includeBox)
{
  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()) {
    if (a_includeBox.contains(ivof.gridIndex())) {
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void
VofUtils::includeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_includeIVS)
{
  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()) {
    if (a_includeIVS.contains(ivof.gridIndex())) {
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void
VofUtils::includeCells(Vector<VolIndex>& a_vofs, const DenseIntVectSet& a_includeIVS)
{
  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()) {
    const IntVect iv = ivof.gridIndex();
    if (a_includeIVS[iv]) {
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void
VofUtils::onlyUnique(Vector<VolIndex>& a_vofs)
{
  std::set<VolIndex> setVofs;
  for (int i = 0; i < a_vofs.size(); i++) {
    setVofs.insert(a_vofs[i]);
  }

  a_vofs.resize(0);
  for (const auto& v : setVofs) {
    a_vofs.push_back(v);
  }
}

Vector<VolIndex>
VofUtils::getFinerVoFs(const VolIndex&   a_coarVoF,
                       const EBISLayout& a_ebislFine,
                       const DataIndex&  a_din,
                       const int         a_refRat)
{
  CH_TIME("VofUtils::getFinerVoFs");

  const Vector<VolIndex> finer = VofUtils::gatherFinerVoFs(a_coarVoF, a_ebislFine, a_din, a_refRat);

  if (finer.size() == 0) {
    return finer;
  }

  return VofUtils::orderFinerVoFs(finer, a_ebislFine, a_din, a_refRat);
}

Vector<FaceIndex>
VofUtils::getFinerFaces(const FaceIndex&  a_coarFace,
                        const EBISLayout& a_ebislFine,
                        const DataIndex&  a_din,
                        const int         a_refRat)
{
  CH_TIME("VofUtils::getFinerFaces");

  const int dir = a_coarFace.direction();

  Vector<FaceIndex> faces;

  if (a_coarFace.isBoundary()) {
    const Side::LoHiSide outside = (a_coarFace.cellIndex(Side::Hi) < 0) ? Side::Hi : Side::Lo;

    const Vector<VolIndex> finer = VofUtils::gatherFinerVoFs(a_coarFace.getVoF(flip(outside)),
                                                             a_ebislFine,
                                                             a_din,
                                                             a_refRat);

    if (finer.size() > 0) {
      VofUtils::appendFinerBoundaryFaces(faces, finer, a_ebislFine, a_din, a_refRat, dir, outside);
    }
  }
  else {
    const Vector<VolIndex> finerLo = VofUtils::gatherFinerVoFs(a_coarFace.getVoF(Side::Lo),
                                                               a_ebislFine,
                                                               a_din,
                                                               a_refRat);
    const Vector<VolIndex> finerHi = VofUtils::gatherFinerVoFs(a_coarFace.getVoF(Side::Hi),
                                                               a_ebislFine,
                                                               a_din,
                                                               a_refRat);

    if (finerLo.size() > 0 && finerHi.size() > 0) {
      VofUtils::appendFinerFaces(faces, finerLo, finerHi, a_ebislFine, a_din, a_refRat, dir);
    }
  }

  return faces;
}

Vector<VolIndex>
VofUtils::gatherFinerVoFs(const VolIndex&   a_coarVoF,
                          const EBISLayout& a_ebislFine,
                          const DataIndex&  a_din,
                          const int         a_refRat)
{
  CH_assert(
    a_ebislFine[a_din].getRegion().contains(refine(Box(a_coarVoF.gridIndex(), a_coarVoF.gridIndex()), a_refRat)));

  return VofUtils::getFinerVoFsInRegion(a_coarVoF, a_ebislFine, a_din, a_refRat);
}

Vector<VolIndex>
VofUtils::getFinerVoFsInRegion(const VolIndex&   a_coarVoF,
                               const EBISLayout& a_ebislFine,
                               const DataIndex&  a_din,
                               const int         a_refRat)
{
  CH_TIME("VofUtils::getFinerVoFsInRegion");

  CH_assert(a_refRat >= 2);
  CH_assert((a_refRat & (a_refRat - 1)) == 0);
  CH_assert(a_ebislFine.getMaxCoarseningRatio() >= a_refRat);

  const EBISBox& ebisBoxFine = a_ebislFine[a_din];
  const IntVect& coarCell    = a_coarVoF.gridIndex();
  const Box      fineCells   = refine(Box(coarCell, coarCell), a_refRat) & ebisBoxFine.getRegion();

  Vector<VolIndex> finer;

  if (fineCells.isEmpty()) {
    return finer;
  }

  for (BoxIterator bit(fineCells); bit.ok(); ++bit) {
    const Vector<VolIndex> vofs = ebisBoxFine.getVoFs(bit());

    for (int i = 0; i < vofs.size(); i++) {
      if (a_ebislFine.coarsen(vofs[i], a_refRat, a_din) == a_coarVoF) {
        finer.push_back(vofs[i]);
      }
    }
  }

  return finer;
}

Vector<Vector<VolIndex>>
VofUtils::orderFinerGroups(const Vector<VolIndex>& a_finest,
                           const EBISLayout&       a_ebislFine,
                           const DataIndex&        a_din,
                           const int               a_ratio)
{
  CH_assert(a_ratio >= 2);
  CH_assert(a_finest.size() > 0);

  const EBISBox& ebisBoxFine = a_ebislFine[a_din];
  const int      childRatio  = a_ratio / 2;

  // The ancestor of each finest VoF one factor of two finer than the VoF they all lie under.
  Vector<VolIndex> ancestors(a_finest.size());

  for (int i = 0; i < a_finest.size(); i++) {
    ancestors[i] = (childRatio == 1) ? a_finest[i] : a_ebislFine.coarsen(a_finest[i], childRatio, a_din);
  }

  // The distinct ancestors are the children, gathered in BoxIterator order of their cells and then by cell index. The
  // children's cells refine one cell by two, so BoxIterator order is the binary number whose bit d is the offset along
  // direction d.
  const IntVect firstCell = 2 * coarsen(ancestors[0].gridIndex(), 2);

#ifndef NDEBUG
  for (int i = 0; i < ancestors.size(); i++) {
    const IntVect offset = ancestors[i].gridIndex() - firstCell;

    CH_assert(offset.min() >= 0 && offset.max() <= 1);
  }
#endif

  const auto gatherKey = [&](const VolIndex& a_child) -> std::pair<int, int> {
    int rank = 0;

    for (int dir = 0; dir < SpaceDim; dir++) {
      rank += (a_child.gridIndex()[dir] - firstCell[dir]) << dir;
    }

    return std::make_pair(rank, a_child.cellIndex());
  };

  Vector<VolIndex> children;

  for (int i = 0; i < ancestors.size(); i++) {
    bool known = false;

    for (int c = 0; c < children.size() && !known; c++) {
      known = children[c] == ancestors[i];
    }

    if (!known) {
      children.push_back(ancestors[i]);
    }
  }

  std::sort(children.stdVector().begin(),
            children.stdVector().end(),
            [&](const VolIndex& a_left, const VolIndex& a_right) -> bool {
              return gatherKey(a_left) < gatherKey(a_right);
            });

  const int numChildren = children.size();

  Vector<Vector<VolIndex>> groups(numChildren);

  for (int i = 0; i < ancestors.size(); i++) {
    for (int c = 0; c < numChildren; c++) {
      if (children[c] == ancestors[i]) {
        groups[c].push_back(a_finest[i]);
      }
    }
  }

  const auto connected = [&](const int a_first, const int a_second) -> bool {
    for (int i = 0; i < groups[a_first].size(); i++) {
      for (int j = 0; j < groups[a_second].size(); j++) {
        if (ebisBoxFine.isConnected(groups[a_first][i], groups[a_second][j])) {
          return true;
        }
      }
    }

    return false;
  };

  // List the children in the order the connected set grows from the first one, as EBGraph::getVoFSets does: sweep the
  // children in gathering order and add every one connected to a child already listed, until a sweep adds none.
  std::vector<int>  order(1, 0);
  std::vector<bool> listed(numChildren, false);

  listed[0] = true;

  bool doneAdding = false;

  while (!doneAdding) {
    const std::size_t oldSize = order.size();

    for (int j = 1; j < numChildren; j++) {
      for (std::size_t k = 0; k < order.size() && !listed[j]; k++) {
        if (connected(order[k], j)) {
          order.push_back(j);

          listed[j] = true;
        }
      }
    }

    doneAdding = order.size() == oldSize;
  }

  // The children are the VoFs one level finer under a single VoF, so a path joins every one of them to the first one.
  // Should one be left over in a graph that breaks this, it is listed last, in gathering order.
  for (int j = 0; j < numChildren; j++) {
    CH_assert(listed[j]);

    if (!listed[j]) {
      order.push_back(j);
    }
  }

  Vector<Vector<VolIndex>> ordered;

  for (std::size_t k = 0; k < order.size(); k++) {
    ordered.push_back(groups[order[k]]);
  }

  return ordered;
}

Vector<VolIndex>
VofUtils::orderFinerVoFs(const Vector<VolIndex>& a_finest,
                         const EBISLayout&       a_ebislFine,
                         const DataIndex&        a_din,
                         const int               a_ratio)
{
  if (a_ratio == 1) {
    CH_assert(a_finest.size() == 1);

    return a_finest;
  }

  const Vector<Vector<VolIndex>> groups = VofUtils::orderFinerGroups(a_finest, a_ebislFine, a_din, a_ratio);

  Vector<VolIndex> ordered;

  for (int g = 0; g < groups.size(); g++) {
    ordered.append(VofUtils::orderFinerVoFs(groups[g], a_ebislFine, a_din, a_ratio / 2));
  }

  return ordered;
}

void
VofUtils::appendFinerFaces(Vector<FaceIndex>&      a_faces,
                           const Vector<VolIndex>& a_finestLo,
                           const Vector<VolIndex>& a_finestHi,
                           const EBISLayout&       a_ebislFine,
                           const DataIndex&        a_din,
                           const int               a_ratio,
                           const int               a_dir)
{
  if (a_ratio == 1) {
    CH_assert(a_finestLo.size() == 1);
    CH_assert(a_finestHi.size() == 1);

    a_faces.push_back(FaceIndex(a_finestLo[0], a_finestHi[0], a_dir));

    return;
  }

  const EBISBox& ebisBoxFine = a_ebislFine[a_din];

  const Vector<Vector<VolIndex>> groupsLo = VofUtils::orderFinerGroups(a_finestLo, a_ebislFine, a_din, a_ratio);
  const Vector<Vector<VolIndex>> groupsHi = VofUtils::orderFinerGroups(a_finestHi, a_ebislFine, a_din, a_ratio);

  // A low and a high child share a face when some pair of their finest VoFs does. Pairs are visited as
  // EBGraph::refine(FaceIndex) visits them: low children outer, high children inner.
  for (int lo = 0; lo < groupsLo.size(); lo++) {
    for (int hi = 0; hi < groupsHi.size(); hi++) {
      bool connected = false;

      for (int i = 0; i < groupsLo[lo].size() && !connected; i++) {
        for (int j = 0; j < groupsHi[hi].size() && !connected; j++) {
          connected = ebisBoxFine.isConnected(groupsLo[lo][i], groupsHi[hi][j]);
        }
      }

      if (connected) {
        VofUtils::appendFinerFaces(a_faces, groupsLo[lo], groupsHi[hi], a_ebislFine, a_din, a_ratio / 2, a_dir);
      }
    }
  }
}

void
VofUtils::appendFinerBoundaryFaces(Vector<FaceIndex>&      a_faces,
                                   const Vector<VolIndex>& a_finest,
                                   const EBISLayout&       a_ebislFine,
                                   const DataIndex&        a_din,
                                   const int               a_ratio,
                                   const int               a_dir,
                                   const Side::LoHiSide    a_side)
{
  const EBISBox& ebisBoxFine = a_ebislFine[a_din];

  if (a_ratio == 1) {
    CH_assert(a_finest.size() == 1);

    a_faces.append(ebisBoxFine.getFaces(a_finest[0], a_dir, a_side));

    return;
  }

  const ProblemDomain& domain = ebisBoxFine.getDomain();
  const IntVect        shift  = sign(a_side) * BASISV(a_dir);

  const Vector<Vector<VolIndex>> groups = VofUtils::orderFinerGroups(a_finest, a_ebislFine, a_din, a_ratio);

  // A child has a boundary face when one of its finest VoFs against the boundary does.
  for (int g = 0; g < groups.size(); g++) {
    bool hasFace = false;

    for (int i = 0; i < groups[g].size() && !hasFace; i++) {
      const VolIndex& vof = groups[g][i];

      hasFace = !domain.contains(vof.gridIndex() + shift) && ebisBoxFine.getFaces(vof, a_dir, a_side).size() > 0;
    }

    if (hasFace) {
      VofUtils::appendFinerBoundaryFaces(a_faces, groups[g], a_ebislFine, a_din, a_ratio / 2, a_dir, a_side);
    }
  }
}

bool
VofUtils::isQuadrantWellDefined(const RealVect& a_normal)
{
  bool ret = true;

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (std::abs(a_normal[dir]) == 1.0 || std::abs(a_normal[dir]) == 0.0) {
      ret = false;
    }
  }

  return ret;
}

std::pair<int, Side::LoHiSide>
VofUtils::getCardinalDirection(const RealVect& a_normal)
{
  std::pair<int, Side::LoHiSide> ret;

  const int cardDir = a_normal.maxDir(true);

  if (a_normal[cardDir] < 0.0) {
    ret = std::make_pair(cardDir, Side::Lo);
  }
  else if (a_normal[cardDir] > 0.0) {
    ret = std::make_pair(cardDir, Side::Hi);
  }
  else {
    MayDay::Error("VofUtils::getCardinalDirection -- I got a zero normal vector so there's no quadrant to be defined!");
  }

  return ret;
}

Box
VofUtils::getQuadrant(const RealVect& a_normal, const VolIndex& a_vof, const EBISBox& a_ebisbox, const Real a_radius)
{

  const IntVect iv = a_vof.gridIndex();

  Box bx(iv, iv);
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (a_normal[dir] < 0.0) {
      bx.growLo(dir, static_cast<int>(a_radius));
    }
    else if (a_normal[dir] > 0.0) {
      bx.growHi(dir, static_cast<int>(a_radius));
    }
    else {
      MayDay::Error(
        "VofUtils::getQuadrant - logic bust. Should not have a normal that aligns with one of the coordinate axes");
    }
  }

  bx &= a_ebisbox.getDomain();

  return bx;
}

Box
VofUtils::getSymmetricQuadrant(const std::pair<int, Side::LoHiSide>& a_cardinal,
                               const VolIndex&                       a_vof,
                               const EBISBox&                        a_ebisbox,
                               const Real                            a_radius)
{
  const IntVect iv = a_vof.gridIndex();

  Box bx(iv, iv);
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir == a_cardinal.first) {
      bx.growDir(dir, a_cardinal.second, static_cast<int>(a_radius));
    }
    else {
      bx.growLo(dir, static_cast<int>(a_radius));
      bx.growHi(dir, static_cast<int>(a_radius));
    }
  }

  bx &= a_ebisbox.getDomain();

  return bx;
}

Vector<VolIndex>
VofUtils::connectedVofsOnly(const VolIndex& a_startVof, const Vector<VolIndex>& a_allVofs, const EBISBox& a_ebisbox)
{

  // If a_allVofs does not contain a_startVof, add it.
  bool                  hasStartingVof = false;
  std::vector<VolIndex> includeVofs;
  for (int i = 0; i < a_allVofs.size(); i++) {
    const VolIndex& ivof = a_allVofs[i];

    if (ivof == a_startVof) {
      hasStartingVof = true;
    }

    includeVofs.emplace_back(ivof);
  }
  if (!hasStartingVof) {
    includeVofs.push_back(a_startVof);
  }

  // Make a map over which Vofs we've added so far
  std::map<VolIndex, bool> vofMap;
  for (const auto& ivof : includeVofs) {
    vofMap.emplace(ivof, false);
  }

  // Insert the first vof
  std::vector<VolIndex> connectedVofs;
  connectedVofs.emplace_back(a_startVof);
  vofMap.at(a_startVof) = true;

  // Add Vofs that are connected to this vof unless they have not already been added.
  // We exit when we can't add more Vofs.
  bool keepGoing = true;
  while (keepGoing) {
    keepGoing = false;

    std::vector<VolIndex> newVofs = connectedVofs;

    for (const auto& ivof : connectedVofs) {

      for (const auto& m : vofMap) {
        const VolIndex otherVof = m.first;

        const bool isConnected  = a_ebisbox.isConnected(ivof, otherVof);
        const bool alreadyAdded = vofMap.at(otherVof);

        if (!alreadyAdded && isConnected) {
          newVofs.emplace_back(otherVof);
          vofMap.at(otherVof) = true;

          keepGoing = true;
        }
      }
    }

    connectedVofs = newVofs;
  }

  Vector<VolIndex> ret;
  for (const auto& ivof : connectedVofs) {
    if (ivof != a_startVof) {
      ret.push_back(ivof);
    }
  }

  return ret;
}

void
VofUtils::getVofsInMonotonePath(Vector<VolIndex>& a_vofList,
                                const VolIndex&   a_startVof,
                                const EBISBox&    a_ebisbox,
                                const int         a_radius,
                                const IntVect&    a_timesMoved,
                                const IntVect&    a_pathSign)
{

  const IntVect        iv     = a_startVof.gridIndex();
  const ProblemDomain& domain = a_ebisbox.getDomain();
  const Box&           region = a_ebisbox.getRegion();

  if (domain.contains(iv) && region.contains(iv)) {

    // Add if not already added
    bool haveStartVof = false;
    for (int ivof = 0; ivof < a_vofList.size(); ivof++) {
      if (a_vofList[ivof] == a_startVof) {
        haveStartVof = true;
      }
    }

    if (!haveStartVof) {
      a_vofList.push_back(a_startVof);

      for (int dir = 0; dir < SpaceDim; dir++) {
        if (a_timesMoved[dir] < a_radius) {
          const IntVect newTimesMoved = a_timesMoved + BASISV(dir);

          // Move in the low direction. If we have already moved in the high direction we are not allowed to "turn
          // back".
          if (a_pathSign[dir] == -1 || a_pathSign[dir] == 0) {
            IntVect newPathSign = a_pathSign;
            newPathSign[dir]    = -1;

            Vector<FaceIndex> faces = a_ebisbox.getFaces(a_startVof, dir, Side::Lo);
            for (const auto& f : faces.stdVector()) {
              const VolIndex& newStartVof = f.getVoF(Side::Lo);

              VofUtils::getVofsInMonotonePath(a_vofList, newStartVof, a_ebisbox, a_radius, newTimesMoved, newPathSign);
            }
          }

          // Move in the high direction. If we have already moved in the low direction we are not allowed to "turn
          // back".
          if (a_pathSign[dir] == 0 || a_pathSign[dir] == 1) {
            IntVect newPathSign = a_pathSign;
            newPathSign[dir]    = 1;

            Vector<FaceIndex> faces = a_ebisbox.getFaces(a_startVof, dir, Side::Hi);
            for (const auto& f : faces.stdVector()) {
              const VolIndex& newStartVof = f.getVoF(Side::Hi);

              VofUtils::getVofsInMonotonePath(a_vofList, newStartVof, a_ebisbox, a_radius, newTimesMoved, newPathSign);
            }
          }
        }
      }
    }
    else {
      // Another path has already visited this vof.
      return;
    }
  }
}

void
VofUtils::getVofsInMonotonePath(std::set<VolIndex>& a_vofSet,
                                const VolIndex&     a_startVof,
                                const EBISBox&      a_ebisbox,
                                const Box&          a_validBox,
                                const int&          a_radius,
                                const IntVect&      a_timesMoved,
                                const IntVect&      a_pathSign)
{

  const IntVect iv = a_startVof.gridIndex();

  if (a_validBox.contains(iv)) {

    // Add if not already added
    bool haveStartVof = false;
    if (a_vofSet.find(a_startVof) != a_vofSet.end()) {
      haveStartVof = true;
    }

    if (!haveStartVof) {
      a_vofSet.insert(a_startVof);

      for (int dir = 0; dir < SpaceDim; dir++) {
        if (a_timesMoved[dir] < a_radius) {
          IntVect newTimesMoved = a_timesMoved;
          newTimesMoved[dir] += 1;

          // Move in the low direction. If we have already moved in the high direction we are not allowed to "turn
          // back".
          if (a_pathSign[dir] <= 0) {
            IntVect newPathSign = a_pathSign;
            newPathSign[dir]    = -1;

            Vector<FaceIndex> faces = a_ebisbox.getFaces(a_startVof, dir, Side::Lo);
            for (const auto& f : faces.stdVector()) {
              const VolIndex& newStartVof = f.getVoF(Side::Lo);

              VofUtils::getVofsInMonotonePath(a_vofSet,
                                              newStartVof,
                                              a_ebisbox,
                                              a_validBox,
                                              a_radius,
                                              newTimesMoved,
                                              newPathSign);
            }
          }

          // Move in the high direction. If we have already moved in the low direction we are not allowed to "turn
          // back".
          if (a_pathSign[dir] >= 0) {
            IntVect newPathSign = a_pathSign;
            newPathSign[dir]    = 1;

            Vector<FaceIndex> faces = a_ebisbox.getFaces(a_startVof, dir, Side::Hi);
            for (const auto& f : faces.stdVector()) {
              const VolIndex& newStartVof = f.getVoF(Side::Hi);

              VofUtils::getVofsInMonotonePath(a_vofSet,
                                              newStartVof,
                                              a_ebisbox,
                                              a_validBox,
                                              a_radius,
                                              newTimesMoved,
                                              newPathSign);
            }
          }
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
