/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_VofUtils.cpp
  @brief  Implementation of CD_VofUtils.H
  @author Robert Marskar
*/

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

  // Can only do this if actually have a normal vecto r
  if (a_normal != RealVect::Zero) {

    // Fetch vofs using radius.
    vofs = VofUtils::getVofsInRadius(a_startVof, a_ebisbox, a_radius, a_connectivity, a_addStartVof);

    // Find the quadrant or "half-plane" where the quadrant vofs live, and restrict to that quadrant
    Box quadBox;
    if (VofUtils::isQuadrantWellDefined(a_normal)) {
      // Box is the quadrant box
      quadBox = VofUtils::getQuadrant(a_normal, a_startVof, a_ebisbox, a_radius);
    }
    else {
      // Box consists of two quadrants.
      const std::pair<int, Side::LoHiSide> cardinal = VofUtils::getCardinalDirection(a_normal);
      quadBox = VofUtils::getSymmetricQuadrant(cardinal, a_startVof, a_ebisbox, a_radius);
    }

    VofUtils::includeCells(vofs, quadBox);
  }

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
    const RealVect x0 = Location::position(a_vofLocation, a_startVof, a_ebisbox, 1.0);

    const Vector<VolIndex> vofsInRadius = VofUtils::getVofsInRadius(a_startVof,
                                                                    a_ebisbox,
                                                                    a_radius,
                                                                    a_connectivity,
                                                                    false);

    for (int i = 0; i < vofsInRadius.size(); i++) {
      const VolIndex& curVoF = vofsInRadius[i];

      const RealVect delta = Location::position(a_cellLocation, curVoF, a_ebisbox, 1.0) - x0;

      if (delta.dotProduct(a_normal) > a_deltaThresh) {
        vofs.push_back(curVoF);
      }
    }

    if (a_addStartVof) {
      vofs.push_back(a_startVof);
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

  auto kernel = [&](const IntVect& iv) -> void {
    const std::vector<VolIndex> vofsInThisCell = a_ebisbox.getVoFs(iv).stdVector();

    for (const auto& ivof : vofsInThisCell) {
      if (ivof != a_startVof) {
        allVofs.push_back(ivof);
      }
    }
  };

  BoxLoops::loop(bx, kernel);

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

bool
VofUtils::isQuadrantWellDefined(const RealVect a_normal)
{
  bool ret = true;

  const RealVect v = a_normal; ///a_normal.vectorLength(); // Idiot guard.

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (std::abs(a_normal[dir]) == 1.0 || std::abs(a_normal[dir]) == 0.0) {
      ret = false;
    }
  }

  return ret;
}

std::pair<int, Side::LoHiSide>
VofUtils::getCardinalDirection(const RealVect a_normal)
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
      bx.growLo(dir, a_radius);
    }
    else if (a_normal[dir] > 0.0) {
      bx.growHi(dir, a_radius);
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
      bx.growDir(dir, a_cardinal.second, a_radius);
    }
    else {
      bx.growLo(dir, a_radius);
      bx.growHi(dir, a_radius);
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
      if (a_vofList[ivof] == a_startVof)
        haveStartVof = true;
    }

    if (!haveStartVof) {
      a_vofList.push_back(a_startVof);

      for (int dir = 0; dir < SpaceDim; dir++) {
        if (a_timesMoved[dir] < a_radius) {
          const IntVect newTimesMoved = a_timesMoved + BASISV(dir);

          // Move in the low direction. If we have already moved in the high direction we are not allowed to "turn back".
          if (a_pathSign[dir] == -1 || a_pathSign[dir] == 0) {
            IntVect newPathSign = a_pathSign;
            newPathSign[dir]    = -1;

            Vector<FaceIndex> faces = a_ebisbox.getFaces(a_startVof, dir, Side::Lo);
            for (const auto& f : faces.stdVector()) {
              const VolIndex& newStartVof = f.getVoF(Side::Lo);

              VofUtils::getVofsInMonotonePath(a_vofList, newStartVof, a_ebisbox, a_radius, newTimesMoved, newPathSign);
            }
          }

          // Move in the high direction. If we have already moved in the low direction we are not allowed to "turn back".
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

          // Move in the low direction. If we have already moved in the high direction we are not allowed to "turn back".
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

          // Move in the high direction. If we have already moved in the low direction we are not allowed to "turn back".
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
