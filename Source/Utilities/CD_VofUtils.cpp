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
#include <CD_NamespaceHeader.H>

Vector<VolIndex> VofUtils::getVofsInRadius(const VolIndex&    a_startVof,
					   const EBISBox&     a_ebisbox,
					   const int          a_radius,
					   const Connectivity a_connectivity,
					   const bool         a_addStartVof){
  Vector<VolIndex> vofs;
  switch(a_connectivity){
  case VofUtils::Connectivity::MonotonePath:
    vofs = VofUtils::getVofsInMonotonePath(a_startVof, a_ebisbox, a_radius, a_addStartVof);
    break;
  case VofUtils::Connectivity::SimplyConnected:
    vofs = VofUtils::getConnectedVofsInRadius(a_startVof, a_ebisbox, a_radius, a_addStartVof);
    break;
  case VofUtils::Connectivity::All:
    vofs = VofUtils::getAllVofsInRadius(a_startVof, a_ebisbox, a_radius, a_addStartVof);
    break;
  }

  return vofs;
}

Vector<VolIndex> VofUtils::getVofsInMonotonePath(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVof){
  Vector<VolIndex> ret;
  Vector<VolIndex> vofList;

  IntVect timesMoved = IntVect::Zero;
  IntVect pathSign   = IntVect::Zero;

  VofUtils::getVofsInMonotonePath(vofList, a_startVoF, a_ebisbox, a_radius, timesMoved, pathSign);

  if(!a_addStartVof){
    for (int i = 0; i < vofList.size(); i++){
      if(vofList[i] != a_startVoF) ret.push_back(vofList[i]);
    }
  }
  else{
    ret = vofList;
  }

  return ret;
}

Vector<VolIndex> VofUtils::getConnectedVofsInRadius(const VolIndex& a_startVof, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVof){
  MayDay::Abort("VofUtils::getConnectedVofsInRadius - not implemented");
}

Vector<VolIndex> VofUtils::getAllVofsInRadius(const VolIndex& a_startVof, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVof){
  Vector<VolIndex> allVofs;
  Box bx(a_startVof.gridIndex(), a_startVof.gridIndex());
  bx.grow(a_radius);
  bx &= a_ebisbox.getDomain();
  for (BoxIterator bit(bx); bit.ok(); ++bit){

    const std::vector<VolIndex> vofsInThisCell = a_ebisbox.getVoFs(bit()).stdVector();

    for (const auto& ivof : vofsInThisCell){
      if(ivof != a_startVof){
	allVofs.push_back(ivof);
      }
    }
  }

  if(a_addStartVof){
    allVofs.push_back(a_startVof);
  }

  return allVofs;
}

void VofUtils::excludeCells(Vector<VolIndex>& a_vofs, const Box& a_excludeBox){

  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()){
    if(!a_excludeBox.contains(ivof.gridIndex())){
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void VofUtils::excludeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_excludeIVS){

  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()){
    if(!a_excludeIVS.contains(ivof.gridIndex())){
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void VofUtils::includeCells(Vector<VolIndex>& a_vofs, const Box& a_includeBox){
  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()){
    if(a_includeBox.contains(ivof.gridIndex())){
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}

void VofUtils::includeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_includeIVS){
  Vector<VolIndex> newVofs;

  for (auto& ivof : a_vofs.stdVector()){
    if(a_includeIVS.contains(ivof.gridIndex())){
      newVofs.push_back(ivof);
    }
  }

  a_vofs = newVofs;
}


// STUFF BELOW HERE SHOULD BE PROTECTED!
bool VofUtils::isQuadrantWellDefined(const RealVect a_normal){
  bool ret = true;

  const RealVect v = a_normal;///a_normal.vectorLength(); // Idiot guard. 
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if (std::abs(a_normal[dir]) == 1.0 || std::abs(a_normal[dir]) == 0.0) {
      ret = false;
    }
  }

  return ret;
}

std::pair<int, Side::LoHiSide> VofUtils::getCardinalDirection(const RealVect a_normal){ 
  std::pair<int, Side::LoHiSide> ret;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if (a_normal[dir] == 1.0) {
      ret = std::make_pair(dir, Side::Hi);
    }
    else if(a_normal[dir] == -1.0){
      ret = std::make_pair(dir, Side::Lo);
    }
  }

  return ret;
}

IntVect VofUtils::getQuadrant(const RealVect a_normal){

  IntVect quadrant = IntVect::Zero;

  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_normal[dir] < 0.0){
      quadrant[dir] = -1;
    }
    else if (a_normal[dir] > 0.0 ) {
      quadrant[dir] = 1;
    }
  }

  return quadrant;
}



Vector<VolIndex> VofUtils::getAllVofsInBox(const Box& a_box, const EBISBox& a_ebisbox){
  Vector<VolIndex> ret;

  Box bx = a_box;
  bx    &= a_ebisbox.getDomain();
  for (BoxIterator bit(bx); bit.ok(); ++bit){
    ret.append(a_ebisbox.getVoFs(bit()));
  }

  return ret;
}



										       
Vector<VolIndex> VofUtils::connectedVofsOnly(const VolIndex& a_startVof, const Vector<VolIndex>& a_allVofs, const EBISBox& a_ebisbox){

  // If a_allVofs does not contain a_startVof, add it.
  bool hasStartingVof = false;
  std::vector<VolIndex> includeVofs;
  for (int i = 0; i < a_allVofs.size(); i++){
    const VolIndex& ivof = a_allVofs[i];
    
    if(ivof == a_startVof){
      hasStartingVof = true;
    }

    includeVofs.emplace_back(ivof);
  }
  if(!hasStartingVof){
    includeVofs.push_back(a_startVof);
  }


  // Make a map over which Vofs we've added so far
  std::map<VolIndex, bool> vofMap;
  for (const auto& ivof : includeVofs){
    vofMap.emplace(ivof, false);
  }

  // Insert the first vof
  std::vector<VolIndex> connectedVofs;
  connectedVofs.emplace_back(a_startVof);
  vofMap.at(a_startVof) = true;


  // Add Vofs that are connected to this vof unless they have not already been added.
  // We exit when we can't add more Vofs. 
  bool keepGoing = true;
  while(keepGoing){
    keepGoing = false;

    std::vector<VolIndex> newVofs = connectedVofs;
    
    for (const auto& ivof : connectedVofs){

      for (const auto& m : vofMap){
	const VolIndex otherVof = m.first;
	
	const bool isConnected  = a_ebisbox.isConnected(ivof, otherVof);
	const bool alreadyAdded = vofMap.at(otherVof);

	if(!alreadyAdded && isConnected){
	  newVofs.emplace_back(otherVof);
	  vofMap.at(otherVof) = true;

	  keepGoing = true;
	}
      }
    }

    connectedVofs = newVofs;
  }


  Vector<VolIndex> ret;
  for (const auto& ivof : connectedVofs){
    if(ivof != a_startVof){
      ret.push_back(ivof);
    }
  }

  return ret;

}

Vector<VolIndex> VofUtils::getConnectedVofsInRadius(const VolIndex&   a_startVof,
						    const EBISBox&    a_ebisbox,
						    const int         a_radius,
						    const IntVectSet& a_excludeIVS){
  Vector<VolIndex> allVofs       = VofUtils::getAllVofsInRadius(a_startVof, a_ebisbox, a_radius, true);
  Vector<VolIndex> connectedVofs = VofUtils::connectedVofsOnly(a_startVof, allVofs, a_ebisbox);

  VofUtils::excludeCells(connectedVofs, a_excludeIVS);

  return connectedVofs;
}

Vector<VolIndex> VofUtils::getConnectedVofsInBox(const VolIndex&   a_startVof,
						 const EBISBox&    a_ebisbox,
						 const Box&        a_box,
						 const IntVectSet& a_excludeIVS){

  // Get all Vofs in box and exclude them so they don't include cells in a_excludeIVS
  Vector<VolIndex> allVofs       = VofUtils::getAllVofsInBox(a_box, a_ebisbox);
  Vector<VolIndex> connectedVofs = VofUtils::connectedVofsOnly(a_startVof, allVofs, a_ebisbox);

  VofUtils::excludeCells(connectedVofs, a_excludeIVS);

  return connectedVofs;  
}

Vector<VolIndex> VofUtils::getConnectedVofsInQuadrant(const VolIndex& a_startVof,
						      const EBISBox&  a_ebisbox,
						      const RealVect& a_normal,
						      const int       a_radius,
						      const bool      a_addStartVof) {

  Vector<VolIndex> ret;
  
  if(VofUtils::isQuadrantWellDefined(a_normal)){

    // Define the quadrant
    const IntVect quadrant = VofUtils::getQuadrant(a_normal);

    ret = VofUtils::getConnectedVofsInQuadrant(a_startVof, a_ebisbox, quadrant, a_radius, a_addStartVof);
  }

  return ret;
}

Vector<VolIndex> VofUtils::getConnectedVofsInQuadrant(const VolIndex& a_startVof,
						      const EBISBox&  a_ebisbox,
						      const IntVect   a_quadrant,
						      const int       a_radius,
						      const bool      a_addStartVof) {

  Vector<VolIndex> ret;
  
  // Grow the box in the direction defined by the quadrant
  Box bx(a_startVof.gridIndex(), a_startVof.gridIndex());
  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_quadrant[dir] > 0){
      bx.growHi(dir, a_radius); 
    }
    else if(a_quadrant[dir] < 0){
      bx.growLo(dir, a_radius);
    }
  }
  bx &= a_ebisbox.getDomain().domainBox();


  ret = VofUtils::getConnectedVofsInBox(a_startVof, a_ebisbox, bx, IntVectSet());

  if(a_addStartVof){
    ret.push_back(a_startVof);
  }

  return ret;
}

Vector<VolIndex> VofUtils::getConnectedVofsSymmetric(const VolIndex&      a_startVof,
						     const EBISBox&       a_ebisbox,
						     const int            a_cardinal,
						     const Side::LoHiSide a_side,
						     const int            a_radius,
						     const bool           a_addStartVof){
  Vector<VolIndex> ret;

  // Create candiate box where we allow ourself to fetch Vofs. 
  Box bx(a_startVof.gridIndex(), a_startVof.gridIndex());
  for (int dir = 0; dir < SpaceDim; dir++){
    if(dir == a_cardinal){
      bx.growDir(dir, a_side, a_radius);
    }
    else{
      bx.growLo(dir, a_radius);
      bx.growHi(dir, a_radius);
    }
  }
  bx &= a_ebisbox.getDomain().domainBox();

  ret = VofUtils::getConnectedVofsInBox(a_startVof, a_ebisbox, bx, IntVectSet());

  if(a_addStartVof){
    ret.push_back(a_startVof);
  }

  return ret;
}

void VofUtils::getVofsInMonotonePath(Vector<VolIndex>& a_vofList,
				     const VolIndex&   a_startVof,
				     const EBISBox&    a_ebisbox,
				     const int         a_radius,
				     const IntVect&    a_timesMoved,
				     const IntVect&    a_pathSign){

  const ProblemDomain& a_domain = a_ebisbox.getDomain();

  if(a_domain.contains(a_startVof.gridIndex())){

    // Add if not already added
    bool haveStartVof = false;
    for (int ivof = 0; ivof < a_vofList.size(); ivof++){
      if(a_vofList[ivof] == a_startVof) haveStartVof = true;
    }

    if(!haveStartVof) a_vofList.push_back(a_startVof);

    for (int dir = 0; dir < SpaceDim; dir++){
      if(a_timesMoved[dir] < a_radius){
	const IntVect newTimesMoved = a_timesMoved + BASISV(dir);

	Side::LoHiSide whichSide;
	
	// Move in the low direction. If we have already moved in the high direction we are not allowed to "turn back".
	if(a_pathSign[dir] == -1 || a_pathSign[dir] == 0){
	  IntVect newPathSign = a_pathSign;
	  newPathSign[dir] = -1;

	  Vector<FaceIndex> faces = a_ebisbox.getFaces(a_startVof, dir, Side::Lo);
	  for (const auto& f : faces.stdVector()){
	    const VolIndex& newStartVof = f.getVoF(Side::Lo);

	    VofUtils::getVofsInMonotonePath(a_vofList, newStartVof, a_ebisbox, a_radius, newTimesMoved, newPathSign);
	  }
	}

	// Move in the high direction. If we have already moved in the low direction we are not allowed to "turn back".
	if(a_pathSign[dir] == 0 || a_pathSign[dir] == 1){
	  IntVect newPathSign = a_pathSign;
	  newPathSign[dir] = 1;

	  Vector<FaceIndex> faces = a_ebisbox.getFaces(a_startVof, dir, Side::Hi);
	  for (const auto& f : faces.stdVector()){
	    const VolIndex& newStartVof = f.getVoF(Side::Hi);	    

	    VofUtils::getVofsInMonotonePath(a_vofList, newStartVof, a_ebisbox, a_radius, newTimesMoved, newPathSign);
	  }
	}
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
