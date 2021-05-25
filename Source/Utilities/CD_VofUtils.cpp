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

Vector<VolIndex> VofUtils::getAllVofsInBox(const Box& a_box, const EBISBox& a_ebisbox){
  Vector<VolIndex> ret;

  Box bx = a_box;
  bx    &= a_ebisbox.getDomain();
  for (BoxIterator bit(bx); bit.ok(); ++bit){
    ret.append(a_ebisbox.getVoFs(bit()));
  }

  return ret;
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

										       
Vector<VolIndex> VofUtils::getConnectedVofs(const VolIndex& a_startVof, const Vector<VolIndex>& a_allVofs, const EBISBox& a_ebisbox){

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

Vector<VolIndex> VofUtils::getAllConnectedVofsInRadius(const VolIndex&   a_startVof,
						       const EBISBox&    a_ebisbox,
						       const int         a_radius,
						       const IntVectSet& a_excludeIVS){
  Vector<VolIndex> allVofs       = VofUtils::getAllVofsInRadius(a_startVof, a_ebisbox, a_radius, true);
  Vector<VolIndex> connectedVofs = VofUtils::getConnectedVofs(a_startVof, allVofs, a_ebisbox);

  VofUtils::excludeCells(connectedVofs, a_excludeIVS);

  return connectedVofs;
}

Vector<VolIndex> VofUtils::getAllConnectedVofsInBox(const VolIndex&   a_startVof,
						    const EBISBox&    a_ebisbox,
						    const Box&        a_box,
						    const IntVectSet& a_excludeIVS){

  // Get all Vofs in box and exclude them so they don't include cells in a_excludeIVS
  Vector<VolIndex> allVofs       = VofUtils::getAllVofsInBox(a_box, a_ebisbox);
  Vector<VolIndex> connectedVofs = VofUtils::getConnectedVofs(a_startVof, allVofs, a_ebisbox);

  VofUtils::excludeCells(connectedVofs, a_excludeIVS);

  return connectedVofs;  
}

Vector<VolIndex> VofUtils::getAllVofsInQuadrant(const VolIndex& a_startVof, const EBISBox& a_ebisbox, const RealVect& a_normal, const int a_radius, const bool a_addStartVof) {

  Vector<VolIndex> ret;
  
  if(VofUtils::isQuadrantWellDefined(a_normal)){

    // Define the quadrant
    const IntVect quadrant = VofUtils::getQuadrant(a_normal);

    ret = VofUtils::getAllVofsInQuadrant(a_startVof, a_ebisbox, quadrant, a_radius, a_addStartVof);
  }

  return ret;
}

Vector<VolIndex> VofUtils::getAllVofsInQuadrant(const VolIndex& a_startVof, const EBISBox& a_ebisbox, const IntVect a_quadrant, const int a_radius, const bool a_addStartVof) {

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


  ret = VofUtils::getAllConnectedVofsInBox(a_startVof, a_ebisbox, bx, IntVectSet());

  if(a_addStartVof){
    ret.push_back(a_startVof);
  }

  return ret;
}

Vector<VolIndex> VofUtils::getAllVofsSymmetric(const VolIndex&      a_startVof,
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

  ret = VofUtils::getAllConnectedVofsInBox(a_startVof, a_ebisbox, bx, IntVectSet());

  if(a_addStartVof){
    ret.push_back(a_startVof);
  }

  return ret;
}

#include <CD_NamespaceFooter.H>
