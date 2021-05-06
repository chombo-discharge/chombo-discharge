/*!
  @file   VoFUtils.H
  @brief  Various functions for getting VoFs near cut-cells
  @author Robert Marskar
  @date   May 2021
  @todo   getAllVoFsInQuadrant is hard-coded for 2D right now. 
*/

#include "VoFUtils.H"

#include "LoHiSide.H"

bool VoFUtils::isQuadrantWellDefined(const RealVect a_normal){
  bool ret = true;

  const RealVect v = a_normal;///a_normal.vectorLength(); // Idiot guard. 
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if (std::abs(a_normal[dir]) == 1.0 || std::abs(a_normal[dir]) == 0.0) {
      ret = false;
    }
  }

  return ret;
}

std::pair<int, Side::LoHiSide> VoFUtils::getCardinalDirection(const RealVect a_normal){ 
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

IntVect VoFUtils::getQuadrant(const RealVect a_normal){

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

Vector<VolIndex> VoFUtils::getAllVoFsInRadius(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVoF){
  Vector<VolIndex> allVoFs;
  Box bx(a_startVoF.gridIndex(), a_startVoF.gridIndex());
  bx.grow(a_radius);
  bx &= a_ebisbox.getDomain();
  for (BoxIterator bit(bx); bit.ok(); ++bit){

    const std::vector<VolIndex> vofsInThisCell = a_ebisbox.getVoFs(bit()).stdVector();

    for (const auto& ivof : vofsInThisCell){
      if(ivof != a_startVoF){
	allVoFs.push_back(ivof);
      }
    }
  }

  if(a_addStartVoF){
    allVoFs.push_back(a_startVoF);
  }

  return allVoFs;
}

Vector<VolIndex> VoFUtils::getAllVoFsInBox(const Box& a_box, const EBISBox& a_ebisbox){
  Vector<VolIndex> ret;

  Box bx = a_box;
  bx    &= a_ebisbox.getDomain();
  for (BoxIterator bit(bx); bit.ok(); ++bit){
    ret.append(a_ebisbox.getVoFs(bit()));
  }

  return ret;
}

void VoFUtils::excludeCells(Vector<VolIndex>& a_vofs, const Box& a_excludeBox){

  Vector<VolIndex> newVoFs;

  for (auto& ivof : a_vofs.stdVector()){
    if(!a_excludeBox.contains(ivof.gridIndex())){
      newVoFs.push_back(ivof);
    }
  }

  a_vofs = newVoFs;
}

void VoFUtils::excludeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_excludeIVS){

  Vector<VolIndex> newVoFs;

  for (auto& ivof : a_vofs.stdVector()){
    if(!a_excludeIVS.contains(ivof.gridIndex())){
      newVoFs.push_back(ivof);
    }
  }

  a_vofs = newVoFs;
}

void VoFUtils::includeCells(Vector<VolIndex>& a_vofs, const Box& a_includeBox){
  Vector<VolIndex> newVoFs;

  for (auto& ivof : a_vofs.stdVector()){
    if(a_includeBox.contains(ivof.gridIndex())){
      newVoFs.push_back(ivof);
    }
  }

  a_vofs = newVoFs;
}

void VoFUtils::includeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_includeIVS){
  Vector<VolIndex> newVoFs;

  for (auto& ivof : a_vofs.stdVector()){
    if(a_includeIVS.contains(ivof.gridIndex())){
      newVoFs.push_back(ivof);
    }
  }

  a_vofs = newVoFs;
}

										       
Vector<VolIndex> VoFUtils::getConnectedVoFs(const VolIndex& a_startVoF, const Vector<VolIndex>& a_allVoFs, const EBISBox& a_ebisbox){

  // If a_allVoFs does not contain a_startVoF, add it.
  bool hasStartingVoF = false;
  std::vector<VolIndex> includeVoFs;
  for (int i = 0; i < a_allVoFs.size(); i++){
    const VolIndex& ivof = a_allVoFs[i];
    
    if(ivof == a_startVoF){
      hasStartingVoF = true;
    }

    includeVoFs.emplace_back(ivof);
  }
  if(!hasStartingVoF){
    includeVoFs.push_back(a_startVoF);
  }


  // Make a map over which VoFs we've added so far
  std::map<VolIndex, bool> vofMap;
  for (const auto& ivof : includeVoFs){
    vofMap.emplace(ivof, false);
  }

  // Insert the first vof
  std::vector<VolIndex> connectedVoFs;
  connectedVoFs.emplace_back(a_startVoF);
  vofMap.at(a_startVoF) = true;


  // Add VoFs that are connected to this vof unless they have not already been added.
  // We exit when we can't add more VoFs. 
  bool keepGoing = true;
  while(keepGoing){
    keepGoing = false;

    std::vector<VolIndex> newVoFs = connectedVoFs;
    
    for (const auto& ivof : connectedVoFs){

      for (const auto& m : vofMap){
	const VolIndex otherVof = m.first;
	
	const bool isConnected  = a_ebisbox.isConnected(ivof, otherVof);
	const bool alreadyAdded = vofMap.at(otherVof);

	if(!alreadyAdded && isConnected){
	  newVoFs.emplace_back(otherVof);
	  vofMap.at(otherVof) = true;

	  keepGoing = true;
	}
      }
    }

    connectedVoFs = newVoFs;
  }


  Vector<VolIndex> ret;
  for (const auto& ivof : connectedVoFs){
    if(ivof != a_startVoF){
      ret.push_back(ivof);
    }
  }

  return ret;

}

Vector<VolIndex> VoFUtils::getAllConnectedVoFsInRadius(const VolIndex&   a_startVoF,
						       const EBISBox&    a_ebisbox,
						       const int         a_radius,
						       const IntVectSet& a_excludeIVS){
  Vector<VolIndex> allVoFs       = VoFUtils::getAllVoFsInRadius(a_startVoF, a_ebisbox, a_radius, true);
  Vector<VolIndex> connectedVoFs = VoFUtils::getConnectedVoFs(a_startVoF, allVoFs, a_ebisbox);

  VoFUtils::excludeCells(connectedVoFs, a_excludeIVS);

  return connectedVoFs;
}

Vector<VolIndex> VoFUtils::getAllConnectedVoFsInBox(const VolIndex&   a_startVoF,
						    const EBISBox&    a_ebisbox,
						    const Box&        a_box,
						    const IntVectSet& a_excludeIVS){

  // Get all Vofs in box and exclude them so they don't include cells in a_excludeIVS
  Vector<VolIndex> allVoFs       = VoFUtils::getAllVoFsInBox(a_box, a_ebisbox);
  Vector<VolIndex> connectedVoFs = VoFUtils::getConnectedVoFs(a_startVoF, allVoFs, a_ebisbox);

  VoFUtils::excludeCells(connectedVoFs, a_excludeIVS);

  return connectedVoFs;  
}

Vector<VolIndex> VoFUtils::getAllVoFsInQuadrant(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const RealVect& a_normal, const int a_radius, const bool a_addStartVoF) {

  Vector<VolIndex> ret;
  
  if(VoFUtils::isQuadrantWellDefined(a_normal)){

    // Define the quadrant
    const IntVect quadrant = VoFUtils::getQuadrant(a_normal);

    // Grow the box in the direction defined by the quadrant
    Box bx(a_startVoF.gridIndex(), a_startVoF.gridIndex());
    for (int dir = 0; dir < SpaceDim; dir++){
      if(quadrant[dir] > 0){
	bx.growHi(dir, a_radius); 
      }
      else if(quadrant[dir] < 0){
	bx.growLo(dir, a_radius);
      }
      else{
	MayDay::Abort("VoFUtils::getAllVoFsInQuadrant - logic bust");
      }
    }
    bx &= a_ebisbox.getDomain().domainBox();


    ret = VoFUtils::getAllConnectedVoFsInBox(a_startVoF, a_ebisbox, bx, IntVectSet());

    if(a_addStartVoF){
      ret.push_back(a_startVoF);
    }

  }

  return ret;
}

Vector<VolIndex> VoFUtils::getAllVoFsSymmetric(const VolIndex&      a_startVoF,
					       const EBISBox&       a_ebisbox,
					       const int            a_cardinal,
					       const Side::LoHiSide a_side,
					       const int            a_radius,
					       const bool           a_addStartVoF){
  Vector<VolIndex> ret;

  // Create candiate box where we allow ourself to fetch VoFs. 
  Box bx(a_startVoF.gridIndex(), a_startVoF.gridIndex());
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

  ret = VoFUtils::getAllConnectedVoFsInBox(a_startVoF, a_ebisbox, bx, IntVectSet());

  if(a_addStartVoF){
    ret.push_back(a_startVoF);
  }

  return ret;
}
