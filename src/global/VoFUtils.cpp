/*!
  @file   VoFUtils.H
  @brief  Various functions for getting VoFs near cut-cells
  @author Robert Marskar
  @date   May 2021
  @todo   getAllVoFsInQuadrant is hard-coded for 2D right now. 
*/

#include "VoFUtils.H"

#include "EBArith.H"

bool VoFUtils::isQuadrantWellDefined(const RealVect a_normal){
  bool ret = true;

  const RealVect v = a_normal/a_normal.vectorLength(); // Idiot guard. 
  
  for (int dir = 0; dir < SpaceDim; dir++){
    if (v[dir] == 1.0 || v[dir] == 0.0) ret = false;
  }

  return ret;
}

IntVect VoFUtils::getQuadrant(const RealVect a_normal){

  IntVect quadrant = IntVect::Zero;

  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_normal[dir] < 0.0){
      quadrant[dir] = -1;
    }
    else if(a_normal[dir] > 0.0){
      quadrant[dir] = 1;
    }
  }

  return quadrant;
}

Vector<VolIndex> VoFUtils::getAllVoFsInQuadrant(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const RealVect a_normal, const int a_radius, const bool a_addStartVoF) {

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

    // Get all VoFs that can be reached with a monotone path and add them if they fall within the quadrant. 
    Vector<VolIndex> monoVoFs;
    EBArith::getAllVoFsInMonotonePath(monoVoFs, a_startVoF, a_ebisbox, a_radius);
    for (auto& ivof : monoVoFs.stdVector()){
      if(bx.contains(ivof.gridIndex()) && ivof != a_startVoF){
	ret.push_back(ivof);
      }
    }
    
    if(a_addStartVoF){
      ret.push_back(a_startVoF);
    }

    return ret;
  }

  

  IntVect quadrant;
  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_normal[dir] <= 0.){
      quadrant[dir] = -1;
    }
    else{
      quadrant[dir] = 1;
    }
  }


  
  IntVect iv0 = a_startVoF.gridIndex();
  IntVect iv1 = iv0 + BASISV(0)*quadrant[0]                         ;
  IntVect iv2 = iv0                         + BASISV(1)*quadrant[1] ;
  IntVect iv3 = iv0 + BASISV(0)*quadrant[0] + BASISV(1)*quadrant[1] ;

  VolIndex vof1, vof2, vof3;

  if(a_ebisbox.getVoFs(iv1).size() > 0) ret.push_back(VolIndex(iv1, 0));
  if(a_ebisbox.getVoFs(iv2).size() > 0) ret.push_back(VolIndex(iv2, 0));
  if(a_ebisbox.getVoFs(iv3).size() > 0) ret.push_back(VolIndex(iv3, 0));

  if(a_addStartVoF){
    ret.push_back(a_startVoF);
  }
  
  return ret;
}

Vector<VolIndex> VoFUtils::getAllVoFsInRadius(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVoF){

  Vector<VolIndex> ret;
  Vector<IntVect>  ivs;
  
  EBArith::getAllVoFsWithinRadiusExcludingStartVoF(ret, ivs, a_startVoF, a_ebisbox, a_radius);

  if(a_addStartVoF){
    ret.push_back(a_startVoF);
  }

  return ret;
}
