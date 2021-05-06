/*!
  @file   VoFUtils.H
  @brief  Various functions for getting VoFs near cut-cells
  @author Robert Marskar
  @date   May 2021
  @todo   getAllVoFsInQuadrant is hard-coded for 2D right now. 
*/

#include "VoFUtils.H"

#include "EBArith.H"

Vector<VolIndex> VoFUtils::getAllVoFsInQuadrant(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const RealVect a_normal, const bool a_addStartVoF) {
  IntVect quadrant;

  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_normal[dir] < 0){
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

  Vector<VolIndex> ret;
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
