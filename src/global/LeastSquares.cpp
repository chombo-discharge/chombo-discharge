/*!
  @file   LeastSquares.cpp
  @brief  Implementation of LeastSquares.H
  @author Robert Marskar
  @date   May 2016
*/

#include "LeastSquares.H"
#include "LaPackUtils.H"
#include "VoFUtils.H"

#include "EBArith.H"



RealVect LeastSquares::position(const CellPosition a_position,
				const VolIndex&    a_vof,
				const EBISBox&     a_ebisbox,
				const Real&        a_dx){

  RealVect ret;
  switch(a_position){
  case CellPosition::Center:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit;
    break;
  case CellPosition::Centroid:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit + a_ebisbox.centroid(a_vof);
    break;
  case CellPosition::Boundary:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit + a_ebisbox.bndryCentroid(a_vof);
    break;
  }

  ret *= a_dx;

  return ret;
}

RealVect LeastSquares::displacement(const CellPosition a_from,
				    const CellPosition a_to,
				    const VolIndex&    a_fromVoF,
				    const VolIndex&    a_toVoF,
				    const EBISBox&     a_ebisbox,
				    const Real&        a_dx){

  const RealVect a = LeastSquares::position(a_from, a_fromVoF, a_ebisbox, a_dx);
  const RealVect b = LeastSquares::position(a_to,   a_toVoF,   a_ebisbox, a_dx);

  return (b-a);
}

RealVect LeastSquares::displacement(const CellPosition a_from,
				    const CellPosition a_to,
				    const VolIndex&    a_fromVoF,
				    const VolIndex&    a_toVoF,
				    const EBISBox&     a_ebisboxFrom,
				    const EBISBox&     a_ebisboxTo,
				    const Real&        a_dxFrom,
				    const Real&        a_dxTo){

  const RealVect a = LeastSquares::position(a_from, a_fromVoF, a_ebisboxFrom, a_dxFrom);
  const RealVect b = LeastSquares::position(a_to,   a_toVoF,   a_ebisboxTo,   a_dxTo);

  return (b-a);
}

Vector<RealVect> LeastSquares::getDisplacements(const CellPosition      a_from,
						const CellPosition      a_to,
						const VolIndex&         a_fromVoF,
						const Vector<VolIndex>& a_toVoFs,
						const EBISBox&          a_ebisbox,
						const Real&             a_dx){

  Vector<RealVect> ret;

  for (int i = 0; i < a_toVoFs.size(); i++){
    const RealVect d = LeastSquares::displacement(a_from, a_to, a_fromVoF, a_toVoFs[i], a_ebisbox, a_dx);

    ret.push_back(d);
  }

  return ret;
}

Vector<Real> LeastSquares::makeDiagWeights(const Vector<RealVect>& a_displacements, const int a_pow){
  Vector<Real> ret(a_displacements.size(), 1.0);
  
  if(a_pow > 0){
    for (int i = 0; i < a_displacements.size(); i++){
      const RealVect& d = a_displacements[i];
      const Real len    = d.vectorLength();
      const Real w      = 1./std::pow(len, a_pow);

      ret[i] = w;
    }
  }

  return ret;
}

void LeastSquares::removeEquations(Vector<VolIndex>& a_allVoFs, Vector<RealVect>& a_displacements, const Real a_tolerance){

  Vector<VolIndex> newVoFs;
  Vector<RealVect> newDeltas;

  for (int i = 0; i < a_allVoFs.size(); i++){
    if(a_displacements[i].vectorLength() > a_tolerance){
      newVoFs.push_back(a_allVoFs[i]);
      newDeltas.push_back(a_displacements[i]);
    }
  }

  a_allVoFs       = newVoFs;
  a_displacements = newDeltas;
}

VoFStencil LeastSquares::getGradStenOrderOne(const Vector<VolIndex>& a_allVoFs,
					     const Vector<RealVect>& a_displacements,
					     const Vector<Real>&     a_weights){

  VoFStencil sten;

  const int D = SpaceDim;
  const int K = a_displacements.size();

  // Make sure we have enough equations. 
  if(K < D) MayDay::Abort("LeastSquares::getGradStenOrderOne -- not enough equations, must have K >= 3");
  
  // Build A, which is a KxD matrix. 
  Vector<Real> linA(K*D, 0.0);
  for (int k = 0; k < K; k++){ // Row
    for (int d = 0; d < D; d++){ // Column
      const int idx = LaPackUtils::linearIndex(k, d, K, D);
      
      const Real delta = a_displacements[k][d];
      
      linA[idx] = a_weights[k]*delta;
    }
  }

  // Compute pseudoinverse
  Vector<Real> linAplus(D*K, 0.0);
  const bool foundSVD = LaPackUtils::computePseudoInverse(linAplus.stdVector(), linA.stdVector(), K, D);

  if(foundSVD){
    
    sten.clear();

    // Pseudoinverse is a DxK matrix. 
    for (int dir = 0; dir < D; dir++){ // Row
      for (int k = 0; k < K; k++){ // Column
	
	const int idx       = LaPackUtils::linearIndex(dir, k, D, K);
	
	const Real weight   = a_weights[k]*linAplus[idx];
	const VolIndex& vof = a_allVoFs[k];        
	
	sten.add(vof, weight, dir);
      }
    }
  }

  return sten;

}

VoFStencil LeastSquares::getGradStenOrderOne(const Vector<VolIndex>& a_allVoFs,
					     const Vector<RealVect>& a_displacements,
					     const int&              a_p){

  // Build weights and call the other version. 
  Vector<Real> weights = LeastSquares::makeDiagWeights(a_displacements, a_p);

  return LeastSquares::getGradStenOrderOne(a_allVoFs, a_displacements, weights);
}

VoFStencil LeastSquares::getBndryGradStenOrderOne(const VolIndex& a_vof,
						  const EBISBox&  a_ebisbox,
						  const Real&     a_dx,
						  const int       a_p){
  VoFStencil bndrySten;

  const bool addStartingVoF = false;
  const RealVect normal     = a_ebisbox.normal(a_vof);

  if(normal != RealVect::Zero){ // Don't support zero normal vector (yet)
    const int radius = 1;
    const int order  = 1;
    const int numTaylorTerms = LeastSquares::getTaylorExpansionSize(order);

    Vector<VolIndex> allVoFs;
    if(VoFUtils::isQuadrantWellDefined(normal)){ // Try to use quadrants. 
      allVoFs = VoFUtils::getAllVoFsInQuadrant(a_vof, a_ebisbox, normal, radius, addStartingVoF);
    }
    else{
      const std::pair<int, Side::LoHiSide> cardinal = VoFUtils::getCardinalDirection(normal);

      allVoFs = VoFUtils::getAllVoFsSymmetric(a_vof, a_ebisbox, cardinal.first, cardinal.second, radius, addStartingVoF);
    }

    // Build the stencil if we can. 
    if(allVoFs.size() >= numTaylorTerms){
      const Vector<RealVect> displacements = LeastSquares::getDisplacements(CellPosition::Boundary,
									    CellPosition::Center,
									    a_vof,
									    allVoFs,
									    a_ebisbox,
									    a_dx);
    
      bndrySten = LeastSquares::getGradStenOrderOne(allVoFs, displacements, a_p);
    }
  }
  else{
    MayDay::Warning("LeastSquares::getBndryGradStenOrderOne -- got zero normal vector but not sure what to do about that. You probably shouldn't use this routine.");
  }

  return bndrySten;
}

VoFStencil LeastSquares::projectGradSten(const VoFStencil& a_stencil, const RealVect& a_projection) {
  VoFStencil sten;

  for (int i = 0; i < a_stencil.size(); i++){
    const VolIndex& vof = a_stencil.vof(i);
    const Real& weight  = a_stencil.weight(i);
    const int dir       = a_stencil.variable(i);

    const Real p = a_projection[dir];

    sten.add(vof, p*weight);
  }

  return sten;
}

Real LeastSquares::sumWeights(const VoFStencil& a_stencil, const int a_variable){

  Real ret = 0.0;

  for (int i = 0; i < a_stencil.size(); i++){
    const int var = a_stencil.variable(i);
    if(var == a_variable){
      ret += a_stencil.weight(i);
    }
  }

  return ret;
}

Real LeastSquares::sumAllWeights(const VoFStencil& a_stencil){
  Real ret = 0.0;

  for (int i = 0; i < a_stencil.size(); i++){
    ret += a_stencil.weight(i);
  }

  return ret;
}

int LeastSquares::getTaylorExpansionSize(const int a_Q){

  int nTerms = 0;
  for (MultiIndex cur(IntVect::Zero); cur <= a_Q; cur.next(a_Q)){
    nTerms ++;
  }
  
  return nTerms;
}

Vector<MultiIndex> LeastSquares::getMultiIndicesLexiOrder(const int a_Q){
  Vector<MultiIndex> ret;
  for (MultiIndex cur = IntVect::Zero; cur <= a_Q; cur.next(a_Q)){
    ret.push_back(cur);
  }

  return ret;
}
