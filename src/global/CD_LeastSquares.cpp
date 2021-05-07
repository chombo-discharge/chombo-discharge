/*!
  @file   CD_LeastSquares.cpp
  @brief  Implementation of CD_LeastSquares.H
  @author Robert Marskar
  @date   May 2016
*/

#include "CD_LaPackUtils.H"
#include "CD_LeastSquares.H"
#include "CD_VofUtils.H"
#include "CD_MultiIndex.H"

#include "EBArith.H"

namespace ChomboDischarge{

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
				      const VolIndex&    a_fromVof,
				      const VolIndex&    a_toVof,
				      const EBISBox&     a_ebisbox,
				      const Real&        a_dx){

    const RealVect a = LeastSquares::position(a_from, a_fromVof, a_ebisbox, a_dx);
    const RealVect b = LeastSquares::position(a_to,   a_toVof,   a_ebisbox, a_dx);

    return (b-a);
  }

  RealVect LeastSquares::displacement(const CellPosition a_from,
				      const CellPosition a_to,
				      const VolIndex&    a_fromVof,
				      const VolIndex&    a_toVof,
				      const EBISBox&     a_ebisboxFrom,
				      const EBISBox&     a_ebisboxTo,
				      const Real&        a_dxFrom,
				      const Real&        a_dxTo){

    const RealVect a = LeastSquares::position(a_from, a_fromVof, a_ebisboxFrom, a_dxFrom);
    const RealVect b = LeastSquares::position(a_to,   a_toVof,   a_ebisboxTo,   a_dxTo);

    return (b-a);
  }

  Vector<RealVect> LeastSquares::getDisplacements(const CellPosition      a_from,
						  const CellPosition      a_to,
						  const VolIndex&         a_fromVof,
						  const Vector<VolIndex>& a_toVofs,
						  const EBISBox&          a_ebisbox,
						  const Real&             a_dx){

    Vector<RealVect> ret;

    for (int i = 0; i < a_toVofs.size(); i++){
      const RealVect d = LeastSquares::displacement(a_from, a_to, a_fromVof, a_toVofs[i], a_ebisbox, a_dx);

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

  void LeastSquares::removeEquations(Vector<VolIndex>& a_allVofs, Vector<RealVect>& a_displacements, const Real a_tolerance){

    Vector<VolIndex> newVofs;
    Vector<RealVect> newDeltas;

    for (int i = 0; i < a_allVofs.size(); i++){
      if(a_displacements[i].vectorLength() > a_tolerance){
	newVofs.push_back(a_allVofs[i]);
	newDeltas.push_back(a_displacements[i]);
      }
    }

    a_allVofs       = newVofs;
    a_displacements = newDeltas;
  }

  VoFStencil LeastSquares::computeGradStenOrderOne(const Vector<VolIndex>& a_allVofs,
						   const Vector<RealVect>& a_displacements,
						   const int&              a_p){
    Vector<Real> weights = LeastSquares::makeDiagWeights(a_displacements, a_p);

    return LeastSquares::computeGradStenOrderOne(a_allVofs, a_displacements, weights);
  }

  VoFStencil LeastSquares::computeGradStenOrderOne(const Vector<VolIndex>& a_allVofs,
						   const Vector<RealVect>& a_displacements,
						   const Vector<Real>&     a_weights){

    VoFStencil sten;

    const int D = SpaceDim;
    const int K = a_displacements.size();

    // Make sure we have enough equations. 
    if(K < D) MayDay::Abort("LeastSquares::computeGradStenOrderOne -- not enough equations, must have K >= 3");
  
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
	  const VolIndex& vof = a_allVofs[k];        
	
	  sten.add(vof, weight, dir);
	}
      }
    }

    return sten;

  }

  VoFStencil LeastSquares::getBndryGradStenOrderOne(const VolIndex& a_vof,
						    const EBISBox&  a_ebisbox,
						    const Real&     a_dx,
						    const int       a_p){
    VoFStencil bndrySten;

    const bool addStartingVof = false;
    const RealVect normal     = a_ebisbox.normal(a_vof);

    if(normal != RealVect::Zero){ // Can only do this if we actual have a normal vector. 
      const int radius = 1;
      const int order  = 1;
      const int numTaylorTerms = SpaceDim; // Must have at least SpaceDim equations to be able to get the stencil. 

      // Get Vofs, try to use quadrants but if the normal is aligned with the grid just get a "symmetric" stencil. 
      Vector<VolIndex> allVofs;
      if(VofUtils::isQuadrantWellDefined(normal)){ // Try to use quadrants. 
	allVofs = VofUtils::getAllVofsInQuadrant(a_vof, a_ebisbox, normal, radius, addStartingVof);
      }
      else{
	const std::pair<int, Side::LoHiSide> cardinal = VofUtils::getCardinalDirection(normal);
	allVofs = VofUtils::getAllVofsSymmetric(a_vof, a_ebisbox, cardinal.first, cardinal.second, radius, addStartingVof);
      }

      // Now build the stencil. 
      if(allVofs.size() >= numTaylorTerms){
	const Vector<RealVect> displacements = LeastSquares::getDisplacements(CellPosition::Boundary,
									      CellPosition::Center,
									      a_vof,
									      allVofs,
									      a_ebisbox,
									      a_dx);
    
	bndrySten = LeastSquares::computeGradStenOrderOne(allVofs, displacements, a_p);
      }
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

  int LeastSquares::getTaylorExpansionSize(const int a_order){

    MultiIndex mi(a_order);
  
    return mi.getNumIndices();
  }

  VoFStencil LeastSquares::computeInterpolationStencil(const Vector<VolIndex>& a_allVofs,
						       const Vector<RealVect>& a_displacements,
						       const int               a_pow,
						       const int               a_order){

    const Vector<Real> weights = LeastSquares::makeDiagWeights(a_displacements, a_pow);

    return LeastSquares::computeInterpolationStencil(a_allVofs, a_displacements, weights, a_order);
  }
  
  VoFStencil LeastSquares::computeInterpolationStencil(const Vector<VolIndex>& a_allVofs,
						       const Vector<RealVect>& a_displacements,
						       const Vector<Real>&     a_weights,
						       const int               a_order){

    VoFStencil ret;
  
    const IntVect deriv = IntVect::Zero;
    const IntVectSet ivs(deriv);

    std::map<IntVect, VoFStencil> allStens = LeastSquares::computeInterpolationStencil(ivs, a_allVofs, a_displacements, a_weights, a_order);

    ret = allStens.at(deriv);

    return ret;
  }

  std::map<IntVect, VoFStencil> LeastSquares::computeInterpolationStencil(const IntVectSet&       a_derivs,
									  const Vector<VolIndex>& a_allVofs,
									  const Vector<RealVect>& a_displacements,
									  const int               a_p,
									  const int               a_order){

    const Vector<Real> weights = LeastSquares::makeDiagWeights(a_displacements, a_p);

    return LeastSquares::computeInterpolationStencil(a_derivs, a_allVofs, a_displacements, weights, a_order);
    
  }

  std::map<IntVect, VoFStencil> LeastSquares::computeInterpolationStencil(const IntVectSet&       a_derivs,
									  const Vector<VolIndex>& a_allVofs,
									  const Vector<RealVect>& a_displacements,
									  const Vector<Real>&     a_weights,
									  const int               a_order){
    std::map<IntVect, VoFStencil> ret;

    // Initialize return. 
    for (IVSIterator ivsIt(a_derivs); ivsIt.ok(); ++ivsIt){
      ret.emplace(ivsIt(), VoFStencil());
    }
  
    if(a_derivs.numPts() > 0){

      const int M = LeastSquares::getTaylorExpansionSize(a_order);
      const int K = a_displacements.size();

      if(K < M) MayDay::Abort("LeastSquares::computeInterpolationStencil -- not enough equations to achieve desired order!");

      // Build the A-matrix in column major order so we can use LaPackUtils::computePseudoInverse.
      // Use of multi-indices makes higher-order Taylor series a walk in the park. 
      int i = 0;
      Vector<Real> linA    (K*M, 0.0);
      Vector<Real> linAplus(M*K, 0.0);
      for (MultiIndex mi(a_order); mi.ok(); ++mi){
    
	for (int k = 0; k < K; k++){
	  linA[i] = a_weights[k]*mi.pow(a_displacements[k])/mi.factorial(); i++;
	}
      }
 
      // Compute the pseudo-inverse.
      const bool foundSVD = LaPackUtils::computePseudoInverse(linAplus.stdVector(), linA.stdVector(), K, M);

      if(foundSVD){
	const MultiIndex mi(a_order);

	// Recall that linAplus is M*K so the stride is always M, starting at some specified row.
	for (IVSIterator ivsIt(a_derivs); ivsIt.ok(); ++ivsIt){
	  const IntVect deriv = ivsIt();
	
	  const int row = mi.getLinearIndex(ivsIt());

	  VoFStencil& sten = ret.at(deriv);

	  sten.clear();
	  for (int k = 0; k < K; k++){
	    const int idx = row + k*M;
	    sten.add(a_allVofs[k], a_weights[k]*linAplus[idx]);
	  }
	}
      }
    }

    return ret;
  }

  VoFStencil LeastSquares::getInterpolationStencilUsingAllVofsInRadius(const CellPosition a_cellPos,
								       const CellPosition a_otherCellsPos,
								       const VolIndex&    a_vof,
								       const EBISBox&     a_ebisbox,
								       const Real         a_dx,
								       const int          a_p,
								       const int          a_radius,
								       const int          a_order){

    // Get all Vofs in a radius, and then compute the displacement vectors. 
    Vector<VolIndex> vofs = VofUtils::getAllConnectedVofsInRadius(a_vof, a_ebisbox, a_radius, IntVectSet());
    vofs.push_back(a_vof);
    const Vector<RealVect> displacements = LeastSquares::getDisplacements(a_cellPos, a_otherCellsPos, a_vof, vofs, a_ebisbox, a_dx);

    const int M = LeastSquares::getTaylorExpansionSize(a_order);
    const int K = displacements.size();

    VoFStencil ret;
    if(K >= M) { // Have enough equations to compute.
      ret = LeastSquares::computeInterpolationStencil(vofs, displacements, a_p, a_order);
    }

    return ret;
  }
}
