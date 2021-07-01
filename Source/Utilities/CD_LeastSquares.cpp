/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LeastSquares.cpp
  @brief  Implementation of CD_LeastSquares.H
  @author Robert Marskar
*/

// Std includes
#include <vector>
#include <set>

// Our includes
#include <CD_LaPackUtils.H>
#include <CD_LeastSquares.H>
#include <CD_MultiIndex.H>
#include <CD_NamespaceHeader.H>

VoFStencil LeastSquares::getInterpolationStencil(const CellLocation a_cellPos,
						 const CellLocation a_otherCellsPos,
						 const Connectivity a_connectivity,						 
						 const VolIndex&    a_startVof,
						 const EBISBox&     a_ebisbox,
						 const Real         a_dx,
						 const int          a_p,
						 const int          a_radius,
						 const int          a_order,
						 const bool         a_addStartingVof){

  if(a_p > 0 && a_addStartingVof){
    MayDay::Error("LeastSquares::getInterpolationStencilUsingAllConnectedVofsInRadius - can't use a_p > 0 && a_addStartingVof (it's an ill-conditioned system of equations)");
  }

  // Get all Vofs in a radius, and then compute the displacement vectors. 
  Vector<VolIndex> vofs = VofUtils::getVofsInRadius(a_startVof, a_ebisbox, a_radius, a_connectivity, a_addStartingVof);
  const Vector<RealVect> displacements = LeastSquares::getDisplacements(a_cellPos, a_otherCellsPos, a_startVof, vofs, a_ebisbox, a_dx);

  const int M = LeastSquares::getTaylorExpansionSize(a_order);
  const int K = displacements.size();

  VoFStencil ret;
  if(K >= M) { // Have enough equations to compute.
    ret = LeastSquares::computeInterpolationStencil(vofs, displacements, a_p, a_order);
  }

  return ret;
}

VoFStencil LeastSquares::getGradSten(const VolIndex&    a_vof,
				     const CellLocation a_gradLocation,
				     const CellLocation a_cellLocation,
				     const EBISBox&     a_ebisbox,
				     const Real         a_dx,
				     const int          a_radius,
				     const int          a_p,
				     const int          a_order,
				     const IntVectSet   a_knownTerms){

  VoFStencil gradSten;


  const Vector<VolIndex> allVofs = VofUtils::getVofsInRadius(a_vof, a_ebisbox, a_radius, VofUtils::Connectivity::MonotonePath, true);

  const int numUnknowns = LeastSquares::getTaylorExpansionSize(a_order) - a_knownTerms.numPts(); // Can be smaller than order if you've eliminated some equations. 

  // Build the stencil if we can. 
  if(allVofs.size() >= numUnknowns){

    const Vector<RealVect> displacements = LeastSquares::getDisplacements(a_gradLocation,
									  a_cellLocation,
									  a_vof,
									  allVofs,
									  a_ebisbox,
									  a_dx);

      
    gradSten = LeastSquares::computeGradSten(allVofs, displacements, a_p, a_order, a_knownTerms); 
  }


  return gradSten;
}

VoFStencil LeastSquares::getGradSten(const FaceIndex&   a_face,
				     const FaceLocation a_gradLocation,
				     const CellLocation a_cellLocation,
				     const EBISBox&     a_ebisbox,
				     const Real         a_dx,
				     const int          a_radius,
				     const int          a_p,
				     const int          a_order,
				     const IntVectSet   a_knownTerms){
  VoFStencil gradSten;

  if(!a_face.isBoundary()){
    const int faceDir = a_face.direction();

    const VolIndex vofLo = a_face.getVoF(Side::Lo);
    const VolIndex vofHi = a_face.getVoF(Side::Hi);

    // Get Vofs in monotone path for both lo/hi. The doing that will fetch duplicates (starting vof, for example) which we discard
    // using std::set below. After that, just get the distances and solve the least squares system. 
    Vector<VolIndex> allVofs;
    Vector<VolIndex> loVofs = VofUtils::getVofsInRadius(vofLo, a_ebisbox, a_radius, VofUtils::Connectivity::MonotonePath, true);
    Vector<VolIndex> hiVofs = VofUtils::getVofsInRadius(vofHi, a_ebisbox, a_radius, VofUtils::Connectivity::MonotonePath, true);

    std::set<VolIndex> setVofs;
    for (const auto& v : loVofs.stdVector()) setVofs.emplace(v);
    for (const auto& v : hiVofs.stdVector()) setVofs.emplace(v);
    for (const auto& v : setVofs           ) allVofs.push_back(v);

    const int numUnknowns = LeastSquares::getTaylorExpansionSize(a_order) - a_knownTerms.numPts(); // Can be smaller than order if you've eliminated some equations. 

    if(allVofs.size() >= numUnknowns){
      const Vector<RealVect> displacements = LeastSquares::getDisplacements(a_gradLocation, a_cellLocation, a_face, allVofs, a_ebisbox, a_dx);

      gradSten = LeastSquares::computeGradSten(allVofs, displacements, a_p, a_order, a_knownTerms);
    }
  }

  return gradSten;
}

VoFStencil LeastSquares::getBndryGradSten(const VolIndex&    a_vof,
					  const Neighborhood a_neighborhood,
					  const CellLocation a_cellPositions,
					  const EBISBox&     a_ebisbox,
					  const Real         a_dx,
					  const int          a_radius,
					  const int          a_p,
					  const int          a_order,
					  const bool         a_addStartingVof){
  VoFStencil bndrySten;

  const RealVect normal     = a_ebisbox.normal(a_vof);

  if(normal != RealVect::Zero){ // Can only do this if we actual have a normal vector. 
    const int numUnknowns = LeastSquares::getTaylorExpansionSize(a_order) - 1; 

    Vector<VolIndex> allVofs;

    switch(a_neighborhood){
    case Neighborhood::Quadrant:
      allVofs = VofUtils::getVofsInQuadrant(a_vof, a_ebisbox, normal, a_radius, VofUtils::Connectivity::MonotonePath, a_addStartingVof);
      break;
    case Neighborhood::Radius:
      allVofs = VofUtils::getVofsInRadius(a_vof, a_ebisbox, a_radius, VofUtils::Connectivity::MonotonePath, a_addStartingVof);
      break;
    default:
      break;
    }

    // Now build the stencil. 
    if(allVofs.size() >= numUnknowns){

      const Vector<RealVect> displacements = LeastSquares::getDisplacements(Location::Cell::Boundary,
									    a_cellPositions,
									    a_vof,
									    allVofs,
									    a_ebisbox,
									    a_dx);

      IntVectSet knownTerms;
      knownTerms |= IntVect::Zero;
      bndrySten = LeastSquares::computeGradSten(allVofs, displacements, a_p, a_order, knownTerms); // This routine eliminates a_vof from the system of equations!
    }
  }

  return bndrySten;
}

RealVect LeastSquares::position(const CellLocation a_position,
				const VolIndex&    a_vof,
				const EBISBox&     a_ebisbox,
				const Real&        a_dx){

  RealVect ret;
  
  switch(a_position){
  case Location::Cell::Center:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit;
    break;
  case Location::Cell::Centroid:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit + a_ebisbox.centroid(a_vof);
    break;
  case Location::Cell::Boundary:
    ret = RealVect(a_vof.gridIndex()) + 0.5*RealVect::Unit + a_ebisbox.bndryCentroid(a_vof);
    break;
  default:
    MayDay::Error("LeastSquares::position -- bad location input.");
    break;
  }

  ret *= a_dx;

  return ret;
}

RealVect LeastSquares::position(const FaceLocation a_position,
				const FaceIndex&   a_face,
				const EBISBox&     a_ebisbox,
				const Real&        a_dx){
  RealVect ret;
  const IntVect iv = a_face.gridIndex(Side::Hi);
  
  switch(a_position){
  case Location::Face::Center:
    ret = (RealVect(iv) + 0.5*RealVect::Unit) - 0.5*BASISREALV(a_face.direction());
    break;
  case Location::Face::Centroid:
    ret = (RealVect(iv) + 0.5*RealVect::Unit) - 0.5*BASISREALV(a_face.direction()) + a_ebisbox.centroid(a_face);
    break;
  default:
    MayDay::Error("LeastSquares::positon - bad location input.");
    break;
  }

  ret *= a_dx;

  return ret;
}

RealVect LeastSquares::displacement(const CellLocation a_from,
				    const CellLocation a_to,
				    const VolIndex&    a_fromVof,
				    const VolIndex&    a_toVof,
				    const EBISBox&     a_ebisbox,
				    const Real&        a_dx){

  const RealVect a = LeastSquares::position(a_from, a_fromVof, a_ebisbox, a_dx);
  const RealVect b = LeastSquares::position(a_to,   a_toVof,   a_ebisbox, a_dx);

  return (b-a);
}

RealVect LeastSquares::displacement(const FaceLocation a_fromLoc,
				    const CellLocation a_toLoc,
				    const FaceIndex&   a_fromFace,
				    const VolIndex&    a_toVof,
				    const EBISBox&     a_ebisbox,
				    const Real&        a_dx){
  const RealVect a = LeastSquares::position(a_fromLoc, a_fromFace, a_ebisbox, a_dx);
  const RealVect b = LeastSquares::position(a_toLoc,   a_toVof,    a_ebisbox, a_dx);

  return (b-a);
}

RealVect LeastSquares::displacement(const CellLocation a_from,
				    const CellLocation a_to,
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

Vector<RealVect> LeastSquares::getDisplacements(const CellLocation      a_from,
						const CellLocation      a_to,
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

Vector<RealVect> LeastSquares::getDisplacements(const FaceLocation      a_fromLoc,
						const CellLocation      a_toLoc,
						const FaceIndex&        a_fromFace,
						const Vector<VolIndex>& a_toVofs,
						const EBISBox&          a_ebisbox,
						const Real&             a_dx){
  Vector<RealVect> ret;

  for (int i = 0; i < a_toVofs.size(); i++){
    const RealVect d = LeastSquares::displacement(a_fromLoc, a_toLoc, a_fromFace, a_toVofs[i], a_ebisbox, a_dx);
    
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
      const Real w      = std::pow(len, -a_pow);

      ret[i] = w;
    }
  }

  return ret;
}

VoFStencil LeastSquares::computeGradSten(const Vector<VolIndex>& a_allVofs,
					 const Vector<RealVect>& a_displacements,
					 const int               a_p,
					 const int               a_order,
					 const IntVectSet        a_knownTerms){
  Vector<Real> weights = LeastSquares::makeDiagWeights(a_displacements, a_p);

  return LeastSquares::computeGradSten(a_allVofs, a_displacements, weights, a_order, a_knownTerms);
}

VoFStencil LeastSquares::computeGradSten(const Vector<VolIndex>& a_allVofs,
					 const Vector<RealVect>& a_displacements,
					 const Vector<Real>&     a_weights,
					 const int               a_order,
					 const IntVectSet        a_knownTerms){


  // TLDR: This routine 
  IntVectSet derivs;
  for (int dir = 0; dir < SpaceDim; dir++){
    derivs |= BASISV(dir);
  }

  std::map<IntVect, VoFStencil> taylorTerms = LeastSquares::computeSingleLevelStencils(derivs, a_knownTerms, a_allVofs, a_displacements, a_weights, a_order);

  VoFStencil sten;
  for (const auto& m : taylorTerms){
    for (int dir = 0; dir < SpaceDim; dir++){
      if (m.first == BASISV(dir)) {
	  
	const VoFStencil& iSten = m.second;
	  
	for (int i = 0; i < iSten.size(); i++){
	  const VolIndex& vof = iSten.vof(i);
	  const Real& weight  = iSten.weight(i);

	  sten.add(vof, weight, dir);
	}
      }
    }
  }

  return sten;

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
  const IntVectSet derivs(deriv);

  const IntVectSet knownTerms = IntVectSet();

  std::map<IntVect, VoFStencil> allStens = LeastSquares::computeSingleLevelStencils(derivs, knownTerms, a_allVofs, a_displacements, a_weights, a_order);

  ret = allStens.at(deriv);

  return ret;
}

std::map<IntVect, VoFStencil> LeastSquares::computeSingleLevelStencils(const IntVectSet&       a_derivs,
								       const IntVectSet&       a_knownTerms,
								       const Vector<VolIndex>& a_allVofs,
								       const Vector<RealVect>& a_displacements,
								       const int               a_p,
								       const int               a_order){

  const Vector<Real> weights = LeastSquares::makeDiagWeights(a_displacements, a_p);

  return LeastSquares::computeSingleLevelStencils(a_derivs, a_knownTerms, a_allVofs, a_displacements, weights, a_order);
    
}

std::map<IntVect, VoFStencil> LeastSquares::computeSingleLevelStencils(const IntVectSet&       a_derivs,
								       const IntVectSet&       a_knownTerms,
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

    const int M = LeastSquares::getTaylorExpansionSize(a_order) - a_knownTerms.numPts(); // This is because some unknowns have been eliminated. 
    const int K = a_displacements.size();

    const IntVectSet isect = a_derivs & a_knownTerms;

    if(K < M)            MayDay::Abort("LeastSquares::computeSingleLevelStencils -- not enough equations to achieve desired order!");
    if(!isect.isEmpty()) MayDay::Abort("LeastSquares::computeSingleLevelStencils - you have specified the same terms as both unknown and known");


    // Build the A-matrix in column major order so we can use LaPackUtils::computePseudoInverse.
    // Use of multi-indices makes higher-order Taylor series a walk in the park. 
    int i = 0;
    Vector<Real> linA    (K*M, 0.0);
    Vector<Real> linAplus(M*K, 0.0);
    for (MultiIndex mi(a_order); mi.ok(); ++mi){ // Loop over column

      if(!a_knownTerms.contains(mi.getCurrentIndex())){ // Add column if it is an unknown in the lsq system. 
	for (int k = 0; k < K; k++){ 
	  linA[i] = a_weights[k]*mi.pow(a_displacements[k])/mi.factorial(); i++;
	}
      }
    }
 
    // Compute the pseudo-inverse.
    const bool foundSVD = LaPackUtils::computePseudoInverse(linAplus.stdVector(), linA.stdVector(), K, M);

    if(foundSVD){
      const MultiIndex mi(a_order);

      // When we have eliminated rows in the linear system we can't use MultiIndex to map directly. 
      // Rather, we need to first map the rows of A to indices and use that to fetch
      // the row in Aplus corresponding to a specific unknown in the Taylor series. 
      std::map<IntVect, int> rowMap;
      int row = 0;
      for (MultiIndex mi(a_order); mi.ok(); ++mi){
	if(!a_knownTerms.contains(mi.getCurrentIndex())){ // This is the order in which A was built. 
	  rowMap.emplace(mi.getCurrentIndex(), row);
	  row++;
	}
      }

      // Recall that linAplus is M*K so the stride is always M, starting at some specified row. 
      for (IVSIterator ivsIt(a_derivs); ivsIt.ok(); ++ivsIt){
	const IntVect deriv = ivsIt();

	int irow;
	if(rowMap.find(ivsIt()) != rowMap.end()){
	  row = rowMap.at(ivsIt());
	}
	else{
	  MayDay::Abort("LeastSquares::computeSingleLevelStencils -- map is out of range but this shouldn't happen!");
	}

	VoFStencil& sten = ret.at(deriv);

	sten.clear();
	for (int k = 0; k < K; k++){
	  const int idx = row + k*M;
	  sten.add(a_allVofs[k], a_weights[k]*linAplus[idx]);
	}
      }
    }
    else{
      MayDay::Warning("LeastSquares::computeSingleLevelStencils - could not perform singular value decomposition");
    }
  }

  return ret;
}

std::map<IntVect, std::pair<VoFStencil, VoFStencil> > LeastSquares::computeDualLevelStencils(const IntVectSet&       a_derivs,
											     const IntVectSet&       a_knownTerms,
											     const Vector<VolIndex>& a_fineVofs,
											     const Vector<VolIndex>& a_coarVofs,
											     const Vector<RealVect>& a_fineDisplacements,
											     const Vector<RealVect>& a_coarDisplacements,
											     const int               a_p,           
											     const int               a_order){
  const Vector<Real> fineWeights = LeastSquares::makeDiagWeights(a_fineDisplacements, a_p);
  const Vector<Real> coarWeights = LeastSquares::makeDiagWeights(a_coarDisplacements, a_p);

  return LeastSquares::computeDualLevelStencils(a_derivs,
						a_knownTerms,
						a_fineVofs,
						a_coarVofs,
						a_fineDisplacements,
						a_coarDisplacements,
						fineWeights,
						coarWeights,
						a_order);
}

std::map<IntVect, std::pair<VoFStencil, VoFStencil> > LeastSquares::computeDualLevelStencils(const IntVectSet&       a_derivs,
											     const IntVectSet&       a_knownTerms,
											     const Vector<VolIndex>& a_fineVofs,
											     const Vector<VolIndex>& a_coarVofs,
											     const Vector<RealVect>& a_fineDisplacements,
											     const Vector<RealVect>& a_coarDisplacements,
											     const Vector<Real>&     a_fineWeights,
											     const Vector<Real>&     a_coarWeights,
											     const int               a_order){
  // Initialize return stuff
  std::map<IntVect, std::pair<VoFStencil, VoFStencil> > ret;

  // Initialize return. 
  for (IVSIterator ivsIt(a_derivs); ivsIt.ok(); ++ivsIt){
    ret.emplace(ivsIt(), std::make_pair(VoFStencil(), VoFStencil()));
  }

  if(a_derivs.numPts() > 0){

    const int M     = LeastSquares::getTaylorExpansionSize(a_order) - a_knownTerms.numPts(); // This is because some unknowns have been eliminated. 
    const int Kfine = a_fineDisplacements.size();
    const int Kcoar = a_coarDisplacements.size();
    const int K     = Kfine + Kcoar;

    const IntVectSet isect = a_derivs & a_knownTerms;

    if(K < M)            MayDay::Abort("LeastSquares::computeSingleLevelStencils -- not enough equations to achieve desired order!");
    if(!isect.isEmpty()) MayDay::Abort("LeastSquares::computeSingleLevelStencils - you have specified the same terms as both unknown and known");


    // Build the A-matrix in column major order so we can use LaPackUtils::computePseudoInverse.
    // Use of multi-indices makes higher-order Taylor series a walk in the park. 
    int i = 0;
    Vector<Real> linA    (K*M, 0.0);
    Vector<Real> linAplus(M*K, 0.0);
    for (MultiIndex mi(a_order); mi.ok(); ++mi){ // Loop over column

      if(!a_knownTerms.contains(mi.getCurrentIndex())){ // Add column if it is an unknown in the lsq system. 
	for (int k = 0; k < K; k++){
	  if(k < Kfine){
	    linA[i] = a_fineWeights[k]*mi.pow(a_fineDisplacements[k])/mi.factorial(); i++;
	  }
	  else{
	    linA[i] = a_coarWeights[k]*mi.pow(a_coarDisplacements[k])/mi.factorial(); i++;
	  }
	}
      }
    }
 
    // Compute the pseudo-inverse.
    const bool foundSVD = LaPackUtils::computePseudoInverse(linAplus.stdVector(), linA.stdVector(), K, M);

    if(foundSVD){
      const MultiIndex mi(a_order);

      // When we have eliminated rows in the linear system we can't use MultiIndex to map directly. 
      // Rather, we need to first map the rows of A to indices and use that to fetch
      // the row in Aplus corresponding to a specific unknown in the Taylor series. 
      std::map<IntVect, int> rowMap;
      int row = 0;
      for (MultiIndex mi(a_order); mi.ok(); ++mi){
	if(!a_knownTerms.contains(mi.getCurrentIndex())){ // This is the order in which A was built. 
	  rowMap.emplace(mi.getCurrentIndex(), row);
	  row++;
	}
      }

      // Recall that linAplus is M*K so the stride is always M, starting at some specified row. 
      for (IVSIterator ivsIt(a_derivs); ivsIt.ok(); ++ivsIt){
	const IntVect deriv = ivsIt();

	int irow;
	if(rowMap.find(ivsIt()) != rowMap.end()){
	  row = rowMap.at(ivsIt());
	}
	else{
	  MayDay::Abort("LeastSquares::computeSingleLevelStencils -- map is out of range but this shouldn't happen!");
	}

	std::pair<VoFStencil, VoFStencil>& sten = ret.at(deriv);

	sten.first. clear();
	sten.second.clear();
	
	for (int k = 0; k < K; k++){
	  const int idx = row + k*M;
	  if(k < Kfine){
	    sten.first.add(a_fineVofs[k], a_fineWeights[k]*linAplus[idx]);
	  }
	  else{
	    sten.second.add(a_coarVofs[k], a_coarWeights[k]*linAplus[idx]);
	  }
	}
      }
    }
    else{
      MayDay::Warning("LeastSquares::computeDualLevelStencils - could not perform singular value decomposition");
    }
  }

  return ret;
}

#include <CD_NamespaceFooter.H>
