/*!
  @file LeastSquares.cpp
  @brief Implementation of LeastSquares.hpp
  @author Robert Marskar
  @date May 2016
*/

#include "LeastSquares.H"
#include "LaPackUtils.H"

#include "EBArith.H"

RealVect LeastSquares::position(const CellPosition a_position,
				const VolIndex&    a_vof,
				const EBISBox&     a_ebisbox,
				const Real&        a_dx){

  RealVect ret;
  switch(a_position){
  case CellPosition::Center:
    ret = RealVect(a_vof.gridIndex());
    break;
  case CellPosition::Centroid:
    ret = RealVect(a_vof.gridIndex()) + a_ebisbox.centroid(a_vof);;
    break;
  case CellPosition::Boundary:
    ret = RealVect(a_vof.gridIndex()) + a_ebisbox.bndryCentroid(a_vof);;
    break;
  }

  return ret;
}

RealVect LeastSquares::displacement(const CellPosition      a_from,
				    const CellPosition      a_to,
				    const VolIndex&         a_fromVoF,
				    const VolIndex&         a_toVoF,
				    const EBISBox&          a_ebisbox,
				    const Real&             a_dx){

  const RealVect a = LeastSquares::position(a_from, a_fromVoF, a_ebisbox, a_dx);
  const RealVect b = LeastSquares::position(a_to,   a_toVoF,   a_ebisbox, a_dx);

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

Vector<RealVect> LeastSquares::getCenterToCentroidDisplacements(const VolIndex&         a_curVoF,
								const Vector<VolIndex>& a_allVoFs,
								const EBISBox&          a_ebisbox,
								const Real&             a_dx){
  Vector<RealVect> displacements(a_allVoFs.size());
  const RealVect centroid = RealVect(a_curVoF.gridIndex()) + a_ebisbox.centroid(a_curVoF);
  
  for (int i = 0; i < a_allVoFs.size(); i++){
    const RealVect center = RealVect(a_allVoFs[i].gridIndex());
    displacements[i] = (center - centroid)*a_dx;
  }

  return displacements;
}

Vector<RealVect> LeastSquares::getCenterToEBDisplacements(const VolIndex&         a_vof,
							  const Vector<VolIndex>& a_monoVoFs,
							  const EBISBox&          a_ebisbox,
							  const Real&             a_dx){
  Vector<RealVect> displacements(a_monoVoFs.size());
  const RealVect ebCentroid = RealVect(a_vof.gridIndex()) + a_ebisbox.bndryCentroid(a_vof);
  for (int i = 0; i < a_monoVoFs.size(); i++){
    const RealVect center = RealVect(a_monoVoFs[i].gridIndex());
    displacements[i] = (center - ebCentroid)*a_dx;;
  }

  return displacements;
}

Vector<RealVect> LeastSquares::getCentroidToCenterDisplacements(const VolIndex&         a_vof,
								const Vector<VolIndex>& a_monoVoFs,
								const EBISBox&          a_ebisbox,
								const Real&             a_dx){

  // Build a vector of distances from the cell centroid
  Vector<RealVect> displacements(a_monoVoFs.size());
  const RealVect center = RealVect(a_vof.gridIndex());
  
  for (int i = 0; i < a_monoVoFs.size(); i++){
    const RealVect centroid = RealVect(a_monoVoFs[i].gridIndex()) + a_ebisbox.centroid(a_monoVoFs[i]);
    displacements[i] = (centroid - center)*a_dx;
  }

  return displacements;
}

bool LeastSquares::getGradSten(VoFStencil&             a_stencil,
			       const Vector<VolIndex>& a_monoVoFs,
			       const Vector<RealVect>& a_displ,
			       const int               a_Q,
			       const int               a_weight){
  bool foundStencil;

  // Build a vector of weights for the least squares method. The weights are given by the inverse distance
  Vector<Real> wi = makeDiagWeights(a_displ, a_weight);

  // Size of the linear system
  Vector<MultiIndex> indices = LeastSquares::getMultiIndicesLexiOrder(a_Q);
  int k;
  int M = LeastSquares::getTaylorExpansionSize(a_Q);
  int N = LeastSquares::getTaylorExpansionSize(a_Q);
  int K = a_displ.size();

  // Build the A-matrix so we can use LaPackUtils::computePseudoInverse
  k = 0;
  Vector<Real> linA(M*N);
  for (MultiIndex q(IntVect::Zero); q <= a_Q; q.next(a_Q)){
    for (MultiIndex p(IntVect::Zero); p <= a_Q; p.next(a_Q)){
      
      linA[k] = 0.;
      for (int i = 0; i < K; i++){
	linA[k] += wi[i]*wi[i]*MultiIndex::pow(a_displ[i], q)*MultiIndex::pow(a_displ[i], p)/(q.factorial()*p.factorial());
      }

      k++;
    }
  }

  // Build the B-matrix
  k = 0;
  double B[N*K];
  for (int i = 0; i < K; i++){
    for (MultiIndex p(IntVect::Zero); p <= a_Q; p.next(a_Q)){
      B[k] = wi[i]*wi[i]*MultiIndex::pow(a_displ[i], p)/p.factorial();

      k++;
    }
  }

  // Let LaPack compute the Moore-Penrose pseudoinverse
  Vector<Real> linAPlus(M*N, 0.0);
  const bool foundSVD = LaPackUtils::computePseudoInverse(linAPlus.stdVector(), linA.stdVector(), M, N);

  if(foundSVD){

    // Must have the pseudoinverse in something useable by Fortran
    double Aplus[M*N];
    for (int i = 0; i < M*N; i++){
      Aplus[i] = linAPlus[i];
    }

    // Now compute Aplus*B which will be an M*K matrix
    double C[M*K];
    char TRANSA = 'N';
    char TRANSB = 'N';
    double ALPHA = 1.0;
    double BETA  = 0.0;
    dgemm_(&TRANSA,
	   &TRANSB,
	   &M,
	   &K,
	   &M,
	   &ALPHA,
	   Aplus,
	   &M,
	   B,
	   &M,
	   &BETA,
	   C,
	   &M);

    if(a_Q == 1){
      // The stencil is now found in the 2nd, 3rd, and 4th rows of C, which is an (M x K) matrix
      a_stencil.clear();
      for (int dir = 0; dir < SpaceDim; dir++){
	for (int j = 0; j < K; j++){
	  const VolIndex& jvof = a_monoVoFs[j]; 
	  const int k = M*j + (dir + 1);
	  a_stencil.add(jvof, C[k], dir);
	}
	foundStencil = true;
      }
    }
    else {
      MayDay::Abort("LeastSquares::getLSqGradStencil only works for a_Q = 1. This is a work in progress...");
    }
  }
  else{
    foundStencil = false;
  }

  return foundStencil;
}

bool LeastSquares::getBndryGradStencil(VoFStencil&     a_stencil,
				       const VolIndex& a_vof,
				       const EBISBox&  a_ebisbox,
				       const Real&     a_dx,
				       const int       a_order,
				       const int       a_radius,
				       const int       a_weight){
  int radius = a_radius;
  Vector<VolIndex> monoVoFs;
  const int numTaylorTerms = LeastSquares::getTaylorExpansionSize(a_order);

  // Can or must we drop order?
  const bool dropOrder = LeastSquares::getVoFsRadius(monoVoFs, a_vof, a_ebisbox, numTaylorTerms, a_order, a_radius);

  bool foundStencil = false;
  
  if(!dropOrder){
    // Get the displacement vectors
    Vector<RealVect> displacements = LeastSquares::getCenterToEBDisplacements(a_vof, monoVoFs, a_ebisbox, a_dx);

    // Build least squarse system and solve it
    foundStencil = LeastSquares::getGradSten(a_stencil, monoVoFs, displacements, a_order, a_weight);

    // Drop order if we could not find stencil
    if(!foundStencil){
      pout() << "LeastSquares::getBndryGradStencil - could not find stencil" << endl;
      a_stencil.clear();
    }
  }

  return foundStencil;
}

bool LeastSquares::getCenterToCentroidInterpStencil(VoFStencil&     a_stencil,
						    const VolIndex& a_vof,
						    const EBISBox&  a_ebisbox,
						    const Real&     a_dx,
						    const int       a_order,
						    const int       a_radius){

  bool foundStencil;
  int radius = a_radius;
  Vector<VolIndex> monoVoFs;
  const int numTaylorTerms = LeastSquares::getTaylorExpansionSize(a_order);

  // Can or must we drop order?
  bool dropOrder = false;
  const RealVect vofCentroid = a_ebisbox.centroid(a_vof);
  if(vofCentroid.vectorLength() < LeastSquares::s_tol){ 
    dropOrder = true;
  }
  else{
    dropOrder = LeastSquares::getVoFsRadius(monoVoFs, a_vof, a_ebisbox, numTaylorTerms, a_order, a_radius);
  }

  if(!dropOrder){
    // Get the displacement vectors
    Vector<RealVect> displacements = LeastSquares::getCenterToCentroidDisplacements(a_vof, monoVoFs, a_ebisbox, a_dx);

    // Build least squarse system and solve it
    foundStencil = LeastSquares::getLSqInterpStencil(a_stencil, monoVoFs, displacements, a_order);

    // Drop order if we could not find stencil
    if(!foundStencil){
      pout() << "LeastSquares::getCenterToCentroidInterpStencil - could not find stencil, dropping order" << endl;
      a_stencil.clear();
      a_stencil.add(a_vof, 1.0);
    }
  }
  else{

    // If the centroid is at the center, get the stencil directly. 
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
    foundStencil = true;
  }

  return foundStencil;
}

bool LeastSquares::getCenterToEBInterpStencil(VoFStencil&     a_stencil,
					      const VolIndex& a_vof,
					      const EBISBox&  a_ebisbox,
					      const Real&     a_dx,
					      const int       a_order,
					      const int       a_radius){
  bool foundStencil;
  int radius = a_radius;
  Vector<VolIndex> monoVoFs;
  const int numTaylorTerms = LeastSquares::getTaylorExpansionSize(a_order);

  // Can or must we drop order?
  bool dropOrder = false;
  const RealVect vofCentroid = a_ebisbox.bndryCentroid(a_vof);
  if(vofCentroid.vectorLength() < LeastSquares::s_tol){ 
    dropOrder = true;
  }
  else{
    dropOrder = LeastSquares::getVoFsRadius(monoVoFs, a_vof, a_ebisbox, numTaylorTerms, a_order, a_radius);
  }

  if(!dropOrder){
    // Get the displacement vectors
    Vector<RealVect> displacements = LeastSquares::getCenterToEBDisplacements(a_vof, monoVoFs, a_ebisbox, a_dx);

    // Build least squarse system and solve it
    foundStencil = LeastSquares::getLSqInterpStencil(a_stencil, monoVoFs, displacements, a_order);

    // Drop order if we could not find stencil
    if(!foundStencil){
      pout() << "LeastSquares::getCenterToEB - could not find stencil, dropping order" << endl;
      a_stencil.clear();
      a_stencil.add(a_vof, 1.0);
    }
  }
  else{

    // If the centroid is at the center, get the stencil directly. 
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
    foundStencil = true;
  }

  return foundStencil;
}



bool LeastSquares::getCentroidToCenterInterpStencil(VoFStencil&     a_stencil,
						    const VolIndex& a_vof,
						    const EBISBox&  a_ebisbox,
						    const Real&     a_dx,
						    const int       a_order,
						    const int       a_radius){
  bool foundStencil;
  int radius = a_radius;
  Vector<VolIndex> monoVoFs;
  const int numTaylorTerms = LeastSquares::getTaylorExpansionSize(a_order);

  // Can or must we drop order?
  bool dropOrder = false;
  const RealVect vofCentroid = a_ebisbox.centroid(a_vof);
  if(vofCentroid.vectorLength() < LeastSquares::s_tol){ 
    dropOrder = true;
  }
  else{
    dropOrder = LeastSquares::getVoFsRadius(monoVoFs, a_vof, a_ebisbox, numTaylorTerms, a_order, a_radius);
  }

  if(!dropOrder){

    // Get the displacement vectors
    Vector<RealVect> displacements = LeastSquares::getCentroidToCenterDisplacements(a_vof, monoVoFs, a_ebisbox, a_dx);

    // Build least squarse system and solve it
    foundStencil = LeastSquares::getLSqInterpStencil(a_stencil, monoVoFs, displacements, a_order);

    // Drop order if we could not find stencil
    if(!foundStencil){
      pout() << "LeastSquares::getCentroidToCenterInterpStencil - could not find stencil, dropping order" << endl;
      a_stencil.clear();
      a_stencil.add(a_vof, 1.0);
    }
  }
  else{

    // If the centroid is at the center, get the stencil directly. 
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
    foundStencil = true;
  }

  return foundStencil;
}

Vector<MultiIndex> LeastSquares::getMultiIndicesLexiOrder(const int a_Q){
  Vector<MultiIndex> ret;
  for (MultiIndex cur = IntVect::Zero; cur <= a_Q; cur.next(a_Q)){
    ret.push_back(cur);
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

bool LeastSquares::getLSqInterpStencil(VoFStencil&             a_stencil,
				       const Vector<VolIndex>& a_monoVoFs,
				       const Vector<RealVect>& a_displ,
				       const int               a_Q){
  bool foundStencil;

  // Build a vector of weights for the least squares method. The weights are given by the inverse distance
  Vector<Real> wi;
  for (int i = 0; i < a_displ.size(); i++){
    wi.push_back(1.0/sqrt(a_displ[i].vectorLength()));
  }

  // Size of the linear system
  Vector<MultiIndex> indices = LeastSquares::getMultiIndicesLexiOrder(a_Q);
  int k;
  int M = LeastSquares::getTaylorExpansionSize(a_Q);
  int N = LeastSquares::getTaylorExpansionSize(a_Q);
  int K = a_displ.size();
  CH_assert(K >= M);

  // Build the A-matrix so we can use LaPackUtils::computePseudoInverse
  k = 0;
  Vector<Real> linA(M*N);
  for (MultiIndex q(IntVect::Zero); q <= a_Q; q.next(a_Q)){
    for (MultiIndex p(IntVect::Zero); p <= a_Q; p.next(a_Q)){
      
      //
      linA[k] = 0.;
      for (int i = 0; i < K; i++){
	linA[k] += wi[i]*wi[i]*MultiIndex::pow(a_displ[i], q)*MultiIndex::pow(a_displ[i], p)/(q.factorial()*p.factorial());
      }

      k++;
    }
  }

  // Build the B-matrix
  k = 0;
  double B[N*K];
  for (int i = 0; i < K; i++){
    for (MultiIndex p(IntVect::Zero); p <= a_Q; p.next(a_Q)){
      B[k] = wi[i]*wi[i]*MultiIndex::pow(a_displ[i], p)/p.factorial();

      k++;
    }
  }

  // Let LaPack compute the Moore-Penrose pseudoinverse
  Vector<Real> linAPlus(M*N, 0.0);
  const bool foundSVD = LaPackUtils::computePseudoInverse(linAPlus.stdVector(), linA.stdVector(), M, N);

  if(foundSVD){

    // Must have the pseudoinverse in something useable by Fortran
    double Aplus[M*N];
    for (int i = 0; i < M*N; i++){
      Aplus[i] = linAPlus[i];
    }

    // Now compute Aplus*B which will be an M*K matrix
    double C[M*K];
    char TRANSA = 'N';
    char TRANSB = 'N';
    double ALPHA = 1.0;
    double BETA  = 0.0;
    dgemm_(&TRANSA,
	   &TRANSB,
	   &M,
	   &K,
	   &M,
	   &ALPHA,
	   Aplus,
	   &M,
	   B,
	   &M,
	   &BETA,
	   C,
	   &M);
    

    // The stencil is now found in the first row of C.
    a_stencil.clear();
    for (int j = 0; j < K; j++){
      const VolIndex& jvof = a_monoVoFs[j]; 
      const int k = j*M;
      a_stencil.add(jvof, C[k]);
    }
    foundStencil = true;
  }
  else{
    foundStencil = false;
  }

  return foundStencil;
}



bool LeastSquares::getVoFsRadius(Vector<VolIndex>& a_allVoFsRad,
				 const VolIndex&   a_vof,
				 const EBISBox&    a_ebisbox,
				 const int         a_minVoFs,
				 const int         a_order, 
				 const int         a_radius){
  // Get VoFs
  int radius     = a_radius;
  bool dropOrder = false;

  if(a_radius == -1){
    radius = a_order - 1;
    bool enoughVoFs = false;
    while(!enoughVoFs){

      EBArith::getAllVoFsInMonotonePath(a_allVoFsRad, a_vof, a_ebisbox, radius);
      a_allVoFsRad.push_back(a_vof);

      if(a_allVoFsRad.size() < a_minVoFs){
	radius += 1;
      }
      else{
	enoughVoFs = true;
	dropOrder  = false;
      }
    }
  }
  else{
    EBArith::getAllVoFsInMonotonePath(a_allVoFsRad, a_vof, a_ebisbox, radius);
    a_allVoFsRad.push_back(a_vof);

    dropOrder = a_allVoFsRad.size() < a_minVoFs ? true : false;
  }

  return dropOrder;
}
