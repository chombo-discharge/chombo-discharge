/*!
  @file LeastSquares.cpp
  @brief Implementation of LeastSquares.hpp
  @author Robert Marskar
  @date May 2016
*/

#include "LeastSquares.H"
#include "LaPackUtils.H"

#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "EBArith.H"

Vector<Real> LeastSquares::makeDiagWeights(const Vector<RealVect>& a_displacements, const Real a_pow){
  Vector<Real> ret(a_displacements.size(), 1.0);

  for (int i = 0; i < a_displacements.size(); i++){
    const Real len = a_displacements[i].vectorLength();
    ret[i] = 1./std::pow(len, a_pow);
  }

  return ret;
}

void LeastSquares::getCenterToCentroidDisplacements(Vector<RealVect>&       a_dist,
						    const VolIndex&         a_vof,
						    const Vector<VolIndex>& a_monoVoFs,
						    const EBISBox&          a_ebisbox,
						    const Real&             a_dx){
  // Build a vector of distances from the cell centroid
  a_dist.resize(a_monoVoFs.size());
  const RealVect centroid = RealVect(a_vof.gridIndex()) + a_ebisbox.centroid(a_vof);
  for (int i = 0; i < a_monoVoFs.size(); i++){
    const RealVect center = RealVect(a_monoVoFs[i].gridIndex());
    a_dist[i] = (center - centroid)*a_dx;
  }
}

void LeastSquares::getCenterToEbCentroidDisplacements(Vector<RealVect>&       a_dist,
						      const VolIndex&         a_vof,
						      const Vector<VolIndex>& a_monoVoFs,
						      const EBISBox&          a_ebisbox,
						      const Real&             a_dx){
  // Build a vector of distances from the cell centroid
  a_dist.resize(a_monoVoFs.size());
  const RealVect ebCentroid = RealVect(a_vof.gridIndex()) + a_ebisbox.bndryCentroid(a_vof);
  for (int i = 0; i < a_monoVoFs.size(); i++){
    const RealVect center = RealVect(a_monoVoFs[i].gridIndex());
    a_dist[i] = (center - ebCentroid)*a_dx;;
  }
}

void LeastSquares::getCentroidToCenterDisplacements(Vector<RealVect>&       a_dist,
						    const VolIndex&         a_vof,
						    const Vector<VolIndex>& a_monoVoFs,
						    const EBISBox&          a_ebisbox,
						    const Real&             a_dx){

  // Build a vector of distances from the cell centroid
  a_dist.resize(a_monoVoFs.size());
  const RealVect center = RealVect(a_vof.gridIndex());
  
  for (int i = 0; i < a_monoVoFs.size(); i++){
    const VolIndex& ivof = a_monoVoFs[i];
    
    const RealVect centroid = RealVect(ivof.gridIndex()) + a_ebisbox.centroid(ivof);
    a_dist[i] = (centroid - center)*a_dx;
  }
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
    Vector<RealVect> displacements;
    LeastSquares::getCenterToCentroidDisplacements(displacements, a_vof, monoVoFs, a_ebisbox, a_dx);

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

bool LeastSquares::getCenterToEbCentroidInterpStencil(VoFStencil&     a_stencil,
						      const VolIndex& a_vof,
						      const EBISBox&  a_ebisbox,
						      const Real&     a_dx,
						      const int       a_order,
						      const int       a_radius){
  CH_TIME("LeastSquares::getCenterToEbCentroidInterpStencil");

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
    Vector<RealVect> displacements;
    LeastSquares::getCenterToEbCentroidDisplacements(displacements, a_vof, monoVoFs, a_ebisbox, a_dx);

    // Build least squarse system and solve it
    foundStencil = LeastSquares::getLSqInterpStencil(a_stencil, monoVoFs, displacements, a_order);

    // Drop order if we could not find stencil
    if(!foundStencil){
      pout() << "LeastSquares::getCenterToEbCentroid - could not find stencil, dropping order" << endl;
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

bool LeastSquares::getEbCentroidGradientStencil(VoFStencil&     a_stencil,
						const VolIndex& a_vof,
						const EBISBox&  a_ebisbox,
						const Real&     a_dx,
						const int       a_order,
						const int       a_radius){
  CH_TIME("LeastSquares::getEbCentroidGradientStencil");

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

  bool foundStencil = false;
  if(!dropOrder){
    // Get the displacement vectors
    Vector<RealVect> displacements;
    LeastSquares::getCenterToEbCentroidDisplacements(displacements, a_vof, monoVoFs, a_ebisbox, a_dx);

    // Build least squarse system and solve it
    foundStencil = LeastSquares::getLSqGradientStencil(a_stencil, monoVoFs, displacements, a_order);

    // Drop order if we could not find stencil
    if(!foundStencil){
      pout() << "LeastSquares::getEbCentroidGradientStencil - could not find stencil" << endl;
      a_stencil.clear();
    }
  }

  return foundStencil;
}

bool LeastSquares::getCentroidToCenterInterpStencil(VoFStencil&     a_stencil,
						    const VolIndex& a_vof,
						    const EBISBox&  a_ebisbox,
						    const Real&     a_dx,
						    const int       a_order,
						    const int       a_radius){
  CH_TIME("LeastSquares::getCentroidToCenterInterpStencil");

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
    Vector<RealVect> displacements;
    LeastSquares::getCentroidToCenterDisplacements(displacements, a_vof, monoVoFs, a_ebisbox, a_dx);

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
  CH_TIME("LeastSquares::getMultiIndicesLexiOrder");


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
  CH_TIME("LeastSquares::getLSqInterpStencil");

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
  const bool foundSVD = LaPackUtils::computePseudoInverse(linAPlus, linA, M, N);

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

bool LeastSquares::getLSqGradientStencil(VoFStencil&             a_stencil,
					 const Vector<VolIndex>& a_monoVoFs,
					 const Vector<RealVect>& a_displ,
					 const int               a_Q){
  CH_TIME("LeastSquares::getLSqGradientStencil");

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
  const bool foundSVD = LaPackUtils::computePseudoInverse(linAPlus, linA, M, N);

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

bool LeastSquares::getVoFsRadius(Vector<VolIndex>& a_allVoFsRad,
				 const VolIndex&   a_vof,
				 const EBISBox&    a_ebisbox,
				 const int         a_minVoFs,
				 const int         a_order, 
				 const int         a_radius){
  CH_TIME("LeastSquares::getVoFsRadius");

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
