/*!
  @file LeastSquares.cpp
  @brief Implementation of LeastSquares.hpp
  @author Robert Marskar
  @date May 2016
*/

#include "LeastSquares.H"
#include "LaPackUtils.H"

#include "EBArith.H"

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
    ret = RealVect(a_vof.gridIndex()) + a_ebisbox.centroid(a_vof);
    break;
  case CellPosition::Boundary:
    ret = RealVect(a_vof.gridIndex()) + a_ebisbox.bndryCentroid(a_vof);
    break;
  }

  ret *= a_dx;

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

VoFStencil LeastSquares::getGradStenOrderOne(const Vector<VolIndex>& a_allVoFs,
					     const Vector<RealVect>& a_displacements,
					     const int&              a_p){

  VoFStencil sten;

  // Make sure we have enough equations
  const int K = a_displacements.size();
  if(K < 3) MayDay::Abort("LeastSquares::getGradStenOrderOne -- not enough equations, must have K >= 3");

  // Build the matrix for which we need the pseudoinverse. Build it in column-major order.
  const Vector<Real> wi = LeastSquares::makeDiagWeights(a_displacements, a_p);
  const int D = SpaceDim;
  
  Vector<Real> linA(K*D);
  for (int d = 0; d < D; d++){
    for (int k = 0; k < K; k++){
      const Real& delta = a_displacements[k][d];
      
      const int idx = k + d*K;

      linA[idx] = wi[k]*delta;
    }
  }

  // Compute pseudoinverse
  Vector<Real> linAplus(K*D, 0.0);
  const bool foundSVD = LaPackUtils::computePseudoInverse(linAplus.stdVector(), linA.stdVector(), K, D);

  if(foundSVD){

    double Aplus[D*K];
    for (int i = 0; i < linAplus.size(); i++){
      Aplus[i] = linAplus[i];
    }
    
    // Now build the diagonal W matrix so we can call dgemm (maybe unnecessary...?)
    double W[K*K];
    for (int i = 0; i < K; i++){
      for (int j = 0; j < K; j++){
	const int k = i + j*K;
	if(i == j){
	  W[k] = wi[i];
	}
	else{
	  W[k] = 0.0;
	}
      }
    }

    // Now compute Aplus*B which will be an M*K matrix
    double C[D*K];
    char TRANSA = 'N';
    char TRANSB = 'N';
    double ALPHA = 1.0;
    double BETA  = 0.0;
    int M = D;
    int N = K;
    dgemm_(&TRANSA,
	   &TRANSB,
	   &M,
	   &N,
	   &N,
	   &ALPHA,
	   Aplus,
	   &M,
	   W,
	   &N,
	   &BETA,
	   C,
	   &M);

    
    sten.clear();
    for (int k = 0; k < K; k++){
	for (int dir = 0; dir < D; dir++){
	const Real weight   = C[k*D + dir];
	const VolIndex& vof = a_allVoFs[k]; // Columns are VoFs
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
  const int order  = 1;
  const int radius = 1;

  VoFStencil sten;
  const int numTaylorTerms = LeastSquares::getTaylorExpansionSize(order);

  // Find VoFs in radius
  Vector<VolIndex> monoVoFs;
  const bool dropOrder = LeastSquares::getVoFsRadius(monoVoFs, a_vof, a_ebisbox, numTaylorTerms, 1, 1);

  if(!dropOrder){
    Vector<RealVect> displacements = LeastSquares::getDisplacements(CellPosition::Boundary,
								    CellPosition::Center,
								    a_vof,
								    monoVoFs,
								    a_ebisbox,
								    a_dx);
    
    // Build least squarse system and solve it
    sten = LeastSquares::getGradStenOrderOne(monoVoFs, displacements, a_p);
  }

  return sten;
}

bool LeastSquares::getGradSten(VoFStencil&             a_stencil,
			       const Vector<VolIndex>& a_monoVoFs,
			       const Vector<RealVect>& a_displ,
			       const int               a_Q,
			       const int               a_weight){
  bool foundStencil;

  // Build a vector of weights for the least squares method. The weights are given by the inverse distance
  Vector<Real> wi = LeastSquares::makeDiagWeights(a_displ, a_weight);

  // Size of the linear system
  Vector<MultiIndex> indices = LeastSquares::getMultiIndicesLexiOrder(a_Q);
  int M = LeastSquares::getTaylorExpansionSize(a_Q);
  int N = LeastSquares::getTaylorExpansionSize(a_Q);
  int K = a_displ.size();

  // Build the A-matrix so we can use LaPackUtils::computePseudoInverse
  int k = 0;
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

  const int numTaylorTerms = LeastSquares::getTaylorExpansionSize(a_order);

  // Find VoFs in radius
  int radius = a_radius;
  Vector<VolIndex> monoVoFs;
  const bool dropOrder = LeastSquares::getVoFsRadius(monoVoFs, a_vof, a_ebisbox, numTaylorTerms, a_order, a_radius);

  bool foundStencil = false;
  
  if(!dropOrder){
    // Get the displacement vectors
    //    Vector<RealVect> displacements = LeastSquares::getCenterToEBDisplacements(a_vof, monoVoFs, a_ebisbox, a_dx);
    Vector<RealVect> displacements = LeastSquares::getDisplacements(CellPosition::Boundary,
								    CellPosition::Center,
								    a_vof,
								    monoVoFs,
								    a_ebisbox,
								    a_dx);

    // Build least squarse system and solve it
    foundStencil = LeastSquares::getGradSten(a_stencil, monoVoFs, displacements, a_order, a_weight);

    // Drop order if we could not find stencil
    if(!foundStencil){
      a_stencil.clear();
    }
#if 0 // debug
    else{
      RealVect sumWeights = RealVect::Zero;
      for (int i = 0; i < a_stencil.size(); i++){
	const Real w  = a_stencil.weight(i);
	const int var = a_stencil.variable(i);

	sumWeights[var] += w;
      }
      std::cout << "sumweights = " << sumWeights << std::endl;
    }
#endif
  }

  return foundStencil;
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
