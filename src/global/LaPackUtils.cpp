/*!
  @file   LaPackUtils.cpp
  @brief  Implementation of LaPackUtils.H
  @author Robert Marskar
  @date   Feb. 2017
*/

#include "LaPackUtils.H"
#include <limits>
#include <algorithm>

#include "CH_Timer.H"

bool LaPackUtils::computeSVD(std::vector<double>&       a_linU,
			     std::vector<double>&       a_linSigma,
			     std::vector<double>&       a_linVT,
			     const std::vector<double>& a_linA,
			     const int&                 a_M,
			     const int&                 a_N){
  int M = a_M;
  int N = a_N;

  // linA must be stored as a double array in order to interface into LaPack
  double A[M*N];
  for (int i = 0; i < a_linA.size(); i++){
    A[i] = a_linA[i];
  }

  const int mx = std::max(M,N);
  const int mn = std::min(M,N);

  // This is the stuff that is needed for LaPack
  char JOBZ  = 'A';
  int LDA    = M;
  int LDU    = M;
  int LDVT   = N;
  int INFO   = 0;
  int LWORK  = 4*mn*mn + 6*mn + mx;
  int IWORK[8*std::min(M,N)];
  double S[std::min(M, N)];
  double U[M*M];
  double VT[N*N];
  double WORK[std::max(1, LWORK)];

  // Do the SVD decomposition
  dgesdd_(&JOBZ,
	  &M,
	  &N,
	  A,
	  &LDA,
	  S,
	  U,
	  &LDU,
	  VT,
	  &LDVT,
	  WORK,
	  &LWORK,
	  IWORK,
	  &INFO);

  const bool foundSVD = (INFO == 0);

  if(foundSVD){
    //
    a_linU.resize(M*M, 0.0);
    a_linSigma.resize(M*N, 0.0);
    a_linVT.resize(N*N, 0.0);
    
    // Build the linearized matrices
    for (int i = 0; i < M*M; i++){
      a_linU[i]  = U[i];
    }
    for (int i = 0; i < N*N; i++){
      a_linVT[i]  = VT[i];
    }
    for (int i = 0; i < std::min(M, N); i ++){
      a_linSigma[i*(M+1)] = S[i];
    }
  }

  return foundSVD;
}

bool LaPackUtils::computePseudoInverse(std::vector<double>&       a_linAplus,
				       const std::vector<double>& a_linA,
				       const int&                 a_M,
				       const int&                 a_N){
  // Compute the singular value decomposition
  std::vector<double> linU, linVT, linSigma;
  const bool foundSVD = computeSVD(linU, linSigma, linVT, a_linA, a_M, a_N);

  if(foundSVD){

    int M = a_M;
    int N = a_N;

    // Pseudoinversion of Sigma can be nasty for singular values close to zero. Define a tolerance for this. 
    double maxS = -1.E99;
    for (int i = 0; i < linSigma.size(); i++){
      maxS = linSigma[i] > maxS ? linSigma[i] : maxS;
    }
    const double eps = std::numeric_limits<double>::epsilon();
    const double tol = eps*std::max(M,N)*maxS;
    
    // Need to storage the matrices in a form useable by LaPack, and then use dgemm to multiply them. 
    double U[M*M];
    double SigmaReciprocal[M*N];
    double VT[N*N];

    for (int i = 0; i < linU.size(); i++){
      U[i] = linU[i];
    }
    for (int i = 0; i < linVT.size(); i++){
      VT[i] = linVT[i];
    }
    for (int i = 0; i < linSigma.size(); i++){
      if(std::abs(linSigma[i]) > tol){
	SigmaReciprocal[i] = 1./linSigma[i];
      }
      else{
	SigmaReciprocal[i] = 0.0;
      }
    }

    // Compute C = V*Transpose(SigmaReciprocal) by using dgemm
    double C[N*M];
    char TRANSA = 'T';
    char TRANSB = 'T';
    double ALPHA = 1.0;
    double BETA  = 0.0;
    dgemm_(&TRANSA,
	   &TRANSB,
	   &N,
	   &M,
	   &N,
	   &ALPHA,
	   VT,
	   &N,
	   SigmaReciprocal,
	   &M,
	   &BETA,
	   C,
	   &N);

    // Now compute Aplus = C*Transpose(U)
    double Aplus[N*M];
    TRANSA = 'N';
    TRANSB = 'T';
    dgemm_(&TRANSA,
	   &TRANSB,
	   &N,
	   &M,
	   &M,
	   &ALPHA,
	   C,
	   &N,
	   U,
	   &M,
	   &BETA,
	   Aplus,
	   &N);


    // Put Aplus into the output vector
    a_linAplus.resize(N*M);
    for (int i = 0; i < N*M; i++){
      a_linAplus[i] = Aplus[i];
    }
  }

  return foundSVD;
}

void LaPackUtils::linearizeColumnMajorMatrix(std::vector<double>&                     a_linA,
					     int&                                     a_M,
					     int&                                     a_N,
					     const std::vector<std::vector<double> >& a_A){
  // Number of rows and columns
  a_M = a_A[0].size();  // Number of rows of the actual matrix
  a_N = a_A.size();     // Number of columns of the actual matrix

  // Linearize the input matrix. 
  a_linA.resize(a_M*a_N);
  for (int j = 0; j < a_N; j++){   // j = Matrix column
    for (int i = 0; i < a_M; i++){ // i = Matrix row 
      const int k = j*a_M + i;
      a_linA[k] = a_A[j][i];
    }
  }
}

void LaPackUtils::linearizeRowMajorMatrix(std::vector<double>&                     a_linA,
					  int&                                     a_M,
					  int&                                     a_N,
					  const std::vector<std::vector<double> >& a_A){
  // Number of rows and columns
  a_M = a_A.size();    
  a_N = a_A[0].size(); 

  // Linearize the input matrix
  a_linA.resize(a_M*a_N);
  for (int j = 0; j < a_N; j++){   // j = Matrix column
    for (int i = 0; i < a_M; i++){ // i = Matrix row 
      const int k = j*a_M + i;
      a_linA[k] = a_A[i][j];
    }
  }
}

void LaPackUtils::linearizeMatrix(std::vector<double>&                     a_linA,
				  int&                                     a_M,
				  int&                                     a_N,
				  const std::vector<std::vector<double> >& a_A,
				  const char&                              a_format){
  if(a_format == 'C'){
    linearizeColumnMajorMatrix(a_linA, a_M, a_N, a_A);
  }
  else if(a_format == 'R'){
    linearizeRowMajorMatrix(a_linA, a_M, a_N, a_A);
  }
  else{
    std::cerr << "LaPackUtils::deLinearizeMatrix - unknown specification of matrix storage\n";
  }
}

void LaPackUtils::deLinearizeColumnMajorMatrix(std::vector<std::vector<double> >& a_A,
					       const int&                         a_M,
					       const int&                         a_N,
					       const std::vector<double>&         a_linA){
  // Resize the matrix
  a_A.resize(a_N);
  for (int i = 0; i < a_N; i++){
    a_A[i].resize(a_M);
  }
  
  // Delinearize it
  for (int j = 0; j < a_N; j++){   // j = Matrix column
    for (int i = 0; i < a_M; i++){ // i = Matrix row 
      const int k = j*a_M + i;
      a_A[j][i] = a_linA[k];
    }
  } 
}

void LaPackUtils::deLinearizeRowMajorMatrix(std::vector<std::vector<double> >& a_A,
					    const int&                         a_M,
					    const int&                         a_N,
					    const std::vector<double>&         a_linA){
  // Resize the matrix
  a_A.resize(a_M);
  for (int i = 0; i < a_M; i++){
    a_A[i].resize(a_N);
  }
  
  // Delinearize it
  for (int j = 0; j < a_N; j++){   // j = Matrix column
    for (int i = 0; i < a_M; i++){ // i = Matrix row 
      const int k = j*a_M + i;
      a_A[i][j] = a_linA[k];
    }
  }  
}

void LaPackUtils::deLinearizeMatrix(std::vector<std::vector<double> >& a_A,
				    const int&                         a_M,
				    const int&                         a_N,
				    const std::vector<double>&         a_linA,
				    const char&                        a_format){
  if(a_format == 'C'){
    deLinearizeColumnMajorMatrix(a_A, a_M, a_N, a_linA);
  }
  else if(a_format == 'R'){
    deLinearizeRowMajorMatrix(a_A, a_M, a_N, a_linA);
  }
  else{
    std::cerr << "LaPackUtils::deLinearizeMatrix - unknown specification of matrix storage\n";
  }
}
