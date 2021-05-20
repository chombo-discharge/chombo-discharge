/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_LaPackUtils.cpp
  @brief  Implementation of CD_LaPackUtils.H
  @author Robert Marskar
*/

// Std includes
#include <limits>
#include <algorithm>
#include <iostream>

// Our includes
#include <CD_LaPackUtils.H>
#include <CD_NamespaceHeader.H>

int LaPackUtils::linearIndex(const int irow, const int jcol, const int M, const int N){
  return irow + jcol*M;
}

bool LaPackUtils::computeSVD(std::vector<double>&       a_linU,
			     std::vector<double>&       a_linSigma,
			     std::vector<double>&       a_linVT,
			     const std::vector<double>& a_linA,
			     const int&                 a_M,
			     const int&                 a_N){


  // linA must be stored as a double array in order to interface into LaPack
  double A[a_M*a_N];
  for (int i = 0; i < a_linA.size(); i++){
    A[i] = a_linA[i];
  }

  const int mx = std::max(a_M,a_N);
  const int mn = std::min(a_M,a_N);

  // This is the stuff that is needed for LaPack
  char JOBZ  = 'A';
  int M      = a_M;
  int N      = a_N;
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
    for (int i = 0; i < std::min(M, N); i++){
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
  std::vector<double> linU, linSigma, linVT;
  const bool foundSVD = computeSVD(linU, linSigma, linVT, a_linA, a_M, a_N);

  if(foundSVD){

    // Pseudoinversion of Sigma can be nasty for singular values close to zero. Define a tolerance for this. 
    double maxS = std::numeric_limits<double>::min();
    for (const auto& s : linSigma){
      maxS = std::max(maxS, s);
    }
    const double eps = std::numeric_limits<double>::epsilon();
    const double tol = eps*std::max(a_M,a_N)*maxS;
    
    // Need to storage the matrices in a form useable by LaPack, and then use dgemm to multiply them. 
    double U[a_M*a_M];
    double SigmaReciprocal[a_M*a_N];
    double VT[a_N*a_N];

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

    // Compute C = V*Transpose(SigmaReciprocal) by using dgemm. The dimensions are:
    //
    // VT = N*N               => op(A) is N*N
    // SigmaReciprocal = M*N  => op(B) is N*M
    //
    // We want to compute V*Transpose(SigmaReciprocal) onto C, which is N*M big.
    //
    double C[a_N*a_M];
    {
      char TRANSA  = 'T';
      char TRANSB  = 'T';
      double ALPHA = 1.0;
      double BETA  = 0.0;
      int M        = a_N; // VT is N*N;
      int N        = a_M; // B = SigmaReciprocal which is M*N so op(B) = B**T is N*M
      int K        = a_N;
      int LDA      = K;
      int LDB      = N;
      int LDC      = M;
      dgemm_(&TRANSA,
	     &TRANSB,
	     &M,
	     &N,
	     &K,
	     &ALPHA,
	     VT,
	     &LDA,
	     SigmaReciprocal,
	     &LDB,
	     &BETA,
	     C,
	     &LDC);
    }

    // Now compute Aplus = C*Transpose(U). C is a N*M matrix. Here, U is an M*M matrix and C is an N*M matrix.
    // So we have
    //
    // op(A) => M*M
    // op(B) => N*M
    double Aplus[a_N*a_M];
    {
      char TRANSA  = 'N';
      char TRANSB  = 'T';
      double ALPHA = 1.0;
      double BETA  = 0.0;
      int M        = a_N;
      int N        = a_M;
      int K        = a_M;
      int LDA      = M;
      int LDB      = K;
      int LDC      = M;
      dgemm_(&TRANSA,
	     &TRANSB,
	     &M,
	     &N,
	     &K,
	     &ALPHA,
	     C,
	     &LDA,
	     U,
	     &LDB,
	     &BETA,
	     Aplus,
	     &LDC);
    }


    // Put Aplus into the output vector
    a_linAplus.resize(a_N*a_M);
    for (int i = 0; i < a_N*a_M; i++){
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

#include <CD_NamespaceFooter.H>
