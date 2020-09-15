#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LAPACK_H_
#define _LAPACK_H_
#include "FortranNameMacro.H"
#include "BaseNamespaceHeader.H"

/// Forward declaration of lapack routines
extern "C" 
{
#ifdef CH_USE_LAPACK
  /// Double General Solve
  void LAPACK(GESV,gesv)(int*, int*, Real*, int*, int*, Real*, int*, int*);

  /// Double Least Squares 
  void LAPACK(GELS,gels)(char* TRANS, int* M, int* N, int* NRHS, Real* A, 
              int* LDA, Real* B, int* LDB, 
              Real* WORK, int* LWORK, int* INFO);

  /// Double Least Squares transpose via SVD
  void LAPACK(GELSD,gelsd)(int* N, int* M, int* NRHS, Real* At, 
               int* LDA, Real* B, int* LDB, Real* S,
               Real* RCOND, int* RANK, 
               Real* WORK, int* LWORK, int* IWORK, int* INFO);

  /// Matrix multiply
  void LAPACK(GEMV,gemv)(char* TRANS, int* M, int* N, Real* ALPHA, 
              Real* A, int* LDA, Real* X, int* INCX,
              Real* BETA, Real* Y, int* INCY);

  /// Double Equality-Constrained Least Squares 
  void LAPACK(GGLSE,gglse)(int* M, int* N, int* P, Real* A, int* LDA, Real* B, int* LDB, 
               Real* C, Real* D, Real* X, 
               Real* WORK, int* LWORK, int* INFO);

  /// QR factorization with Pivoting
  void LAPACK(GEQP3,geqp3)(int* M, int* N, Real* A, int* LDA, int* JPVT,
               Real* TAU, Real* WORK, int* LWORK, int* INFO);

  /// QR factorization
  void LAPACK(GEQRF,geqrf)(int* M, int* N, Real* A, int* LDA, 
               Real* TAU, Real* WORK, int* LWORK, int* INFO);

  /// Multiply by Q from QR factorization
  void LAPACK(ORMQR,ormqr)(char* SIDE, char* TRANS, int* M, int* N, int* K, 
               Real* A, int* LDA, Real* TAU, Real* C, int* LDC,
               Real* WORK, int* LWORK, int* INFO);

  /// Generate Q from QR factorization
  void LAPACK(ORGQR,orgqr)(int* M, int* N, int* K, Real* A, int* LDA, 
               Real* TAU, Real* WORK, int* LWORK, int* INFO);

  /// Triangular solve
  void LAPACK(TRTRS,trtrs)(char* UPLO, char* TRANS, char* DIAG, int* N, int* NRHS, 
               Real* A, int* LDA, Real* B, int* LDB, int* INFO);

  /// LU decomoposition of a general matrix
  void LAPACK(GETRF,getrf)(int* M, int *N, Real* A, int* lda, int* IPIV, int* INFO);

  /// generate inverse of a matrix given its LU decomposition
  void LAPACK(GETRI,getri)(int* N, Real* A, int* lda, int* IPIV, Real* WORK, int* lwork, int* INFO);

  /// get condition number of the matrix
  void LAPACK(GECON,gecon)(char* TRANS, int* N,  Real* A, int* LDA, Real* ANORM, Real* rcond, 
              Real* WORK, int* LWORK, int* INFO);

  /// get condition number of triangular matrix
  void LAPACK(TRCON,trcon)(const char* TRANS, const char* UPLO, const char* DIAG,
                    const int* N,  const Real* A, const int* LDA,
                    Real* rcond, Real* WORK, int* IWORK, int* INFO);

#endif
                    
  /// multiply two matrices (BLAS routine)
  void LAPACK(GEMM,gemm)(const char* TRANSA, const char* TRANSB,
                   const int* M, const int* N, const int* K, const Real* ALPHA,
                   const Real* A, const int* LDA,
                   const Real* B, const int* LDB,
                   const Real* BETA, Real* C, const int* LDC);

  /// solve matrix equations with triangular matrix (BLAS routine)
  void LAPACK(TRSM,trsm)(const char* SIDE, const char* UPLO,
                   const char* TRANSA, const char* DIAG,
                   const int* M, const int* N, const Real* ALPHA,
                   const Real* A, const int* LDA,
                   Real* B, const int* LDB);
}

#include "BaseNamespaceFooter.H"
#endif
