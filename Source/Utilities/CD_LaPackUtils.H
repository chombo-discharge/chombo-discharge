/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_LaPackUtils.H
  @brief  Interface to some LaPack routines. 
  @note   The routines in this file are double precision only. 
  @author Robert Marskar
*/

#ifndef CD_LaPackUtils_H
#define CD_LaPackUtils_H

// Std includes
#include <vector>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief This is an interface to LAPACK for solving least squares problems using singular-value decomposition
*/
extern "C" void dgelss_(int*    a_M, 
			int*    a_N, 
			int*    a_nRHS, 
			double* a_A, 
			int*    a_LDA, 
			double* a_B, 
			int*    a_LDB, 
			double* a_S, 
			double* a_rcond,
			int*    a_RANK,
			double* a_WORK,
			int*    a_LWORK,
			int*    a_INFO);

/*!
  @brief Interface to LaPack for computing the singular value decomposition of a matrix. 
*/
extern "C" void dgesdd_(char*   a_JOBZ,
			int*    a_M,
			int*    a_N,
			double* a_A,
			int*    a_LDA, 
			double* a_S, 
			double* a_U,
			int*    a_LDU,
			double* a_VT,
			int*    a_LDVT,
			double* a_WORK,
			int*    a_LWORK,
			int*    a_IWORK,
			int*    a_INFO);

/*!
  @brief Interface to LaPack for matrix multiplication
*/
extern "C" void dgemm_(char*   a_TRANSA,
		       char*   a_TRANSB,
		       int*    a_M,
		       int*    a_N,
		       int*    a_K,
		       double* a_ALPHA,
		       double* a_A,
		       int*    a_LDA,
		       double* a_B,
		       int*    a_LDB,
		       double* a_BETA,
		       double* a_C,
		       int*    a_LDC);

/*!
  @brief Interface to LaPack for computing solution to Ax=b
*/
extern "C" void dgesv_(int*    N,
		       int*    NRHS,
		       double* A,
		       int*    LDA,
		       int*    IPIV,
		       double* B,
		       int*    LDB,
		       int*    INFO);



namespace LaPackUtils {

  /*!
    @brief Get the index in the linearized matrix of size MxN
    @param[in] irow The row number
    @param[in] jcol The column number
    @param[in] M Number of rows
    @param[in] N Number columns
  */
  int linearIndex(const int irow, const int jcol, const int M, const int N);

  /*!
    @brief Compute the singular value decomposition of a matrix. 
    @param[out] a_linU \f$U\f$-matrix. Linearized in column major Fortran order. 
    @param[out] a_linSigma \f$\Sigma\f$-matrix. Linearized in column major Fortran order. 
    @param[out] a_linVtran \f$V^\intercal\f$-matrix. Linearized in column major Fortran order. 
    @param[in]  a_linA \f$A\f$-matrix. Stored in column-major Fortran form
    @param[in]  a_M Number of rows in A
    @param[in]  a_N Number of columns in A
    @details This routine will compute the singular value decomposition of the input matrix \f$A\f$ in the form
    \f[
    A = U\Sigma V^\intercal.
    \f]
    The input and output matrices are assumed to be stored column-major so that a_A[0] is the first column of \f$A\f$. 
  */
  bool computeSVD(std::vector<double>&       a_linU,
		  std::vector<double>&       a_linSigma,
		  std::vector<double>&       a_linVtran,
		  const std::vector<double>& a_linA,
		  const int&                 a_M,
		  const int&                 a_N);

  /*!
    @brief Compute the pseudoinverse of matrix through singular value decomposition
    @param[out] a_linAplus Pseudoinverse in column major Fortran order. 
    @param[in]  a_linA Input matrix in column major Fortran order
    @param[in]  a_M Number of rows in A
    @param[in]  a_N Number of colums in A
  */
  bool computePseudoInverse(std::vector<double>&       a_linAplus,
			    const std::vector<double>& a_linA,
			    const int&                 a_M,
			    const int&                 a_N);
  /*!
    @brief Linearize a matrix to column major Fortran form by assuming row or major colum format of the input matrix. 
    @param[out] a_linA Linearized matrix
    @param[out] a_M Number of rows
    @param[out] a_N Number of columns
    @param[in]  a_A Input matrix
    @param[in]  a_format Format of input matrix, either 'C' or 'R'. If 'C', then a_A[0] is the first column of the matrix. 
    Likewise, if 'R' then a_A[0] is the first row of the matrix. 
  */
  void linearizeMatrix(std::vector<double>& a_linA, int& a_M, int& a_N, const std::vector<std::vector<double> >& a_A, const char& a_format);

  /*!
    @brief Linearize a matrix to column major Fortran form by column major format of the input matrix. 
    @param[out] a_linA Linearized matrix
    @param[out] a_M Number of rows
    @param[out] a_N Number of columns
    @param[in]  a_A Input matrix
  */
  void linearizeColumnMajorMatrix(std::vector<double>& a_linA, int& a_M, int& a_N, const std::vector<std::vector<double> >& a_A);
  
  /*!
    @brief Linearize a matrix to column major Fortran form by row major format of the input matrix. 
    @param[out] a_linA Linearized matrix
    @param[out] a_M Number of rows
    @param[out] a_N Number of columns
    @param[in]  a_A Input matrix
  */
  void linearizeRowMajorMatrix(std::vector<double>& a_linA, int& a_M, int& a_N, const std::vector<std::vector<double> >& a_A);

  /*!
    @brief Delinearize a linearized matrix from column major Fortran to row or major colum matrix format.
    @param[out] a_A Output matrix
    @param[in]  a_M Number of rows
    @param[in]  a_N Number of columns
    @param[in]  a_linA Linearized matrix
    @param[in]  a_format Format of input matrix, either 'C' or 'R'. If 'C', then a_A[0] is the first column of the matrix. 
    Likewise, if 'R' then a_A[0] is the first row of the matrix. 
  */
  void deLinearizeMatrix(std::vector<std::vector<double> >& a_A,
			 const int&                         a_M,
			 const int&                         a_N,
			 const std::vector<double>&         a_linA,
			 const char&                        a_format);

  /*!
    @brief Delinearize a linearized matrix from column major Fortran to column major form
    @param[out] a_A Output matrix
    @param[in]  a_M Number of rows
    @param[in]  a_N Number of columns
    @param[in]  a_linA Linearized matrix
  */
  void deLinearizeColumnMajorMatrix(std::vector<std::vector<double> >& a_A,
				    const int&                         a_M,
				    const int&                         a_N,
				    const std::vector<double>&         a_linA);
				   
  /*!
    @brief Delinearize a linearized matrix from column major Fortran to row major form
    @param[out] a_A Output matrix
    @param[in]  a_M Number of rows
    @param[in]  a_N Number of columns
    @param[in]  a_linA Linearized matrix
  */
  void deLinearizeRowMajorMatrix(std::vector<std::vector<double> >& a_A,
				 const int&                         a_M,
				 const int&                         a_N,
				 const std::vector<double>&         a_linA);
}

#include <CD_NamespaceFooter.H>

#endif
