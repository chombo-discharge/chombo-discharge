/*!
  @file   LaPackUtils.H
  @brief  Agglomeration of various useful LaPack calls
  @author Robert Marskar
  @date   Feb. 2017
*/

#ifndef _LAPACKUTILS_
#define _LAPACKUTILS_

#include "Vector.H"
#include "RealVect.H"

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
  @brief Aggregation of various useful LaPack routines
 */  
namespace LaPackUtils {

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
  bool computeSVD(Vector<Real>&       a_linU,
		  Vector<Real>&       a_linSigma,
		  Vector<Real>&       a_linVtran,
		  const Vector<Real>& a_linA,
		  const int&          a_M,
		  const int&          a_N);

  /*!
    @brief Compute the pseudoinverse of matrix through singular value decomposition
    @param[out] a_linAplus Pseudoinverse in column major Fortran order. 
    @param[in] a_linA Input matrix in column major Fortran order
    @param[in] a_M Number of rows in A
    @param[in] a_N Number of colums in A
  */
  bool computePseudoInverse(Vector<Real>&       a_linAplus,
			    const Vector<Real>& a_linA,
			    const int&          a_M,
			    const int&          a_N);
  /*!
    @brief Linearize a matrix to column major Fortran form by assuming row or major colum format of the input matrix. 
    @param[out] a_linA Linearized matrix
    @param[out] a_M Number of rows
    @param[out] a_N Number of columns
    @param[in] a_A Input matrix
    @param[in] a_format Format of input matrix, either 'C' or 'R'. If 'C', then a_A[0] is the first column of the matrix. Likewise, if 'R' then a_A[0] is the first row of the matrix. 
  */
  void linearizeMatrix(Vector<Real>& a_linA, int& a_M, int& a_N, const Vector<Vector<Real> >& a_A, const char& a_format);

  /*!
    @brief Linearize a matrix to column major Fortran form by column major format of the input matrix. 
    @param[out] a_linA Linearized matrix
    @param[out] a_M Number of rows
    @param[out] a_N Number of columns
    @param[in] a_A Input matrix
  */
  void linearizeColumnMajorMatrix(Vector<Real>& a_linA, int& a_M, int& a_N, const Vector<Vector<Real> >& a_A);

  
  /*!
    @brief Linearize a matrix to column major Fortran form by row major format of the input matrix. 
    @param[out] a_linA Linearized matrix
    @param[out] a_M Number of rows
    @param[out] a_N Number of columns
    @param[in] a_A Input matrix
  */
  void linearizeRowMajorMatrix(Vector<Real>& a_linA, int& a_M, int& a_N, const Vector<Vector<Real> >& a_A);



  /*!
    @brief Delinearize a linearized matrix from column major Fortran to row or major colum matrix format.
    @param[out] a_A Output matrix
    @param[in] a_M Number of rows
    @param[in] a_N Number of columns
    @param[in] a_linA Linearized matrix
    @param[in] a_format Format of input matrix, either 'C' or 'R'. If 'C', then a_A[0] is the first column of the matrix. 
    Likewise, if 'R' then a_A[0] is the first row of the matrix. 
  */
  void deLinearizeMatrix(Vector<Vector<Real> >& a_A,
			 const int&             a_M,
			 const int&             a_N,
			 const Vector<Real>&    a_linA,
			 const char&            a_format);

  /*!
    @brief Delinearize a linearized matrix from column major Fortran to column major form
    @param[out] a_A Output matrix
    @param[in] a_M Number of rows
    @param[in] a_N Number of columns
    @param[in] a_linA Linearized matrix
  */
  void deLinearizeColumnMajorMatrix(Vector<Vector<Real> >& a_A,
				    const int&             a_M,
				    const int&             a_N,
				    const Vector<Real>&    a_linA);
				   
  /*!
    @brief Delinearize a linearized matrix from column major Fortran to row major form
    @param[out] a_A Output matrix
    @param[in] a_M Number of rows
    @param[in] a_N Number of columns
    @param[in] a_linA Linearized matrix
  */
  void deLinearizeRowMajorMatrix(Vector<Vector<Real> >& a_A,
				 const int&             a_M,
				 const int&             a_N,
				 const Vector<Real>&    a_linA);
}
#endif
