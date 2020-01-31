#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LAPACKMatrix_H_
#define _LAPACKMatrix_H_

#include "REAL.H"
#include "LoHiSide.H"
#include <utility>
#include "BaseNamespaceHeader.H"

///
/**
   Getting sick of writing the basics here over and over.
   Silly class but it will cut down on the typing.
 */
class LAPACKMatrix
{
public:
  

  ///  main constructor.  Matrix of Chombo Real type.  values are unintitialized
  LAPACKMatrix(int a_nrow, int a_ncol):m_nrow(a_nrow),m_ncol(a_ncol),m_data(new Real[a_nrow*a_ncol]){}

  ///  alias constructor.  Use a_data as m_data, does not delete it in destructor
  LAPACKMatrix( int a_nrow, int a_ncol, Real* a_data)
    :m_nrow(a_nrow),m_ncol(a_ncol),m_data(a_data),m_alias(true) {;}

  /// null constructor   builds matrix with nullptr and size=[0,0]
  LAPACKMatrix():m_nrow(0),m_ncol(0),m_data(nullptr){}

  ///  deep copy constructor
  LAPACKMatrix(const LAPACKMatrix& a_input);

  /// move copy constructor
  LAPACKMatrix(LAPACKMatrix&& a_input)
    :m_nrow(a_input.m_nrow),m_ncol(a_input.m_ncol),m_data(a_input.m_data),m_alias(a_input.m_alias)
  {a_input.m_data=nullptr;}




  ///return sqrt(sum of squares of all values)
  Real normLTwo() const;

  ///
  /**
    if(i==j)  M=1 else M=0   
   */
  void setToIdentity();

  ///
  ~LAPACKMatrix();

  ///
  void  setVal(const Real& a_val);

  Real* dataPtr();

  ///
  const Real* dataPtr() const;

  /// 
  /**
     make nrows = a_nrows etc.   Useful when you want to discard the last few rows of the matrix.
  */
  void truncate(int a_nrow,int a_ncol);


  ///
  /**
     Get the maximum absolute value of the matrix (for iterative solves)
   */
  Real maxNorm() const;

  ///
  const Real& operator() (int irow, int icol) const;
  
  ///
  Real & operator() (int irow, int icol);

  //return row, col
  std::pair<int, int>  dims() const
  {
    std::pair<int, int>  retval;
    retval.first = m_nrow;
    retval.second= m_ncol;
    return retval;
  }

  /// deep assign. 
  LAPACKMatrix& operator=(const LAPACKMatrix& a_matrix);

  /// move assignment
  inline LAPACKMatrix& operator=(LAPACKMatrix&& a_matrix)
  {
    if(&a_matrix != this)
      {
        m_nrow=a_matrix.m_nrow;
        m_ncol=a_matrix.m_ncol;
        std::swap(m_data,a_matrix.m_data);
        std::swap(m_alias,a_matrix.m_alias);
      }
    return *this;
  }
  
    
  ///
  LAPACKMatrix& operator+=(const LAPACKMatrix& a_matrix);

  ///
  LAPACKMatrix& operator-=(const LAPACKMatrix& a_matrix);

  ///
  LAPACKMatrix& operator*=(const Real& a_scalingFactor);


  ///
  int offset(int irow, int icol) const;

  ///
  void poutAll() const;

  ///
  void poutDiag() const;

  ///
  void poutDiagMatlab() const;

  ///
  void poutMatlab() const;
 
  ///inverts this matix 
  /**
     fails if matrix is not square
     if return value != 0, probably a singular matrix
   */
  int invert();

  /// inverts this matrix using SVD
  /**
   *  Get Ainverse using  least squares with svd
   *  if return value != 0, lapack was unhappy somehow 
   */
  int invertUsingSVD(int a_maxiter, Real a_tol);

  ///
  //pseudoinverse = (AT A)^(-1) AT
  int  pseudoInvertUsingSVD(int a_maxiter, Real a_tol);

  // pseudoinverse: X such that R*X = Q^T, where A = Q*R; hence A*X = Q*R*X = Q*Q^T = I
  int pseudoInvertUsingQR();

  /// inverts this matrix using least squares
  /**
   *  Get Ainverse using  least squares with svd
   *  if return value != 0, lapack was unhappy somehow 
   */
  int invertUsingLeastSquares();


  //for tiny cells, set moments to their limits (1 0 0 0...)
  void setSmallCellRow(const int& irow)
  {
    LAPACKMatrix& Mvol = *this;
    Mvol(irow, 0) = 1.0;
    for(int icol = 1; icol < m_ncol; icol++)
      {
        Mvol(irow, icol) = 0.0;
      }
  }


  ///
  void transpose();


  void checkConditionNumber() const;

  void checkUpperTriangularConditionNumber() const;

  ///
  /**
     sets product = a_left* a_right
     fails if a_left.m_col != a_right.m_rows
  */
  friend void multiply(LAPACKMatrix& a_product, 
                       const LAPACKMatrix& a_left,
                       const LAPACKMatrix& a_right);

  ///below stuff is shamelessly stolen from lapackwrapper class

  ///
  /**
     Solves A*X = B using general least squares, for each column of B
  */
  friend int solveLeastSquares(LAPACKMatrix& A, LAPACKMatrix& B);

  ///
  /**
     Solves A'*X = B using least squares, for vector b.   Answer goes 
     back into B I think
  */
  friend int solveLeastSquaresTranspose(LAPACKMatrix& A, LAPACKMatrix& B);

  ///
  /**
   *  Solves A^T X = B using least squares with SVD, for vector b
   */
  friend int solveLSTSVD(LAPACKMatrix& A, LAPACKMatrix& B, int a_maxiter, Real a_tol);

  ///
  /**
   *  Solves A*X = B using least squares with SVD, for X
   */
  friend int solveLSTSVD(LAPACKMatrix& X, const LAPACKMatrix& A, const LAPACKMatrix& B, int a_maxiter, Real a_tol);

  ///
  friend int solveLSTSVDOnce(LAPACKMatrix& X, LAPACKMatrix& B);
  
  /// 
  /**
     Solves A*X = B using least squares with SVD, for X
   */
  friend int solveLSTSVDOnce(LAPACKMatrix& X, const LAPACKMatrix& A, const LAPACKMatrix& B);

  ///
  /**
   *  Solves equality constrained least squares problem
   *    Find x, s.t. min norm(A x - c) with B x = d
   */
  friend int solveEqualityConstrainedLS(LAPACKMatrix& A, LAPACKMatrix& c, LAPACKMatrix& B, LAPACKMatrix& d, LAPACKMatrix& x);

  ///
  /**
   *  Solves A'*X = B using reduced rank least squares, for vector b
   */
  friend int solveReducedRankLS(LAPACKMatrix& A, LAPACKMatrix& b);

  ///
  /**
     Following Lapack, gets inverse of condition number.   Returning a number
     near zero means the matrix is not really solvable.
   */
  friend Real getInverseOfConditionNumber(const LAPACKMatrix& A);

  friend Real getInverseOfUpperTriangularConditionNumber(const LAPACKMatrix& A);
  
  ///  turn on if you want every solve to check the condition number
  static bool s_checkConditionNumber;

  ///
  static bool s_verbose;
  ///
  static bool s_outputStenData;

private:

  int   m_nrow;
  int   m_ncol;
  Real* m_data;
  bool  m_alias=false;
};
///
Real getInverseOfConditionNumber(const LAPACKMatrix& A);

///
/**
   sets product = a_left* a_right
   fails if a_left.m_col != a_right.m_rows
*/
void multiply(LAPACKMatrix& a_product, 
              const LAPACKMatrix& a_left,
              const LAPACKMatrix& a_right);

///below stuff is shamelessly stolen from lapackwrapper class

///
/**
   Solves A*X = B using general least squares, for each column of B
*/
int solveLeastSquares(LAPACKMatrix& A, LAPACKMatrix& B);

///
/**
   Solves A'*X = B using least squares, for vector b
*/
int solveLeastSquaresTranspose(LAPACKMatrix& A, LAPACKMatrix& B);

///
int solveLSTSVDOnce(LAPACKMatrix& X, const LAPACKMatrix& A, const LAPACKMatrix& B);

///
/**
 *  Solves A^T X = B using least squares with SVD, for vector b
 */
int solveLSTSVD(LAPACKMatrix& A, LAPACKMatrix& B, int a_maxiter, Real a_tol);

///
int solveLSTSVDOnce(LAPACKMatrix& X, const LAPACKMatrix& A, const LAPACKMatrix& B);

///
/**
 *  Solves equality constrained least squares problem
 *    Find x, s.t. min norm(A x - c) with B x = d
 */
int solveEqualityConstrainedLS(LAPACKMatrix& A, LAPACKMatrix& c, LAPACKMatrix& B, LAPACKMatrix& d, LAPACKMatrix& x);

///
/**
 *  Solves A'*X = B using reduced rank least squares, for vector b
 */
int solveReducedRankLS(LAPACKMatrix& A, LAPACKMatrix & b);


#include "BaseNamespaceFooter.H"
#endif
