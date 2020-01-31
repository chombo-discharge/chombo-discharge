#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include "parstream.H"
#include "LAPACKMatrix.H"
#include "Misc.H"
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "UsingNamespace.H"

#ifdef CH_USE_DOUBLE
Real g_tol = 1.0e-12;
#else
Real g_tol = 1.0e-6;
#endif
using std::endl;
//shamelessly swiped from Hans' Chombo4 test
int solveLeastSquaresTest()
{
  int size = 3;
  /* matrix A */
  Real cA[3*3] = { 3.1, 1.0 , 3.4,   
                     1.3,-6.9 , 7.2,    
                    -5.7, 5.8 ,-8.8}   ;

  /* matrix B */
  Real cB[3*3] = {-1.3, -0.1, 1.8,   
                    -1.2, -0.3, 1.9,    
                    -1.2, -0.2,  1.8}  ; 

  /* least squares solution from octave, A \ B */
  Real cC[3*3] = {1.,                1.,                1.,
                    0.941236852109999, 0.973843619313142, 0.944531745025979 ,
                    0.938461538461539, 0.953846153846154, 0.938461538461539};


  LAPACKMatrix A(size, size, cA);
  LAPACKMatrix B(size, size, cB);
  LAPACKMatrix C(size, size, cC);

  pout() << "least squares test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << "B = " << endl;
  B.poutAll();
  solveLeastSquares(A, B);
  pout() << "Answer = " << endl;
  B.poutAll();
  pout() << "correct Answer = " << endl;
  C.poutAll();
  

  // Check the answer
  C -= B;
  for(int irow = 0; irow < size; irow++)
    {
      for(int icol = 0; icol < size; icol++)
        {
          if(Abs(C(irow,icol)) > g_tol)
            {
              pout() << "at row = " << irow << ", icol = " << icol << endl;
              pout() << "least square test returned error of "<< C(irow,icol) << endl;
              return -7;
            }
        }
    }

  return 0;
}

int inverseTest()
{
  int n =4;
  LAPACKMatrix A(n, n);
  Real val = 1;
  for(int irow = 0; irow < n; irow++)
    {
      for(int icol = 0; icol < n; icol++)
        {
          bool odd = (irow%2 == 1);
          if(odd && (irow== icol))
            A(irow, icol) = -val;
          else
            A(irow, icol) =  val;
          val += 1;
        }
    }
  LAPACKMatrix Ainv = A;
  int test = Ainv.invert();
  if(test != 0)
    {
      pout() << "we started with a singular matrix" << endl;
      return test;
    }

  LAPACKMatrix AAinv;
  multiply(AAinv, A, Ainv);

  pout() << "inverse test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << endl;

  pout() << "Ainv = " << endl;
  Ainv.poutAll();
  pout() << endl;

  pout() << "A * Ainv = " << endl;
  AAinv.poutAll();
  pout() << endl;

  int retval = 0;
  for(int irow = 0; irow < n; irow++)
    {
      for(int icol = 0; icol < n; icol++)
        {
          Real corrval = 0;
          if(irow == icol)
            {
              corrval = 1.0;
            }
          if(Abs(AAinv(irow, icol)-corrval) > 1.0e-3) retval -= 1;
        }
    }
  return retval;
}
/**/
int inverseWithSVDTest()
{
  int n =4;
  int m = n;
  LAPACKMatrix A(m, n);
  Real val = 1;
  A.setVal(0.);
  for(int irow = 0; irow < m; irow++)
    {
      for(int icol = 0; icol < n; icol++)
        {
          bool odd = (irow%2 == 1);
          if(odd && (irow== icol))
            A(irow, icol) = -val;
          else
            A(irow, icol) =  val;
          val += 1;
        }
    }
  LAPACKMatrix Ainv = A;
  int test = Ainv.invertUsingSVD(1, 1.0e-6);
  if(test != 0)
    {
      pout() << "we started with a singular matrix" << endl;
      return test;
    }

//  A.truncate(n, n);
//  Ainv.truncate(n, n);

  LAPACKMatrix AAinv;
  Ainv.transpose();
  multiply(AAinv, A, Ainv);

  pout() << "inverse with SVD test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << endl;

  pout() << "Ainv = " << endl;
  Ainv.poutAll();
  pout() << endl;

  pout() << "A * Ainv = " << endl;
  AAinv.poutAll();
  pout() << endl;

  int retval = 0;
  for(int irow = 0; irow < n; irow++)
    {
      for(int icol = 0; icol < n; icol++)
        {
          Real corrval = 0;
          if(irow == icol)
            {
              corrval = 1.0;
            }
          if(Abs(AAinv(irow, icol)-corrval) > 1.0e-3) retval -= 1;
        }
    }
  return retval;
}

/**/
int inverseWithLeastSquaresTest()
{
  int n =4;
  LAPACKMatrix A(n, n);
  Real val = 1;
  A.setVal(0.);
  for(int irow = 0; irow < n; irow++)
    {
      for(int icol = 0; icol < n; icol++)
        {
          bool odd = (irow%2 == 1);
          if(odd && (irow== icol))
            A(irow, icol) = -val;
          else
            A(irow, icol) =  val;
          val += 1;
        }
    }
  LAPACKMatrix Ainv = A;
  int test = Ainv.invertUsingLeastSquares();
  if(test != 0)
    {
      pout() << "we started with a singular matrix" << endl;
      return test;
    }

  A.truncate(n, n);
  Ainv.truncate(n, n);

  LAPACKMatrix AAinv;
  multiply(AAinv, A, Ainv);

  pout() << "inverse with Least Squares test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << endl;

  pout() << "Ainv = " << endl;
  Ainv.poutAll();
  pout() << endl;

  pout() << "A * Ainv = " << endl;
  AAinv.poutAll();
  pout() << endl;

  int retval = 0;
  for(int irow = 0; irow < n; irow++)
    {
      for(int icol = 0; icol < n; icol++)
        {
          Real corrval = 0;
          if(irow == icol)
            {
              corrval = 1.0;
            }
          if(Abs(AAinv(irow, icol)-corrval) > 1.0e-3) retval -= 1;
        }
    }
  return retval;
}

/**/
int transposeTest()
{
  Real cA[3*3] = 
    {
      1,2,3,
      4,5,6,
      7,8,9
    };
  Real cAtrancorrect [3*3] = 
    {
      1,4,7,
      2,5,8,
      3,6,9
    };

  pout() << "transpose test" << endl;
  pout() << "A = " << endl;
  LAPACKMatrix A(3, 3, cA); 
  A.poutAll();
  
  LAPACKMatrix Atran = A;
  Atran.transpose();
  pout() << "A transpose = " << endl;
  Atran.poutAll();

  Real eps = 1.0e-6;
  for(int irow = 0; irow < 3; irow++)
    {
      for(int icol = 0; icol < 3; icol++)
        {
          Real compval = Atran(irow,icol);
          int offset   = Atran.offset(irow,icol);
          Real corrval = cAtrancorrect[offset];
          if(Abs(compval - corrval) > eps)
            {
              MayDay::Warning("tranposition error");
              return -7;
            }
        }
    }
  return 0;
}


int testReducedRankLS()
{
  int m = 4;
  int n = 3;

  // A, the overdetermined least squares matrix, is m x n, m >= n
  // (row major order, these are the columns
  Real cA[12] = { 1, 1, 0,
                    0, 1, 0,
                    1, 0, 1,
                    0, 0, 1 };
  Real cB[] = {1, 1, 1, 0};  //last zero implicit in old code
  Real cexactRRLS[] = {1, 0, 0};
  
  LAPACKMatrix A(m, n, cA);
  LAPACKMatrix B(m, 1, cB);
  LAPACKMatrix exac(n, 1, cexactRRLS);

  pout() << "reduced rank test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << "B = " << endl;
  B.poutAll();
  solveReducedRankLS(A, B);

  pout() << "answer = " << endl;
  B.poutAll();
  pout() << "right answer =" << endl;
  exac.poutAll()        ;
  for (int i=0; i < n; ++i) 
    {
      Real calc =    B(i,0);
      Real corr = exac(i,0);
      if(Abs(calc-corr) > g_tol)
        {
          pout() << "reduced rank ls error = "<< calc -corr << endl;
          return -7;
        }
    }
  return 0;
}


///
int testLSTransposeResid()
{
  const int m = 4;
  const int n = 3;

  // A, the underdetermined least squares matrix, is m x n, m < n
  // (row major order, these are the columns
  Real cA[12] = { 1, 1, 0, 0,
                   1, 0, 0, 0,
                   1, 0, 0, 0};
  Real cB[] = {0, 0, 1, 0};

  Real cExactLS[] = {.5, -.5, 0, 0};

  LAPACKMatrix A(m, n, cA);
  LAPACKMatrix B(m, 1, cB);
  LAPACKMatrix exactLS(m, 1, cExactLS);

  int retval = solveLSTSVDOnce(A, B);
  if(retval != 0) return retval;

  for (int i=0; i < m; ++i) 
    {
      if(Abs(B(i, 0) - exactLS(i,0)) > g_tol)
        {
          Real residerr = B(i, 0) - exactLS(i,0);
          pout() << "lst point error  = "<< residerr << endl;
          return -7;
        }
    }

  return 0;
}

/**/
int testLSTransposeCOF()
{
  int m = 4;
  int n = 3;

  // A, the underdetermined least squares matrix, is m x n, m < n

  Real cA[12] = { 1, 1, 0, 0,
                         1, 0, 1, 0,
                         1, 0, 0, 1};
  Real cB[] = {1, 1, 1, 0};

  Real cexactRRLS[] = {.75, .25, .25, .25};

  LAPACKMatrix A(m, n, cA);
  LAPACKMatrix B(m, 1, cB);
  LAPACKMatrix exactRRLS(m, 1, cexactRRLS);

  solveLSTSVDOnce(A, B);

  for (int i=0; i < m; ++i) 
    {
      Real calc = B(i,0);
      Real corr = exactRRLS(i,0);
      if(Abs(calc-corr) > g_tol)
        {
          pout() << "lst transpose  cof error = " << calc-corr << endl;
          return -7;
        }
    }
  return 0;
}
/**/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  LAPACKMatrix::s_checkConditionNumber = true;

  int icode = 0; // to be returned, sum of return values (negative or zero)

  int retval = 0;
  retval  = testLSTransposeResid();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: testLSTTransposeResid test returned with value = " << retval << endl;
    }
  else
    {
      pout() << "LSTTransposeResid test passed" << endl;
    }

  retval = testLSTransposeCOF();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: testLSTransposecof test returned with value = " << retval << endl;
    }
  else
    {
      pout() << "lst transpose cof test passed" << endl;
    }


  /**/
  retval = inverseTest();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: inverse test returned with value = " << retval<< endl;
    }
  else
    {
      pout() << "inverse test passed" << endl;
    }
  /**/

  retval = 0; 
  retval = solveLeastSquaresTest();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: leastSquares test returned with value = " << retval << endl;
    }
  else
    {
      pout() << "least squares test passed" << endl;
    }
  /**/

  retval = testReducedRankLS();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: reduced rank ls test returned with value = " << retval << endl;
    }
  else
    {
      pout() << "reduced rank ls test passed" << endl;
    }

  retval = inverseWithLeastSquaresTest();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: inverse with least squares test returned with value = " << retval<< endl;
    }
  else
    {
      pout() << "inverse with least squares test passed" << endl;
    }
  /**/

  /**/
  retval = inverseWithSVDTest();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: inverse with SVD test returned with value = " << retval<< endl;
    }
  else
    {
      pout() << "inverse with SVD test passed" << endl;
    }
  /**/

  /**/
  retval = 0;
  retval = transposeTest();
  icode += retval;
  if(retval != 0)
    {
      pout() << "Error: transpose test returned with value = " << retval << endl;
    }
  else
    {
      pout() << "transpose test passed" << endl;
    }
  /**/

  if(icode == 0)
    {
      pout() << "all tests passed" << endl;
    }
  else
    {
      pout() << "not all tests passed, returned total value " << icode << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}
