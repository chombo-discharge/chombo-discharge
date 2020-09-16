#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*==============================================================================
 *
 * This test fits a cubic polynomial to a coarse grid and checks the
 * interpolation to a fine mesh.  With simply a linear stretching of the grid,
 * the results should be exact.  Accuracy and conservation is tested.
 *
 *============================================================================*/

#include <iomanip>
#include <limits>

#include "AMRLevel.H"
#include "DebugOut.H"

#include "CartesianCS.H"
#include "FourthOrderUtil.H"
#include "LevelGridMetrics.H"
#include "SingleBlockCSAdaptor.H"

#include "UsingNamespace.H"

//--Prototypes

int
mapped4OInterpTest();

void
parseTestOptions(int argc ,char* argv[]);

//--Global variables for handling output:

static const char *pgmname = "mapped4OInterpTest";
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

/*--------------------------------------------------------------------*
 *  Entry
 *--------------------------------------------------------------------*/

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl ;

  int eekflag = mapped4OInterpTest();
  if (eekflag == 0)
    {
      pout() << indent << pgmname
             << ": test passed." << std::endl;
    }
  else
    {
      pout() << indent << pgmname
             << ": test FAILED with error code "
             << eekflag << std::endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}

/*--------------------------------------------------------------------*
 * Routines for comparing floating point numbers
 *--------------------------------------------------------------------*/

// Comparison with limit tol^2 as x and y -> 0.
// Returns true if not equal
template <typename T>
inline bool compare(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
  return std::fabs(x - y) >
    (std::min(std::fabs(x), std::fabs(y)) + tol)*tol;
}

// Comparison with limit tol as x and y -> 0.
// Return true if not equal
template <typename T>
inline bool compare1(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
  return std::fabs(x - y) >
    std::min(std::fabs(x), std::fabs(y))*tol + tol;
}

/*--------------------------------------------------------------------*
 * Test polynomials
 *--------------------------------------------------------------------*/

#if (CH_SPACEDIM == 2)
Real test_polynomial_4thO(const Real *const a, const RealVect &p)
{
  const Real x = p[0];
  const Real y = p[1];
  const Real x2 = x*x;
  const Real y2 = y*y;
  const Real x3 = x2*x;
  const Real y3 = y2*y;
  return (a[0] + a[1]*x + a[2]*y +
          a[3]*x2 + a[4]*x*y + a[5]*y2 +
          a[6]*x3 + a[7]*x2*y + a[8]*x*y2 + a[9]*y3);
}
#else
Real test_polynomial_4thO(const Real *const a, const RealVect &p)
{
  const Real x = p[0];
  const Real y = p[1];
  const Real z = p[2];
  const Real x2 = x*x;
  const Real y2 = y*y;
  const Real z2 = z*z;
  const Real x3 = x2*x;
  const Real y3 = y2*y;
  const Real z3 = z2*z;
  return (a[0] + a[1]*x + a[2]*y + a[3]*z +
          a[4]*x2 + a[5]*y2 + a[6]*z2 + a[7]*x*y + a[8]*y*z + a[9]*z*x +
          a[10]*x3 + a[11]*y3 + a[12]*z3 +
          a[13]*x2*y + a[14]*x*y2 +
          a[15]*y2*z + a[16]*y*z2 +
          a[17]*z2*x + a[18]*z*x2 +
          a[19]*x*y*z);
}
#endif

/*--------------------------------------------------------------------*
 * Dummy AMR Level
 *--------------------------------------------------------------------*/

class AMRDummy : public AMRLevel
{
public:
  AMRDummy()
  {
  }

  ~AMRDummy()
  {
  }

  Real advance()
  {
    return 0.;
  }

  void postTimeStep()
  {
  }

  void tagCells(IntVectSet& a_tags)
  {
  }

  void tagCellsInit(IntVectSet& a_tags)
  {
  }

  void regrid(const Vector<Box>& a_new_grids)
  {
  }

  void initialGrid(const Vector<Box>& a_new_grids)
  {
  }

  void initialData()
  {
  }

  void postInitialize()
  {
  }

  Real computeDt()
  {
    return 0.;
  }

  Real computeInitialDt()
  {
    return 0.;
  }

  void writeCheckpointHeader (HDF5Handle& a_handle) const
  {
  }

  void writeCheckpointLevel (HDF5Handle& a_handle) const
  {
  }

  void readCheckpointHeader (HDF5Handle& a_handle)
  {
  }

  void readCheckpointLevel (HDF5Handle& a_handle)
  {
  }

  void writePlotHeader (HDF5Handle& a_handle) const
  {
  }

  void writePlotLevel (HDF5Handle& a_handle) const
  {
  }
};

/*--------------------------------------------------------------------*
 * Fourth-order mapped interpolation testing
 *--------------------------------------------------------------------*/

int mapped4OInterpTest()
{
  int status = 0;

#if CH_SPACEDIM >= 2
  // Number of significant digits for comparing the interpolated solution with
  // the analytical solution
#ifdef CH_USE_DOUBLE
  const int interpTol = std::numeric_limits<Real>::digits10 - 4;
  const int sumTol = std::numeric_limits<Real>::digits10 - 4;
#else
  // petermc, 14 Jun 2012, changed tolerances from digits10-2 to digits10-3.
  const int interpTol = std::numeric_limits<Real>::digits10 - 3;
  const int sumTol = std::numeric_limits<Real>::digits10 - 3;
#endif

  // Coefficients of the test polynomial
#if (CH_SPACEDIM == 2)
  const Real a0[10] =
  {
    1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
  };
  const Real a1[10] =
  {
    2., 2., 2., 2., 2., 2., 2., 2., 2., 4.
  };
#else
  const Real a0[20] =
  {
    1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
    1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
  };
  const Real a1[20] =
  {
    2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
    4., 2., 2., 2., 2., 2., 2., 2., 2., 2.
  };
#endif

  // Coordinate system (metrics)
  RealVect origin(RealVect::Zero);
  RealVect stretch(RealVect(D_DECL(2., 0.5, 0.5)));
  // Multiblock adaptor (takes care of deleting the single block CS)
  MultiBlockCoordSysFactory *coordSysFact =
    new SingleBlockCSAdaptorFactory(new CartesianCSFactory(origin, stretch));

  // Problem parameters
  const int nComp = 2;
  const int nGhost = 1;
  IntVect ghostVect(nGhost*IntVect::Unit);

//--1st level

  ProblemDomain probDom1(Box(IntVect(D_DECL(0, 0, 0)),
                             IntVect(D_DECL(4, 4, 4))));
  AMRDummy amrLevel1;
  amrLevel1.define(NULL, probDom1, 0, 4);
  Real dx1 = 1./5.;
  RealVect dx1Vect(D_DECL(dx1, dx1, dx1));
  Vector<int> vproc(1);
  vproc[0] = 0;
  Vector<Box> vbox(1);
  vbox[0].define(IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(4, 4, 4)));
  DisjointBoxLayout dbl1(vbox, vproc);
  LevelData<FArrayBox> ld1U(dbl1, nComp, ghostVect);
  LevelData<FArrayBox> ld1JU(dbl1, nComp, ghostVect);
  LevelGridMetrics lgm1(nComp, 4);
  lgm1.define(&amrLevel1, coordSysFact, 0, dx1Vect, ghostVect);
  lgm1.postInitialGrid(&dbl1);
  DataIterator dit1 = dbl1.dataIterator();
  int i1 = 0;
  RealVect XiOffset = 0.5*dx1Vect;
  for (dit1.begin(); dit1.ok(); ++dit1)
    {
      FArrayBox& fabU = ld1U[dit1];
      Box box = fabU.box();  // With ghosts
      FArrayBox pntXi(box, SpaceDim);
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
          RealVect loc;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              // Get mapped coordinate
              loc[dir] = iv[dir]*dx1Vect[dir] + XiOffset[dir];
              pntXi(iv, dir) = loc[dir];
            }
          // Get pointwise U
          fabU(iv, 0) = test_polynomial_4thO(a0, loc);
          fabU(iv, 1) = test_polynomial_4thO(a1, loc);
//           fabU(iv, 0) = 1.;
//           fabU(iv, 1) = 2.;
        }
      // Get pointwise J
      FArrayBox pntJ(box, 1);
      const NewCoordSys *coordSys = lgm1.getCoordSys(dbl1[dit1]);
      coordSys->pointwiseJ(pntJ, pntXi, box);
      // Multiply the two to get JU
      // See AMRLevelVolTest::initialData
      FArrayBox &fabJU = ld1JU[dit1];
      fabJU.copy(fabU);
      for (int iComp = 0; iComp != nComp; ++iComp)
        {
          fabJU.mult(pntJ, box, 0, iComp, 1);
        }
      ++i1;
    }
  // Obtain <JU> - but this is only on the valid cells
  fourthOrderAverage(ld1JU, probDom1);
  // Need <JU> on a +1 ring of cells so we can compute <U>
  for (dit1.begin(); dit1.ok(); ++dit1)
    {
      Box interiorBox = lgm1.getBoxes()[dit1];
      // Extrapolate <JU> at domain boundaries
      if (!probDom1.domainBox().contains(grow(interiorBox, 1)))
        {
          FArrayBox& JU = ld1JU[dit1];
          secondOrderCellExtrapAtDomainBdry(JU,
                                            interiorBox,
                                            probDom1);
        }
    }
//   dit1.begin();
//   pout() << "<JU0> on coarse\n";
//   dumpFAB2DSlicePretty(&(ld1JU[dit1]), 0,
//                        IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(4, 4, 0)));
//   pout() << "<JU1> on coarse\n";
//   dumpFAB2DSlicePretty(&(ld1JU[dit1]), 1,
//                        IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(4, 4, 0)));
  // Calculate <U>
  for (dit1.begin(); dit1.ok(); ++dit1)
    {
      Box interiorBox = lgm1.getBoxes()[dit1];
      cellFGToCellF(ld1U[dit1],
                    ld1JU[dit1],
                    lgm1.m_J[dit1],
                    interiorBox);
    }
//   dit1.begin();
//   pout() << "<U0> on coarse\n";
//   dumpFAB2DSlicePretty(&(ld1U[dit1]), 0,
//                        IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(4, 4, 0)));
//   pout << "<U1> on coarse\n";
//   dumpFAB2DSlicePretty(&(ld1U[dit1]), 1,
//                        IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(4, 4, 0)));

//--2nd level (*x is the exact solution)

  ProblemDomain probDom2(Box(IntVect(D_DECL(0, 0, 0)),
                             IntVect(D_DECL(19, 19, 19))));
  AMRDummy amrLevel2;
  amrLevel2.define(&amrLevel1, probDom2, 1, 1);
  Real dx2 = 1./20.;
  RealVect dx2Vect(D_DECL(dx2, dx2, dx2));
  vbox[0].define(IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(4, 4, 4)));
//   vbox[1].define(IntVect(6, 6), IntVect(7, 9));
  DisjointBoxLayout dbl2c(vbox, vproc);
  DisjointBoxLayout dbl2;
  refine(dbl2, dbl2c, 4);
  LevelData<FArrayBox> ld2U(dbl2, nComp, ghostVect);
  LevelData<FArrayBox> ld2JU(dbl2, nComp, ghostVect);
  LevelData<FArrayBox> ld2JUx(dbl2, nComp, ghostVect);  // Exact for comparison
  LevelGridMetrics lgm2(nComp, 4);
  lgm2.define(&amrLevel2, coordSysFact, &lgm1, dx2Vect, ghostVect);
  lgm2.postInitialGrid(&dbl2);
  DataIterator dit2 = dbl2.dataIterator();
  int i2 = 0;
  XiOffset = 0.5*dx2Vect;
  for (dit2.begin(); dit2.ok(); ++dit2)
    {
      FArrayBox& fabU = ld2U[dit2];
      Box box = fabU.box();  // With ghosts
      FArrayBox pntXi(box, SpaceDim);
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
          RealVect loc;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              // Get mapped coordinate
              loc[dir] = iv[dir]*dx2Vect[dir] + XiOffset[dir];
              pntXi(iv, dir) = loc[dir];
            }
          // Get pointwise U
          fabU(iv, 0) = test_polynomial_4thO(a0, loc);
          fabU(iv, 1) = test_polynomial_4thO(a1, loc);
//           fabU(iv, 0) = 1.;
//           fabU(iv, 1) = 2.;
        }
      // Get pointwise J
      FArrayBox pntJ(box, 1);
      const NewCoordSys *coordSys = lgm2.getCoordSys(dbl2[dit2]);
      coordSys->pointwiseJ(pntJ, pntXi, box);
      // Multiply the two to get JU
      // See AMRLevelVolTest::initialData
      FArrayBox &fabJUx = ld2JUx[dit2];
      fabJUx.copy(fabU);
      for (int iComp = 0; iComp != nComp; ++iComp)
        {
          fabJUx.mult(pntJ, box, 0, iComp, 1);
        }
      // Default for the computed value
      FArrayBox &fabJU = ld2JU[dit2];
      fabJU.setVal(0.);
      ++i2;
    }
  // Obtain <JU> - (only) on the valid cells
  fourthOrderAverage(ld2JUx, probDom2);

//--Interpolate

  RealVect dx(IntVect::Unit);
  IntVect refRatio(4*IntVect::Unit);
  FourthOrderMappedFineInterp stencil(3, nComp, refRatio, dx);
  dit1.begin();
  dit2.begin();
  Interval dummyIntv;  // Only used for MB
  for (int i = 0; i != 1; ++i)
    {
      stencil.interpToFine(ld2JU,
                           lgm2.m_J,
                           ld1U,
                           ld1JU,
                           &lgm1,
                           dummyIntv);
    }

//--Check the results

  // Loop over the coarse cells and the fine cells therein
  dit2.begin();
  for (dit1.begin(); dit1.ok(); ++dit1)
    {
      pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
      pout().precision(std::numeric_limits<Real>::digits10);
      const int width = std::numeric_limits<Real>::digits10 + 7;
      FArrayBox &fab1JU  = ld1JU[dit1];
      FArrayBox &fab2JUx = ld2JUx[dit2];
      FArrayBox &fab2JU  = ld2JU[dit2];
      BoxIterator bit1(dbl1[dit1]);
      for (bit1.begin(); bit1.ok(); ++bit1)
        {
          bool faultInCoarseCell = false;
          const IntVect coarseIV = bit1();
          // Fine box in a coarse cell
          Box boxFnInCr(coarseIV*refRatio, coarseIV*refRatio + refRatio - 1);
          BoxIterator bit2(boxFnInCr);
          for (int iComp = 0; iComp != nComp; ++iComp)
            {
              bool faultInComponent = false;
//               pout() << "--> For component " << iComp << std::endl;
//               pout() << "Coarse JU: " << fab1JU(coarseIV, iComp)
//                      << " ... at: " << coarseIV << std::endl;
              Real sumFine = 0.;
              for (bit2.begin(); bit2.ok(); ++bit2)
                {
                  const IntVect iv = bit2();
                  if (compare(fab2JUx(iv, iComp), fab2JU(iv, iComp), interpTol))
                    {
                      pout() << "Analytic and computed interpolations do "
                        "not compare at fine cell " << iv << std::endl;
                      pout() << "Exact: " << std::setw(width)
                             << fab2JUx(iv, iComp) << ", Computed: "
                             << std::setw(width) << fab2JU(iv, iComp)
                             << std::endl;
                      faultInCoarseCell = true;
                      faultInComponent = true;
                    }
                  sumFine += fab2JU(iv, iComp);
                }
              sumFine /= refRatio.product();
              if (compare(fab1JU(coarseIV, iComp), sumFine, sumTol))
                {
                  pout() << "Sum of averages in fine cells is not equal to "
                    "coarse average for component " << iComp
                         << ".  Interpolation was not conservative.\n";
                  pout() << "Coarse <JU>: "<< std::setw(width)
                         <<  fab1JU(coarseIV, iComp)
                         << ", fine weighted sum <JU>: " << std::setw(width)
                         << sumFine << std::endl;
                  faultInCoarseCell = true;
                  faultInComponent = true;
                }
//               pout() << "Sum of fine values: " << sumFine << std::endl;
              if (faultInComponent)
                {
                  pout() << "^^^ Failure for component " << iComp << std::endl;
                }
            }
          if (faultInCoarseCell)
            {
              pout() << "^^^ Failure at coarse cell " << coarseIV << std::endl;
              ++status;
            }
        }  // Loop over coarse cells
      pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      pout().precision(6);
      ++dit2;  // Assume fine boxes overlay coarse boxes
    }

//--Cleanup

  delete coordSysFact;
#endif
  if (status)
    {
      pout() << "Failure in " << status << " coarse cells\n";
    }
  return status;
}

/*--------------------------------------------------------------------*
 *  Parse the standard test options (-v -q) out of the command line.
 *--------------------------------------------------------------------*/

void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if ( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              // argv[i] = "" ;
            }
        }
    }
  return ;
}

