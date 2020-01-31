#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::ifstream;
using std::ios;

#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "TestsIBC.H"

#include "TestsIBCF_F.H"
#include "SelfGravF_F.H"
#include "CONSTANTS_C.H"

// Null constructor
TestsIBC::TestsIBC()
{
    m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines
TestsIBC::TestsIBC(const std::string a_problem,
                   const Real&       a_gamma,
                   const Real&       a_ms,
                   const RealVect&   a_center,
                   const Real&       a_size,
                   const RealVect&   a_velocity,
                   const Real&      a_artvisc,
                   const Real&      a_rsTolerance,
                   const Real&      a_maxRsIter,
                   const Real&      a_maxMach)
{
    setFortranCommon(a_gamma,
                     a_ms,
                     a_center,
                     a_size,
                     a_velocity,
                     a_artvisc,
                     a_rsTolerance,
                     a_maxRsIter,
                     a_maxMach);

    setTestProblem(a_problem);
}

// Sets parameters in a common block used by Fortran routines:
/*
     a_gamma    - Gamma for polytropic, gamma-law gas
     a_ms       - Mach shock number of the discontinuity
     a_center   - Center of the explosion
     a_size     - Initial radius of the explosion
     a_velocity - Initial velocity of the gas
     a_artvisc  - Artificial viscosity coefficient

     a_rsTolerance - tolerance for the Riemann solver
     a_max_Rs_iter - max number of iteration in the Riemann solver
     a_Sthreshold  - min Eth/Etot to switch from Etot to entropy (S) eq.
*/
void TestsIBC::setFortranCommon(const Real&       a_gamma,
                                const Real&       a_ms,
                                const RealVect&   a_center,
                                const Real&       a_size,
                                const RealVect&  a_velocity,
                                const Real&      a_artvisc,
                                const Real&      a_rsTolerance,
                                const Real&      a_maxRsIter,
                                const Real&      a_maxMach)
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_SETSELFGRAV(CHF_CONST_REAL(a_gamma),
                   CHF_CONST_REAL(a_artvisc),
                   CHF_CONST_REAL(a_rsTolerance),
                   CHF_CONST_REAL(a_maxRsIter),
                   CHF_CONST_REAL(a_maxMach));

  FORT_SETTESTS(CHF_CONST_REAL(a_ms),
                CHF_CONST_REALVECT(a_center),
                CHF_CONST_REAL(a_size),
                CHF_CONST_REALVECT(a_velocity),
                CHF_CONST_REAL(a_artvisc));

  const Real small_r0 = small* a_gamma;
  const Real small_p0 = small* one;
  const Real small_u0 = small* a_ms;
  const Real small_s0 = small* exp(-a_gamma*log(a_gamma));

  FORT_SETFLOORS(CHF_CONST_REAL(small_r0),
                 CHF_CONST_REAL(small_u0),
                 CHF_CONST_REAL(small_p0),
                 CHF_CONST_REAL(small_s0));

  m_isFortranCommonSet = true;
}

// set the specific test problem
void TestsIBC::setTestProblem(const int& a_problem)
{
  m_problem=a_problem;
  m_isProblemSet = true;
}

//
void TestsIBC::setTestProblem(const std::string a_problem)
{
  //
  pout() << " test problem " << a_problem << '\n';

  if (a_problem == "explosion")
  {
    m_problem = PROBLEM_EXPLOSION;
  }
  else if (a_problem == "windcollision")
  {
    m_problem = PROBLEM_WINDCOLLISION;
  }
  else if (a_problem == "expandingshock")
  {
    m_problem = PROBLEM_EXPANDINGSHOCK;
  }
  else if (a_problem == "jeans")
  {
    m_problem = PROBLEM_JEANS;
  }
  else
  {
    pout() << "Invalid problem, \"" << a_problem <<
      "\", specified in input file" << '\n' << '\n';
    return;
  }

  m_isProblemSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_physIBC() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void TestsIBC::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* TestsIBC::new_physIBC()
{
  TestsIBC* retval = new TestsIBC();

  if (m_isFortranCommonSet == true)
  {
    retval->setFortranCommonSet();
    retval->setTestProblem(m_problem);
  }

  return static_cast<PhysIBC*>(retval);
}

// Set boundary fluxes
void TestsIBC::primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(PhysIBC::m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!PhysIBC::m_domain.isPeriodic(a_dir))
  {
    int lohisign;
    Box WGdnvBox = a_WGdnv.box();
    Box tmp = WGdnvBox;

    // Determine which side and thus shifting directions
    lohisign = sign(a_side);
    tmp.shiftHalf(a_dir,lohisign);

    // Is there a domain boundary next to this grid
    if (!PhysIBC::m_domain.contains(tmp))
    {
      tmp &= PhysIBC::m_domain;

      Box boundaryBox;

      // Find the strip of cells next to the domain boundary
      if (a_side == Side::Lo)
      {
        boundaryBox = bdryLo(tmp,a_dir);
      }
      else
      {
        boundaryBox = bdryHi(tmp,a_dir);
      }

      // Set the boundary fluxes

      /*
      FORT_SOLIDBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_Wextrap),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
      */

      // Set the boundary fluxes
      FORT_INFLOWF(CHF_FRA(a_WGdnv),
                   CHF_CONST_FRA(a_Wextrap),
                   CHF_CONST_FRA(a_W),
                   CHF_CONST_INT(lohisign),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(boundaryBox));
    }
  }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void TestsIBC::setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(PhysIBC::m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!PhysIBC::m_domain.isPeriodic(a_dir))
  {
    Box loBox,hiBox,centerBox,domain;
    int hasLo,hasHi;
    Box slopeBox = a_dW.box()&PhysIBC::m_domain;

    // Generate the domain boundary boxes, loBox and hiBox, if there are
    // domain boundarys there
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
               slopeBox,PhysIBC::m_domain,a_dir);

    // Set the boundary slopes if necessary
    if ((hasLo != 0) || (hasHi != 0))
    {
      FORT_SLOPEBCSF(CHF_FRA(a_dW),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi));
    }
  }
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void TestsIBC::artViscBC(FArrayBox&       a_F,
                         const FArrayBox& a_U,
                         const FArrayBox& a_divVel,
                         const int&       a_dir,
                         const Real&      a_time)
{
}

// Set up initial conditions
void TestsIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isProblemSet == true);
  CH_assert(PhysIBC::m_isDefined == true);
  const Real dx = PhysIBC::m_dx;

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= PhysIBC::m_domain;

    switch (m_problem)
    {
    case PROBLEM_EXPLOSION:
      {
        FORT_EXPLOSIONINITC(CHF_FRA(U),
                            CHF_CONST_REAL(dx),
                            CHF_BOX(uBox));
        break;
      }
    case PROBLEM_WINDCOLLISION:
      {
        FORT_WINDCOLLISION(CHF_FRA(U),
                           CHF_CONST_REAL(dx),
                           CHF_BOX(uBox));
        break;
      }
    case PROBLEM_EXPANDINGSHOCK:
      {
        FORT_EXPANDINGSHOCK(CHF_FRA(U),
                            CHF_CONST_REAL(dx),
                            CHF_BOX(uBox));
        break;
      }
    case PROBLEM_JEANS:
      {
        break;
      }
    }
  }
}

