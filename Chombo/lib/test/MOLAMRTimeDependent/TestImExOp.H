#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TESTIMEXOP_H_
#define _TESTIMEXOP_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "TestSolnData.H"
#include "TestRhsData.H"

#include "NamespaceHeader.H"

// Exact solution is exp((cE + cI)*t)*(1+t)
// So that:
//   explicit op = (1/(1+t) + cE) * phi;
//   implicit op = cI * phi;

class TestImExOp
{
public:
  static const Real s_cE;
  static const Real s_cI;

  /// Default constructor
  TestImExOp();

  /// Destructor
  virtual ~TestImExOp();

  virtual void define(const TestSolnData&   a_state,
                      Real a_dt,
                      Real a_dtscale);

  virtual bool isDefined() const { return m_isDefined; }

  virtual void resetDt(Real a_dt);

  virtual void stage0stage1(int a_stage0,
                            int a_stage1);

  virtual void explicitOp(TestRhsData&             a_result,
                          Real                    a_time,
                          const TestSolnData&  a_state,
                          int                  a_stage);

  virtual void implicitOp(TestRhsData&             a_result,
                          Real                    a_time,
                          const TestSolnData&  a_state,
                          int                  a_stage);

  virtual void solve(TestSolnData&   a_soln,
                     const TestRhsData&   a_rhs,
                     Real               a_time,
                     int               a_stage);

protected:

  bool m_isDefined;

  Real m_dt;

  Real m_dtscale;

  int m_stage0;

  int m_stage1;

private:
  // Disallowed for all the usual reasons
  void operator=(const TestImExOp&);
  TestImExOp(const TestImExOp&);
};

#include "NamespaceFooter.H"

#endif
