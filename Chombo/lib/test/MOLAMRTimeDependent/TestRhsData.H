#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TESTRHSDATA_H_
#define _TESTRHSDATA_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelFluxRegister.H"
#include "RealVect.H"
#include "TestSolnData.H"

#include "NamespaceHeader.H"

// Forward declaration
class TestSolnData;

class TestRhsData
{
public:

  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  TestRhsData();

  /// Destructor
  virtual ~TestRhsData();

  /// Actual constructor.
  void define(/// box layout at this level
              const DisjointBoxLayout&    a_layout,
              /// number of conserved components
              const int                   a_nComp,
              /// ghost vector
              const IntVect&              a_ghostVect
              );

  void define(const TestSolnData& a_state);

  void aliasData(LevelData<FArrayBox>& a_data, 
      LevelData<FArrayBox>* a_accum = NULL);

  void copy(const TestRhsData& a_state);

  void copy(const TestSolnData& a_state);

  void zero();

  void increment(const TestRhsData& a_increment,
                 Real a_factor = 1.,
                 bool a_updateFluxReg = false);

  const LevelData<FArrayBox>& data() const
  { return m_data; }

  LevelData<FArrayBox>& data()
  { return m_data; }

protected:

  /// layout for this level
  DisjointBoxLayout m_grids;

  /// number of components
  int m_nComp;

  /// data on this level
  LevelData<FArrayBox> m_data;

  /// extra data to accumulate fluxes
  LevelData<FArrayBox>* m_accum;

  /// ghost vector
  IntVect m_ghostVect;

  /// whether has been defined
  bool m_isDefined;

private:

  // Disallowed for all the usual reasons
  void operator=(const TestRhsData&);
  TestRhsData(const TestRhsData&);
};

#include "NamespaceFooter.H"

#endif
