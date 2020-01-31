#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// generalFuncs.H
// petermc, 24 March 2003; modified 28 Jan 2009

#ifndef _GENERALFUNCS_H_
#define _GENERALFUNCS_H_
#include  <iostream>
#include  <cstdio>
#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "ParmParse.H"
#ifdef CH_USE_PYTHON
#include "PyParse.H"
#endif

#include "UsingNamespace.H"

namespace GenFuncs
{


/*==============================================================================
 *  Routines for comparing floating point numbers
 *============================================================================*/

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


/*==============================================================================
 *  Routines for generating a box layout
 *============================================================================*/

Box boxFromVector(const Vector<int>&  a_ints,
                  int                 a_start = 0);

int readBoxes(Vector< Vector<Box> >&     a_vectVectBoxes,
              ParmParse&                 a_pp,
              const ProblemDomain&       a_domain,
              const int&                 a_maxGridSize,
              const int&                 a_blockFactor,
              const int&                 a_numlevels,
              const Vector<int>&         a_vectRefRatio,
              const bool&                a_verbose);

#ifdef CH_USE_PYTHON
// Python-savvy version of readBoxes.
int readBoxes(Vector< Vector<Box> >&     a_vectVectBoxes,
              PyParse&                   a_pp,
              const ProblemDomain&       a_domain,
              const int&                 a_maxGridSize,
              const int&                 a_blockFactor,
              const int&                 a_numlevels,
              const Vector<int>&         a_vectRefRatio,
              const bool&                a_verbose);
#endif

int getDiff(Vector<LevelData<FArrayBox>* >& a_vectDiff,
            const Vector<LevelData<FArrayBox>* >& a_vectPhi1,
            const Vector<LevelData<FArrayBox>* >& a_vectPhi2);

int getMagnitude(Vector<LevelData<FArrayBox>* >& a_vectMag,
                 const Vector<LevelData<FArrayBox>* >& a_vectField);

int project2(Vector<LevelData<FArrayBox>* >& a_vectProj,
             const Vector<LevelData<FArrayBox>* >& a_vectPhi,
             const Vector<DisjointBoxLayout>& a_vectGrids);


int printErrorNorms(const string& a_prefix,
                    Real a_dxFine,
                    int a_deginterp,
                    Real a_normMaxerrFine,
                    Real a_norm1errFine,
                    Real a_norm2errFine,
                    Real a_normMaxerrCoarse,
                    Real a_norm1errCoarse,
                    Real a_norm2errCoarse);

int printDiffNorms(const string& a_prefix,
                   Real a_dxFine,
                   int a_deginterp,
                   Real a_normMaxdiff,
                   Real a_norm1diff,
                   Real a_norm2diff);


int getBoxesFromIndices(Vector<Box>&        a_vectBoxes,
                        const Box&          a_domainBox,
                        const int&          a_maxGridSize,
                        const int&          a_blockFactor,
                        const int&          a_boxCount,
                        const Vector<int>&  a_boxIndices,
                        const Vector<int>&  a_domainIndices);

}  // End of namespace GenFuncs

#endif
