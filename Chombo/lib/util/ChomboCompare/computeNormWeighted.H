#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTE_NORM_WEIGHTED_H_
#define _COMPUTE_NORM_WEIGHTED_H_

#include "REAL.H"
#include "Vector.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Interval.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

Real normWeighted(const BoxLayoutData<FArrayBox>& a_layout,
                  const BoxLayoutData<FArrayBox>& a_J,
                  const Interval& a_intvl,
                  const int& p);

/// Returns the volume-weighted (including a_J) p-norm of phi over all valid regions
Real computeNormWeighted(const Vector<LevelData<FArrayBox>* >& a_phi,
                         const Vector<LevelData<FArrayBox>* >& a_J,
                         const Vector<int>&                    a_nRefFine,
                         const Real                            a_dxCrse,
                         const Interval                        a_comps = Interval(0,0),
                         const int                             a_p = 2,
                         const int                             a_lBase = 0);

/// Returns the volume-weighted (including a_J) p-norm of phi over the valid region
Real computeNormWeighted(const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_J,
                         const DisjointBoxLayout*    a_finerGrids,
                         const int                   a_nRefFine,
                         const Real                  a_dx,
                         const Interval              a_comps = Interval(0,0),
                         const int                   a_p = 2);

/// Returns the volume-weighted (including a_J) p-norm of phi over all valid regions: anisotropic version
Real computeNormWeighted(const Vector<LevelData<FArrayBox>* >& a_phi,
                         const Vector<LevelData<FArrayBox>* >& a_J,
                         const Vector<IntVect>&                a_nRefFine,
                         const RealVect&                       a_dxCrse,
                         const Interval                        a_comps = Interval(0,0),
                         const int                             a_p = 2,
                         const int                             a_lBase = 0);

/// Returns the volume-weighted (including a_J) p-norm of phi over the valid region: anisotropic version
Real computeNormWeighted(const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_J,
                         const DisjointBoxLayout*    a_finerGrids,
                         const IntVect&              a_nRefFine,
                         const RealVect&             a_dx,
                         const Interval              a_comps = Interval(0,0),
                         const int                   a_p = 2);

#include "NamespaceFooter.H"
#endif
