#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTEMAPPEDSUM_H_
#define _COMPUTEMAPPEDSUM_H_

#include "REAL.H"
#include "Vector.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "Interval.H"
#include "CoordSys.H"
#include "NamespaceHeader.H"

/// Returns the volume-weighted sum (integral) of phi over all valid regions
/*
 */
Real computeSum(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Vector<CoordSys<FArrayBox,FluxBox>* >& a_coordSysPtr,
                const Interval&                       a_comps = Interval(0,0),
                const int&                            a_lBase = 0);

/// Returns the volume-weighted sum (integral) of phi over all valid regions and the volume
/*
 */
Real computeSum(Real&                                 a_volume,
                const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Vector<CoordSys<FArrayBox,FluxBox>* >& a_coordSysPtr,
                const Interval&                       a_comps = Interval(0,0),
                const int&                            a_lBase = 0);

/// Returns the volume-weighted sum (integral) of phi over the valid region
/*
  This is the single-level version. If a_global is true, then do any MPI
  reductions in this function.  If it's false, don't do any MPI stuff
  (assumption is that it will be done outside this function).  In serial,
  the value of the a_global parameter has no effect.
 */
Real computeSum(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGrids,
                const int&                  a_nRefFine,
                const CoordSys<FArrayBox,FluxBox>*    a_coordSysPtr,
                const Interval&             a_comps = Interval(0,0),
                const bool                  a_global= true);

/// Returns the volume-weighted sum (integral) of phi over the valid region and the volume
/*
  This is the single-level version. If a_global is true, then do any MPI
  reductions in this function.  If it's false, don't do any MPI stuff
  (assumption is that it will be done outside this function).  In serial,
  the value of the a_global parameter has no effect.
 */
Real computeSum(Real&                       a_volume,
                const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGrids,
                const int&                  a_nRefFine,
                const CoordSys<FArrayBox,FluxBox>*    a_coordSysPtr,
                const Interval&             a_comps = Interval(0,0),
                const bool                  a_global= true);

#include "NamespaceFooter.H"
#endif
