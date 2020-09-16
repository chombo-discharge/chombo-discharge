#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTEMAPPEDNORM_H_
#define _COMPUTEMAPPEDNORM_H_

#include "REAL.H"
#include "Vector.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Interval.H"
#include "FluxBox.H"
#include "CoordSys.H"
#include "NamespaceHeader.H"

/// Returns the volume-weighted norm of phi over all valid regions
Real computeNorm(const Vector<LevelData<FArrayBox>* >& a_phi,
                 const Vector<int>&                    a_nRefFine,
                 const Vector<CoordSys<FArrayBox,FluxBox>* >& a_coordSys,
                 const Interval                        a_comps = Interval(0,0),
                 const int                             a_p = 2,
                 const int                             a_lBase = 0);

/// Returns the volume-weighted norm of phi over the valid region
Real computeNorm(const LevelData<FArrayBox>& a_phi,
                 const DisjointBoxLayout*    a_finerGrids,
                 const int                   a_nRefFine,
                 const CoordSys<FArrayBox,FluxBox>* a_coordSys,
                 const Interval              a_comps = Interval(0,0),
                 const int                   a_p = 2);

/// Returns the max value of phi over all valid regions
Real computeMax(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Interval                        a_comps = Interval(0,0),
                const int                             a_lBase = 0);



#include "NamespaceFooter.H"
#endif
