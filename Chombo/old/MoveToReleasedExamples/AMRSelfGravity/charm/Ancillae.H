#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// various simple functions useful for a variety of purposes

#ifndef _ANCILLAE_H_
#define _ANCILLAE_H_

#include "REAL.H"
#include "Vector.H"
#include "List.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Interval.H"
#include "LayoutIterator.H"
#include "LGintegrator.H"
#include "CONSTANTS_C.H"


///
Real globalAverage(const LevelData<FArrayBox>& a_levFab,
                   const int                   a_indx = 0,
                   const bool                  a_perVol = true);

/// need to pass a copy cause MPI_Allreduce does not guarantee it's constantness
Real globalMax(Real a_localVal);

/// need to pass a copy cause MPI_Allreduce does not guarantee it's constantness
Real globalMin(Real a_localVal);

/// return max value of a_comp component of this levelData
Real levelMax(const LevelData<FArrayBox>& a_levFab,
              const int a_comp);

///return min value of a_comp component of this levelData
Real levelMin(const LevelData<FArrayBox>& a_levFab,
              const int a_comp);

/// set to zero
inline void resetToZero(LevelData<FArrayBox>& a_lev);

/// set to val
template <class T>
inline void setToVal(LevelData<BaseFab<T> >& a_lev, const T& a_val);

/// offset a levelData by a given value
inline void offset(LevelData<FArrayBox>& a_lev, const Real& a_offset);

/// rescale a levelData by a given factor
inline void rescale(LevelData<FArrayBox>& a_lev, const Real& a_factor);

///
void computeDeltaU(const LevelData<FArrayBox>& a_UNew,
                   const LevelData<FArrayBox>& a_UOld);

///
void computePlotVars(LevelData<FArrayBox>&       a_UPlot,
                     const LevelData<FArrayBox>& a_UNew);

// returns total mass of particles in a list;
//inline Real massPartList(const List<Particle >& a_partList);

// returns mass of a particle on this level;
//inline Real levParticleMass(const LevelData<PBinFab >& a_levDM);

///
void bufferBoxes(DisjointBoxLayout&       a_bufferGrids,
                 const DisjointBoxLayout& a_grids,
                 const DisjointBoxLayout& a_crseGrids,
                 const int                a_bufferSize,
                 const int                a_refRatio);

/// Interpolate in time
void interpolateInTime(LevelData<FArrayBox>&       a_phi,
                       const LevelData<FArrayBox>& a_phiNew,
                       const LevelData<FArrayBox>& a_phiOld,
                       const Real&                 a_time,
                       const Real&                 a_tNew,
                       const Real&                 a_tOld,
                       const Real&                 a_dt,
                       const Interval&     a_srcInterval=Interval(0,0),
                       const Interval&     a_dstInterval=Interval(0,0),
                       const IntVect&      a_ghost      =IntVect::Zero);

/// set to zero
inline void resetToZero(LevelData<FArrayBox>& a_lev)
{
  for (DataIterator di = a_lev.dataIterator(); di.ok(); ++di)
  {
    a_lev[di].setVal(zero);
  }
}


/// set to val
template<class T>
inline void setToVal(LevelData<BaseFab<T> >& a_lev, const T& a_val)
{
  for (DataIterator dit = a_lev.dataIterator(); dit.ok(); ++dit)
  {
    a_lev[dit()].setVal(a_val);
  }
}

/// offset a levelData by a given value
inline void offset(LevelData<FArrayBox>& a_lev, const Real& a_offset)
{
  for (DataIterator dit = a_lev.dataIterator(); dit.ok(); ++dit)
  {
    a_lev[dit()] -= a_offset;
  }
}

/// rescale a levelData by a given factor
inline void rescale(LevelData<FArrayBox>& a_lev, const Real& a_factor)
{
  for (DataIterator dit = a_lev.dataIterator(); dit.ok(); ++dit)
  {
    a_lev[dit()] *= a_factor;
  }
}


#endif
