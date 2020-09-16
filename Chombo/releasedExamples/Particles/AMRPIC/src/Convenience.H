#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONVENIENCE_H_
#define _CONVENIENCE_H_

#include "LevelData.H"
#include "JointParticle.H"
#include "ParticleData.H"

#include <map>
using std::map;

#include "NamespaceHeader.H"

/// Sort JointParticles into bins based on the mesh spacing. 
template <class P>
void binmapParticles(std::map<IntVect,JointParticle,CompIntVect>& a_mip,
                     const List<P>&                               a_part,
                     const RealVect&                              a_meshSpacing,
                     const RealVect&                              a_origin)
{
  CH_TIMERS("binmapParticles");

  for (ListIterator<P> li(a_part); li; ++li)
    {
      const P& p = li();
      const IntVect bin = locateBin(p.position(), a_meshSpacing, a_origin);

      map<IntVect,JointParticle,CompIntVect>::iterator it = a_mip.find(bin);

      if (it == a_mip.end())
        {
          a_mip[bin] = JointParticle(p.mass(),p.position(),p.velocity(),1);
        }
      else
        {
          it->second.addParticle( &p );
        }
    }
}

/// Set to val defined for an entire BoxLayoutData
template<class T, class V>
inline void setToVal(BoxLayoutData<T>& a_lt, const V a_v)
{
  for (DataIterator di=a_lt.dataIterator();di.ok();++di)
    {
      a_lt[di].setVal(a_v);
    }
}

/// Offset a BoxLayoutData by a given value
template <class T>
inline void offset(BoxLayoutData<T>& a_lt, const Real a_offset)
{
  for (DataIterator di=a_lt.dataIterator(); di.ok(); ++di)
    {
      a_lt[di] -= a_offset;
    }
}

/// Interpolate the data in time, given two endpoints
void interpolateInTime(LevelData<FArrayBox>&       a_phi,
                       const LevelData<FArrayBox>& a_phiNew,
                       const LevelData<FArrayBox>& a_phiOld,
                       const Real&                 a_time,
                       const Real&                 a_tNew,
                       const Real&                 a_tOld,
                       const Real&                 a_dt,
                       const Interval&             a_srcInterval=Interval(0,0),
                       const Interval&             a_dstInterval=Interval(0,0),
                       const IntVect&              a_ghost=IntVect::Zero);

#include "NamespaceFooter.H"

#endif
