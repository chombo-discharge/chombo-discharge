#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _PARTICLEVALIDREGION_H_
#define _PARTICLEVALIDREGION_H_

#include "REAL.H"
#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
 #include "NamespaceHeader.H"

/// This class defines the region in which the particle data is valid.
/// The invalid region is a layer of a_buffer cells abutting the next
/// coarser level.  The buffer is computed in such a way that during
/// two successive synchronization points the particles
/// deposition/interpolation stencil remains fully contained on this level

class ParticleValidRegion
{
public:

  /// Default constructor.  User must subsequently call define().
  ParticleValidRegion()
    : m_isDefined(false) {}

  /// Destructor.
  ~ParticleValidRegion() {}

  /// Full constructor.
  ParticleValidRegion(const DisjointBoxLayout&         a_domain,
                      const LevelData<BaseFab<bool> >* a_crseValid,
                      const int                        a_refRatio,
                      const int                        a_buffer);

  /// Like the full constructor.
  void define(const DisjointBoxLayout&         a_domain,
              const LevelData<BaseFab<bool> >* a_crseValid,
              const int                        a_refRatio,
              const int                        a_buffer);

  /// Return pointer to valid region.
  const LevelData<BaseFab<bool> >* mask() const
  {

    CH_TIME("ParticleValidRegion::mask");

    if (m_isDefined) return &m_valid;
    else return NULL;
  }

protected:

  /// internal use only
  inline void bufferGhostJob(BaseFab<bool>& a_valid,
                             const Box&     a_box,
                             const Box&     a_domain,
                             const int      a_buffer,
                             const int      a_dir);
protected:

  bool m_isDefined;

  /// levelData mask of valid region
  LevelData<BaseFab<bool> > m_valid;

};

#include "NamespaceFooter.H"

#endif
