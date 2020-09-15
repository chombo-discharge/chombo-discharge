#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PARTICLEDATA_H_
#define _PARTICLEDATA_H_

#include <map>
using std::map;

#include "BaseFab.H"
#include "RealVect.H"
#include "List.H"
#include "Box.H"
#include "LevelData.H"
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "NamespaceHeader.H"

struct boxids
{
  boxids(const unsigned id, const unsigned pi)
    : idx(id), pid(pi) {}

  boxids(const boxids& bi)
    : idx(bi.idx), pid(bi.pid) {}

  boxids()
    : idx(0), pid(0) {}

  void operator=(const boxids& bi)
  {
    idx=bi.idx;
    pid=bi.pid;
  }
  unsigned idx;
  unsigned pid;
};

/// compute the cell index containing the physical position of the item
inline IntVect locateBin(const RealVect a_x, const RealVect a_dx, const RealVect a_origin)
{
  CH_TIME("locateBin");
  const IntVect bin(D_DECL6((int)floor((a_x[0] - a_origin[0])/a_dx[0]),
			    (int)floor((a_x[1] - a_origin[1])/a_dx[1]),
			    (int)floor((a_x[2] - a_origin[2])/a_dx[2]),
			    (int)floor((a_x[3] - a_origin[3])/a_dx[3]),
			    (int)floor((a_x[4] - a_origin[4])/a_dx[4]),
			    (int)floor((a_x[5] - a_origin[5])/a_dx[5])));
  return bin;
}

///
/** This is a class for holding and sorting particle items.  Class <T>
    must have a \p RealVect \b <T>::position() const function which is
    used to assign the item in a list.  See class \p BinItem for a
    base class for particles.
*/
///
template<class P>
class ParticleData : public LayoutData< ListBox<P> >
{
public:

  /// Weak Constructor
  ParticleData();

  /// Full Constructor. Arguments:
  /// 
  /// a_dp:           The BoxLayout on which the particle data will be defined.
  ///                 Note that this *must* be a 'fixed size' layout.
  /// a_fixedBoxSize: The (fixed) box size in cells.
  /// a_domain:       the problem domain, used for enforcing periodic boundary conditions.
  /// a_meshSpacing:  the physical mesh spacing, used for binning the particles.
  /// a_origin:       the origin of the coordinate system, also for particle binning.
  ParticleData(const BoxLayout&         a_dp,
               const ProblemDomain&     a_domain,
	       const int&               a_fixedBoxSize,
	       const RealVect&          a_meshSpacing,
	       const RealVect&          a_origin);

  /// Define function. Same as the full constructor
  void define(const BoxLayout&         a_dp,
	      const ProblemDomain&     a_domain,
	      const int&               a_fixedBoxSize,
	      const RealVect&          a_meshSpacing,
	      const RealVect&          a_origin);

  /// Get the BoxLayout on which this ParticleData<P> is defined.
  const BoxLayout& getBoxes() const;

  /// Removes all the items from all the boxes in the container
  void clear();

  /// Return the number of particles in the container, outcast and valid.
  //  Note that this does not do a parallel reduce operation; each process
  /// returns the total number of particles in all the boxes it owns.
  size_t numParticlesLocal() const;

  /// Return the number of particles in the container, outcast and valid.
  /// This version does a parallel reduce to count all the particles on
  /// all processes. 
  size_t numParticles() const;

  /// Return the number of valid particles in the container, not counting outcasts.
  //  Note that this does not a do parallel reduce operation; each process
  /// returns the total number of particles in all the boxes it owns.
  size_t numValidLocal() const;

  /// Return the number of valid particles in the container, not counting outcasts.
  /// This version does a parallel reduce to count all the particles on
  /// all processes. 
  size_t numValid() const;

  /// Return the number of particles in the outcast list.
  //  Note that this does not a do parallel reduce operation; each process
  /// returns the total number of outcasts it owns.
  size_t numOutcastLocal() const;

  /// Return the number of particles in the outcast list.
  /// This version does a parallel reduce to count all the outcasts on
  /// all processes. 
  size_t numOutcast() const;
  
  /// Return a reference to the list of "outcast" particles - a list with
  /// particles not currently associated with any box in the container.
  List<P>& outcast();

  /// Is this ParticleData<P> closed - i.e. is the outcast list empty?
  bool isClosed() const;
  
  /// Get the physical domain associated with this ParticleData.
  const ProblemDomain& physDomain() const;

  /// Get the mesh spacing associated with this ParticleData,
  const RealVect& meshSpacing() const;

  /// Get the origin of the coordinate system 
  const RealVect& origin() const;

  /// Get the fixed Box size
  const int& fixedBoxSize() const;

  /// After the particles have moved, they might not still belong in the same ListBox
  /// they started in. This function collects all the particles that are no longer in
  /// the right box and puts them in the outcast list.
  void gatherOutcast();

  /// Redistribute the particles in the outcast list onto the proper
  /// box and process.
  void remapOutcast();

  /// Routine for filling ghost particles. This function will augment each
  /// ListBox in a_particlesWithGhosts with copies of all the particles that
  /// lie within a_numGhost cells of the box boundaries. If the problem domain 
  /// used to define this ParticleData is periodic, then ghost particles will be
  /// grabbed from the opposite sides of the domain as well. Ghost particles that are 
  /// copied from accross periodic boundaries in this fashion will have their position
  /// coordinates periodically shifted such that computing the distance between them
  /// and the "valid" particles will return the *shortest* distance defined on the 
  /// periodic domain.
  ///
  /// Arguments:
  /// a_particlesWithGhost: this is a ParticleData built on a BoxLayout like this one, 
  ///                       except that each box has been grown by a_numGhosts cells.
  ///                       On exit, each ListBox in a_particlesWithGhost will contain
  ///                       both "valid" and "ghost" particls.
  /// a_numGhost            the number of ghost cells to fill
  void fillGhosts(ParticleData<P>& a_particlesWithGhosts,
		  const int a_numGhost) const;

  bool isDefined() const;

private:

  void allocateVector();

  IntVect enforcePeriodic(const IntVect& a_index) const;

  List<P> m_outcast;
  ProblemDomain m_physDomain;
  RealVect m_meshSpacing;
  RealVect m_origin;
  int m_fixedBoxSize;
  ListBoxFactory<P> m_factory;
  bool m_isDefined;
};

struct CompIntVect
{
  bool operator() (const IntVect& a, const IntVect& b) const
  {
    D_TERM6(
            if      (a[0]>b[0]) return true;
            else if (a[0]<b[0]) return false;
            ,
            if      (a[1]>b[1]) return true;
            else if (a[1]<b[1]) return false;
            ,
            if      (a[2]>b[2]) return true;
            else if (a[2]<b[2]) return false;
            ,
            if      (a[3]>b[3]) return true;
            else if (a[3]<b[3]) return false;
            ,
            if      (a[4]>b[4]) return true;
            else if (a[4]<b[4]) return false;
            ,
            if      (a[5]>b[5]) return true;
            else if (a[5]<b[5]) return false;
            );
    return false;
  }
};

/// Routine for transfering particles between AMR levels. Arguments:
/// a_partValid: the destination to which particles will be sent if 
///              they meet the requirements.
/// a_PD:        the source from which the particles will be transfered.
/// a_mask:      this defines the valid region. If a_flip is false, a 
///              particle will be transfered from the source to the 
///              dest if it lies in a region of the domain where the mask
///              is true. 
/// a_dx:        particle position will be converted to integer indices
///              by binning them by this cell spacing. 
/// a_refRatio:  a_mask will be copied onto a boxlayout that has been 
///              refined or coarsened by this argument (refined if 
///              positive, coarsened if negative.) This is to allow
///              comparison with mask arrays defined on arbitrary 
///              refinement levels
/// a_flip       if true, the sign of mask will be treated as reversed
///              (i.e., this becomes a collectInvalidParticles function).
template <class P>
void collectValidParticles(List<P>&                         a_partValid,
                           ParticleData<P>&                 a_PD,
                           const LevelData<BaseFab<bool> >* a_mask,
                           const RealVect                   a_meshSpacing,
                           const int                        a_refRatio,
                           const bool                       a_flip=false,
                           const RealVect                   a_origin=RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));

/// overloaded version of the above
template <class P>
void collectValidParticles(ParticleData<P>&                 a_PDValid,
                           ParticleData<P>&                 a_PD,
                           const LevelData<BaseFab<bool> >* a_mask,
                           const RealVect                   a_meshSpacing,
                           const int                        a_refRatio,
                           const bool                       a_flip=false,
                           const RealVect                   a_origin=RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));


#include "NamespaceFooter.H"

// Implementation
#include "ParticleDataI.H"

#endif
