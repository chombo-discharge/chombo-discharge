#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKMESHREFINE_H_
#define _MULTIBLOCKMESHREFINE_H_

#include "BRMeshRefine.H"
#include "MultiBlockCoordSys.H"
#include "RefCountedPtr.H"
#include <map>

#include "NamespaceHeader.H"

/**
  A BRMeshRefine algorithm for MultiBlock domains.   Only needs to cleverly intercept the virtual
  makePND operation.

  bvs: March 2012
*/
class MultiBlockMeshRefine : public BRMeshRefine
{
public:

  /// Destructor
  virtual ~MultiBlockMeshRefine();


    /// Full constructor -- leaves object in usable state  
  MultiBlockMeshRefine(
                       /// Level 0 domain.  used to index into factory object for levels
                       const ProblemDomain& a_baseDomain,
                       /// level 0 grid spacing.  needed for mapping information
                       RealVect dxCoarse,
                       const RefCountedPtr<MultiBlockCoordSysFactory>& a_factory,
                       /// Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                       const Vector<int>& a_refRatios,
                       /// Measure of how efficiently tagged cells will be covered
                       const Real a_fillRatio,
                       /// Amount by which grids are guaranteed to be coarsenable
                       const int a_blockFactor,
                       /// Proper nesting buffer amount
                       const int a_bufferSize,
                       /// Maximum grid length in any direction -- 0 means no limit.
                       const int a_maxSize);



protected:


  /// This one is to let Mapped Multiblock grids respect the PND across mapped boundaries
  virtual bool properlyNested(const Box& a_box,
                              const ProblemDomain& a_domain,
                              const IntVectSet& a_pnd,
                              int   a_totalBufferSize) const;

  ///These two are to get the tagging to respect Mapped Boundaries
  virtual void buildSupport(const ProblemDomain& lvldomain, Vector<Box>& lvlboxes, IntVectSet& modifiedTags);
  virtual void clipBox(Box& a_box, const ProblemDomain& a_domain) const ;

  const MultiBlockCoordSys* getCoords(const ProblemDomain& a_lvldomain) const;

  /// Implementation/helper class 
  /// Default constructor -- leaves object in an unusable state.  not callable by outside classes
  /// strong construction only for now.
  MultiBlockMeshRefine() {}
  RefCountedPtr<MultiBlockCoordSysFactory> m_factory;
  mutable std::map<ProblemDomain, MultiBlockCoordSys*> m_coords;
  RealVect m_dxCoarse;

};

#include "NamespaceFooter.H"
#endif
