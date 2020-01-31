#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFREMAPPER_H_
#define _MFREMAPPER_H_

#include "MFCellFAB.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "MFIndexSpace.H"
#include "NamespaceHeader.H"

///
/**
   This class is used by the time-dependent aspects of MFChombo.
 This class handles the process of mapping the state data within
 one LevelData<MFCellFAB> based on one MFIndexSpace, to a whole new
 LevelData<MFCellFAB> based on a different MFIndexSpace.

 */
class MFRemapper
{
public:

  MFRemapper();

  virtual ~MFRemapper();

  ///
  /**
     Remap operation that takes an existing state data the source MFIndexSpace
     and maps it too the same DisjointBoxLayout but with a new MFIndexSpace.

     The new MFIndexSpace is usually the result an interface movement operation, followed
     by some form of LevelSet computation.  At the end of this, a new MFIndexSpace is
     created.  The current state data needs to be mapped into a new LevelData<MFCellFAB>
     based on this new MFIndexSpace. The DisjointBoxLayouts from a_source and a_dest
     are assumed to be identical.  Thus, no coarse-fine AMR operations are performed.
  */
  void remap(const MFIndexSpace& a_sourceMF,
             const LevelData<MFCellFAB>& a_source,
             const MFIndexSpace& a_destMF,
             LevelData<MFCellFAB>& a_dest) const;

  ///
  /**  Regridding version of remapping operation.  Same MFIndexSpace, same
       state data, different grid configurations.  This routine is usually called
       after a MeshRefine operation.
  */
  void remap(const MFIndexSpace& a_MF,
             const ProblemDomain& a_domainCoar,
             const ProblemDomain& a_domainFine,
             const LevelData<MFCellFAB>& a_source,
             const LevelData<MFCellFAB>& a_coarse,
             const int& nref,
             const int& nghost,
             LevelData<MFCellFAB>&  a_dest);
};

#include "NamespaceFooter.H"
#endif
