#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKFLUXREGISTER_H_
#define _MULTIBLOCKFLUXREGISTER_H_

#include "LevelFluxRegister.H"
#include "MultiBlockCoordSys.H"
#include <map>
#include "RemoteCopier.H"
#include "NamespaceHeader.H"

//! \class MultiBlockFluxRegister
//! This subclass of LevelFluxRegister handles multi-block AMR refluxing.
class MultiBlockFluxRegister: public LevelFluxRegister
{
  public:

  //! Creates a flux register that handles multi-block refluxing.
  //! This capability requires a multi-block coordinate system, which
  //! handles block connectivity as well as coordinate mappings.
  explicit MultiBlockFluxRegister(MultiBlockCoordSys* a_coordSysCoarse,
                                  MultiBlockCoordSys* a_coordSysFine);

  //! Destructor.
  ~MultiBlockFluxRegister();

  ///
  /**
     dProblem is at the fine resolution.
     scaleFluxes means the same thing as in levelfluxregister
   */
  void define(const DisjointBoxLayout& a_dbl,
              const DisjointBoxLayout& a_dblCoarse,
              const ProblemDomain&     a_dProblem,
              int                      a_nRefine,
              int                      a_nComp,
              bool                     a_scaleFineFluxes,
              int                      a_unrefinedDirection=-1);

  ///
  /**
     the vector interval is the interval over which the fancy vector transformations
     are done.
   */
  void reflux(LevelData<FArrayBox>& a_uCoarse,
              const Interval&       a_coarseVectorIntv,
              Real                  a_scale);

  ///
  /**
     the coarse interval is the interval over which refluxing is done.
     the vector interval is the interval over which the fancy vector transformations
     are done.
   */
  void reflux(LevelData<FArrayBox>& a_uCoarse,
              const Interval&       a_coarse_interval,
              const Interval&       a_coarseVectorIntv,
              const Interval&       a_flux_interval,
              Real                  a_scale);

  ///
  /**
     Increment fine flux.  For scalars, no change from LevelFluxRegister.
     For vector fluxes, convert to frame of opposite block.
     This is for face-CENTERED fluxes.
   */
  void incrementFineCentered(const FArrayBox& a_fineFlux,
                             Real a_scale,
                             const DataIndex& a_fineDataIndex,
                             const Interval& a_vectorFluxIntv,
                             const Interval& a_srcInterval,
                             const Interval& a_dstInterval,
                             int a_dir,
                             Side::LoHiSide a_sd);

  /// Call with a_sd == Side::Lo and a_sd == Side::Hi.
  void incrementFineCentered(const FArrayBox& a_fineFlux,
                             Real a_scale,
                             const DataIndex& a_fineDataIndex,
                             const Interval& a_vectorFluxIntv,
                             const Interval& a_srcInterval,
                             const Interval& a_dstInterval,
                             int a_dir);

  ///
  /**
     Increment fine flux.  For scalars, no change from LevelFluxRegister.
     For vector fluxes, convert to frame of opposite block.
     This is for face-AVERAGED fluxes.
     The argument a_fineFlux2 is a 2nd-order approximation that
     includes a ghost layer.
   */
  void incrementFineAveraged(const FArrayBox& a_fineFlux,
                             const FArrayBox& a_fineFlux2,
                             Real a_scale,
                             const DataIndex& a_fineDataIndex,
                             const Interval& a_vectorFluxIntv,
                             const Interval& a_srcInterval,
                             const Interval& a_dstInterval,
                             int a_dir,
                             Side::LoHiSide a_sd);

  /// Call with a_sd == Side::Lo and a_sd == Side::Hi.
  void incrementFineAveraged(const FArrayBox& a_fineFlux,
                             const FArrayBox& a_fineFlux2,
                             Real a_scale,
                             const DataIndex& a_fineDataIndex,
                             const Interval& a_vectorFluxIntv,
                             const Interval& a_srcInterval,
                             const Interval& a_dstInterval,
                             int a_dir);

  private:

  //! Multi-block coordinate system at coarser level.
  MultiBlockCoordSys* m_coordSysCoarse;

  //! Multi-block coordinate system at finer level.
  MultiBlockCoordSys* m_coordSysFine;

  //! Copier for moving values from fine fluxes on other blocks to
  //! blocks where the corresponding coarse fluxes are defined.
  //! One copier per dimension.
  RemoteCopier m_remoteCopiers[SpaceDim];

  //! This stores (coarsened) fine fluxes on neighboring blocks.
  //! We can ensure that the "destination" box layout is disjoint by
  //! storing the different flux components separately.
  LevelData<FArrayBox> m_multiblockFineFluxes[SpaceDim];

  //! This stores (coarsened) fine fluxes that have been copied to
  //! the blocks containing the coarse fluxes they will replace.
  //! Once again, components of fluxes in each direction are stored
  //! separately so that the boxes can be disjoint.
  LevelData<FArrayBox> m_neighborFineFluxes[SpaceDim];

  // Index of the block containing each coarsened fine grid.
  LayoutData<int> m_blockNum;

  //! No default constructor.
  MultiBlockFluxRegister()
  { };
};
#include "NamespaceFooter.H"
#endif
