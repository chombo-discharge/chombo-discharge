#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEWMULTIBLOCKFLUXREGISTER_H_
#define _NEWMULTIBLOCKFLUXREGISTER_H_

#include "LevelFluxRegister.H"
#include "BlockRegister.H"
#include "NamespaceHeader.H"
//NewMultiBlockFluxRegister-Derived from LevelFluxRegister to support multiblock grids.
/**
   A LevelFluxRegister  handles all the data choreography
   necessary to create a levels worth of flux registers.
*/
///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers. a_domain is is coarse-level
     domain.

     If a_scaleFineFluxes is true (default then all fine fluxes are
     scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
  */
class NewMultiBlockFluxRegister : public LevelFluxRegister
{
  public:
  ///
  /**
     domainCoar --coordinate system at coarser level
   */
  NewMultiBlockFluxRegister(
                         RefCountedPtr<MultiBlockCoordSys> a_domCoar,
                         const DisjointBoxLayout& a_dblFine,
                         const DisjointBoxLayout& a_dblCoar,
                         int                      a_nRefine,
                         int                      a_nComp,
                         bool                     a_scaleFineFluxes=false,
                         int                      a_unrefinedDirection=-1) 
  {
    define(a_domCoar,
           a_dblFine,
           a_dblCoar,
           a_nRefine,
           a_nComp,
           a_scaleFineFluxes,
           a_unrefinedDirection);
  }

  static int index(int a_dir, Side::LoHiSide a_side);

  NewMultiBlockFluxRegister(const NewMultiBlockFluxRegister&);

  NewMultiBlockFluxRegister()
  { };

  ~NewMultiBlockFluxRegister();
  
  ///
  /**
     domainCoar --coordinate system at coarser level
   */
  virtual void define(RefCountedPtr<MultiBlockCoordSys> a_domainCoar,
                      const DisjointBoxLayout& a_dbl,
                      const DisjointBoxLayout& a_dblCoarse,
                      int                      a_nRefine,
                      int                      a_nComp,
                      bool                     a_scaleFineFluxes=false,
                      int                      a_unrefinedDirection=-1);

  virtual void reflux(
              LevelData<FArrayBox>& a_uCoarse,
              Real a_scale);

  virtual void reflux(LevelData<FArrayBox>& a_uCoarse,
              const Interval&       a_coarseVectorIntv,
              Real                  a_scale);

  virtual void reflux(LevelData<FArrayBox>& a_uCoarse,
                      const Interval&       a_coarseVectorIntv,
                      const Interval&       a_fluxInterval,
                      Real                  a_scale);

  virtual void reflux(LevelData<FArrayBox>& a_uCoarse,
              const Interval&       a_coarse_interval,
              const Interval&       a_coarseVectorIntv,
              const Interval&       a_flux_interval,
              Real                  a_scale);

  virtual void setToZero();

  virtual void incrementFine(
                             const FArrayBox& a_fineFlux,
                             Real a_scale,
                             const DataIndex& a_fineDataIndex,
                             const Interval& a_srcInterval,
                             const Interval& a_dstInterval,
                             int a_dir,
                             Side::LoHiSide a_sd);

  virtual void incrementCoarse(
                               const FArrayBox& a_coarseFlux,
                               Real a_scale,
                               const DataIndex& a_coarseDataIndex,
                               const Interval& a_srcInterval,
                               const Interval& a_dstInterval,
                               int a_dir,
                               Side::LoHiSide a_sd);

private:
  Vector<RefCountedPtr<BlockRegister > > m_coarRegisters;
  Vector<RefCountedPtr<BlockRegister > > m_cofiRegisters;
  LevelData<FluxBox > m_fluxCoFi;
  Vector<RefCountedPtr<BlockRegister > > m_blockCFLocations;
  LayoutData<Vector<Box> > m_fineBlockBoundaries;
  int m_nComp;

  // internal function to make define less  monolithic
  void defineBlockCFLocations(RefCountedPtr<MultiBlockCoordSys> a_domain,
                              const DisjointBoxLayout&          a_dblFine,
                              const DisjointBoxLayout&          a_dblCoar);

  // internal function to make define less  monolithic
  void getCoFiBoxBlockList(std::list<Box>                   & a_blockBoxesCoFi, 
                           RefCountedPtr<MultiBlockCoordSys>& a_domCoar,
                           const DisjointBoxLayout          & a_dblFine);

  // internal function to make define less  monolithic
  void defineBlockBoundaries( RefCountedPtr<MultiBlockCoordSys> a_domain,
                              const DisjointBoxLayout&          a_dblFine);
  

};
#include "NamespaceFooter.H"

#endif // include guard
