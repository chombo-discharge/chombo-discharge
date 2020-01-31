#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewMultiBlockFluxRegister.H"
#include "BlockRegister.H"
#include "LevelFluxRegisterF_F.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

void setZero(const Box& b, int comp, FluxBox& T) { T.setVal(0);}

NewMultiBlockFluxRegister::
~NewMultiBlockFluxRegister()
{

}

void 
NewMultiBlockFluxRegister::
define( RefCountedPtr<MultiBlockCoordSys>  a_domCoar,
        const DisjointBoxLayout&           a_dblFine,
        const DisjointBoxLayout&           a_dblCoar,
        int                                a_nRefine,
        int                                a_nComp,
        bool                               a_scaleFineFluxes,
        int                                a_unrefinedDirection)
 
{
  CH_TIME("NewMultiBlockFluxRegister::define");
  LevelFluxRegister::define(a_dblFine,a_dblCoar,a_dblFine.physDomain(),
                            a_nRefine,a_nComp,a_scaleFineFluxes, a_unrefinedDirection);
  m_coarRegisters.resize(a_nComp);
  m_cofiRegisters.resize(a_nComp);
  for (int k =0; k < a_nComp; k++)
    {      
      RefCountedPtr<BlockRegister >
        nextBR(new BlockRegister(a_domCoar,a_dblCoar));
      m_coarRegisters[k] = nextBR;
    }
  for (int k =0; k < a_nComp; k++)
    {      
      RefCountedPtr<BlockRegister >
        nextBR(new BlockRegister(a_domCoar,a_dblCoar));
      m_cofiRegisters[k] = nextBR;
    }
  
  DisjointBoxLayout dblCoFi;
  coarsen(dblCoFi,a_dblFine,m_nRefine);
  m_fluxCoFi.define(dblCoFi,a_nComp);
  
  m_fluxCoFi.apply(setZero);
  m_nComp = a_nComp;
  defineBlockBoundaries( a_domCoar, a_dblFine);
  m_blockCFLocations.resize(0);
  defineBlockCFLocations(a_domCoar, a_dblFine,  a_dblCoar);
  
}

int
NewMultiBlockFluxRegister::
index(int dir, Side::LoHiSide side)
{
  CH_assert(dir >= 0);
  CH_assert(dir < SpaceDim);
  CH_assert((side == Side::Lo) || (side == Side::Hi));

  int ioffset;
  if (side == Side::Lo)
    ioffset = 0;
  else
    ioffset = 1;

  return ioffset*SpaceDim+dir;
}

void 
NewMultiBlockFluxRegister::
defineBlockCFLocations(RefCountedPtr<MultiBlockCoordSys> a_domCoar,
                       const DisjointBoxLayout&          a_dblFine,
                       const DisjointBoxLayout&          a_dblCoar)
{
  CH_TIME("NewMultiBlockFluxRegister::defineBlockCFLocations");
  std::list<Box> blockBoxesCoFi;
  //get list of boxes which map to fine boxes at  block boundaries.
  //it is OK for there to be a fine grid over it because we are planning to 
  //average down which will overwrite data under finer grids.   
  getCoFiBoxBlockList(blockBoxesCoFi, a_domCoar, a_dblFine);
  m_blockCFLocations.resize(2*SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int iindex = index(idir, sit());
          
          m_blockCFLocations[iindex] = RefCountedPtr<BlockRegister>(new BlockRegister());
          m_blockCFLocations[iindex]->define(a_domCoar,a_dblCoar);
          m_blockCFLocations[iindex]->setToZero(1);

          for (DataIterator dit = a_dblCoar.dataIterator(); dit.ok(); ++dit)
            {
              const Box& coarBox = a_dblCoar[dit()];
              for (std::list<Box>::iterator lbit=blockBoxesCoFi.begin(); lbit != blockBoxesCoFi.end(); ++lbit) 
                {
                  const Box& cofiBoxTrans = *lbit;
                  Box intersectBox = coarBox;
                  intersectBox.shift(idir, sign(sit()));
                  intersectBox &= cofiBoxTrans;
                  if (!intersectBox.isEmpty())
                    {
                      //we have a coarse fine interface.   set the block register on the other side 
                      //of the interface to 1
                      Box faceBox = surroundingNodes(intersectBox, idir);
                      FArrayBox oneFlux(faceBox, 1);
                      oneFlux.setVal(1.0);
                      //high side of the box gets 1 on the high side of the register
                      m_blockCFLocations[iindex]->storeFlux(oneFlux, dit(), idir, sit());
                    }
                }
            }

          m_blockCFLocations[iindex]->close();
        }
    }
}
/////////////////
//get fine boxes at block boundaries
void
NewMultiBlockFluxRegister::
getCoFiBoxBlockList(std::list<Box>                   & a_blockBoxesCoFi, 
                    RefCountedPtr<MultiBlockCoordSys>& a_domCoar,
                    const DisjointBoxLayout          & a_dblFine)
{
  CH_TIME("NewMultiBlockFluxRegister::getCoFiBoxBlockList");
  //loop through all directions and sides so we can get all the boxes 
  //at all block boundaries in one list
  a_blockBoxesCoFi.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          for (LayoutIterator lit = a_dblFine.layoutIterator(); lit.ok(); ++lit)
            {
              const Box& fineBox = a_dblFine[lit()];
              Box boxCoFi = coarsen(fineBox, m_nRefine);
              Box shiftedCurBox = boxCoFi;
              shiftedCurBox.shift(idir, sign(sit()));
              int blockID = a_domCoar->whichBlock(shiftedCurBox);
              if (blockID != -1)
                {
                  // easy case, whole shifted box is in one block, jump to next box
                  continue;
                }
              else
                {
                  blockID = a_domCoar->whichBlockOverlap(shiftedCurBox);
                  if (blockID == -1) 
                    {
                      MayDay::Error("Box + buffer reaches into more than one block");
                    }
                  const Tuple<BlockBoundary, 2*SpaceDim>& boundaries = a_domCoar->boundaries()[blockID];
                  
                  int iindex = index(idir, sit());
                  const BlockBoundary& bb = boundaries[iindex];
                  if (!bb.isDomainBoundary())
                    {
                      IndicesTransformation tr = bb.getTransformation();
                      Box boxTrans = tr.transformFwd(boxCoFi);
                      a_blockBoxesCoFi.push_back(boxTrans);
                    }
                }
            }
        }
    }
}
void 
NewMultiBlockFluxRegister::
defineBlockBoundaries( RefCountedPtr<MultiBlockCoordSys> a_domCoar,
                       const DisjointBoxLayout&          a_dblFine)
{
  CH_TIME("NewMultiBlockFluxRegister::defineBlockBoundaries");
  m_fineBlockBoundaries.define(a_dblFine);
  for (DataIterator dit =a_dblFine.dataIterator();dit.ok();++dit)
    {
      Vector<Box>&  edge = m_fineBlockBoundaries[dit];
      edge.resize(2*SpaceDim); 

      Box curBox = a_dblFine[dit()];
      Box coarsenedCurBox = coarsen(curBox,m_nRefine);
      for (int k = 0; k < a_domCoar->numBlocks();k++)
        {
          const Box& block = (a_domCoar->mappingBlocks())[k];
          if (block.contains(coarsenedCurBox))
            {
              for (SideIterator sd;sd.ok();sd.next())
                {
                  for (int dir = 0 ; dir < CH_SPACEDIM; dir++)
                    {
                      Box shiftedCurBox = coarsenedCurBox;
                      shiftedCurBox.shift(dir,sign(sd()));
                      int iindex = index(dir, sd());
                      if (!block.contains(shiftedCurBox))
                        {
                          edge[iindex] = 
                            bdryBox(curBox,dir,sd(),1);
                        }
                    }
                }
            }
        }
    } 
}

void 
NewMultiBlockFluxRegister::
reflux(LevelData<FArrayBox>& a_uCoarse,
       const Interval&       a_interval,
       Real                  a_scale)
{
  this->reflux(a_uCoarse, a_interval, a_interval, a_scale);
}

void 
NewMultiBlockFluxRegister::
reflux(LevelData<FArrayBox>& a_uCoarse,
       const Interval&       a_coarseVectorIntv,
       const Interval&       a_fluxInterval,
       Real                  a_scale)
{
  CH_TIME("NewMultiBlockFluxRegister::reflux");
  LevelFluxRegister::reflux(a_uCoarse,
                            a_coarseVectorIntv,
                            a_fluxInterval,
                            a_scale);
  
  /**/  
  DisjointBoxLayout dblCoarse = a_uCoarse.getBoxes();
  LevelData<FluxBox > fluxFromFine(dblCoarse, a_uCoarse.nComp());
  //coarsenedFine.setToZero
  for (DataIterator dit = dblCoarse.dataIterator(); dit.ok(); ++dit)
    {
      fluxFromFine[dit].setVal(0);
    }
  m_fluxCoFi.copyTo(fluxFromFine);
  for (int ivar=0; ivar<m_coarRegisters.size(); ivar++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              for (DataIterator dit = dblCoarse.dataIterator(); dit.ok(); ++dit)
                {
                  if (m_coarRegisters[ivar]->hasInterface(dit(),idir,sit()))
                    {
                      FluxBox singleVarFlux(dblCoarse[dit()], 1);
                      singleVarFlux.copy(fluxFromFine[dit()], ivar, 0, 1);
                      m_cofiRegisters[ivar]->storeFlux(singleVarFlux[idir], dit(), idir, sit());
                    }
                }
            }
        }
      m_cofiRegisters[ivar]->close();
      m_coarRegisters[ivar]->close();
    }
      
  Real tol = 1.0e-3;
  for (DataIterator dit = dblCoarse.dataIterator(); dit.ok(); ++dit)
    {
      for (SideIterator sit;sit.ok();sit.next())
        {
          for (int idir = 0;idir < CH_SPACEDIM;idir++)
            {
              if (m_coarRegisters[0]->hasInterface(dit(),idir,sit()))
                {
                  
                  const Box& bx =dblCoarse[dit];
                  Box bxBndryFace = bdryBox(bx, idir, sit(), 1);
                  FArrayBox coarFlux(bxBndryFace,a_uCoarse.nComp());
                  FArrayBox cofiFlux(bxBndryFace,a_uCoarse.nComp());
                  coarFlux.setVal(0.);
                  cofiFlux.setVal(0.);
                  int kstart = a_coarseVectorIntv.begin();
                  int kend   = a_coarseVectorIntv.end();
                  //int ksize = kend - kstart + 1;
                  for (int k2 = kstart;k2 <= kend;k2++)
                    {
                      FArrayBox coarFluxAlias(Interval(k2,k2), coarFlux);
                      FArrayBox cofiFluxAlias(Interval(k2,k2), cofiFlux);
                      m_coarRegisters[k2]->getFlux(coarFluxAlias, dit(),idir, sit(), sit());
                      m_cofiRegisters[k2]->getFlux(cofiFluxAlias, dit(),idir, sit(), sit());
                    }
                  FArrayBox flagsThisSide(bxBndryFace, 1);
                  FArrayBox flagsFlipSide(bxBndryFace, 1);
                  int iindex = index(idir, sit());
                  const BlockRegister& flagblock = *m_blockCFLocations[iindex];
                  //the opposite side of the interface on the high side of the 
                  //box is the high side of the interface
                  flagblock.getFlux(flagsThisSide, dit(), idir, sit(), flip(sit()));
                  flagblock.getFlux(flagsFlipSide, dit(), idir, sit(),      sit());

                  //shift flux to cells
                  coarFlux.shiftHalf(     idir,-sign(sit()));
                  cofiFlux.shiftHalf(     idir,-sign(sit()));
                  flagsThisSide.shiftHalf(idir,-sign(sit()));
                  flagsFlipSide.shiftHalf(idir,-sign(sit()));
                  bxBndryFace.shiftHalf(  idir,-sign(sit()));
                  for (BoxIterator bit(bxBndryFace); bit.ok(); ++bit)
                    {
                      Real coarSide = flagsFlipSide(bit(),0);  
                      Real fineSide = flagsThisSide(bit(),0);  
                      bool reallyCoar = Abs(coarSide    ) < tol;
                      bool reallyFine = Abs(fineSide-1.0) < tol;
                      //only change the solution where there really is a C/F interface
                      if (reallyCoar && reallyFine)
                        {
                          for (int ivar = kstart; ivar <= kend; ivar++)
                            {
                              Real coar =  coarFlux(bit(), ivar);
                              Real cofi =  cofiFlux(bit(), ivar);
                              Real diff =  cofi - coar;
                              diff *= -a_scale*sign(sit());
                              a_uCoarse[dit](bit(), ivar) += diff;
                            }
                        }
                    }
                  //shift back
                  coarFlux.shiftHalf(     idir,sign(sit()));
                  cofiFlux.shiftHalf(     idir,sign(sit()));
                  flagsThisSide.shiftHalf(idir,sign(sit()));
                  flagsFlipSide.shiftHalf(idir,sign(sit()));
                  bxBndryFace.shiftHalf(  idir,sign(sit()));
                }
            }
        }
    }
}

void
NewMultiBlockFluxRegister::
reflux(LevelData<FArrayBox>& a_uCoarse,
       Real a_scale)
{
  Interval interval(0, a_uCoarse.nComp()-1);
  reflux(a_uCoarse, interval, interval, a_scale);
}
void 
NewMultiBlockFluxRegister::
reflux(LevelData<FArrayBox>& a_uCoarse,
       const Interval&       a_coarse_interval,
       const Interval&       a_coarseVectorIntv,
       const Interval&       a_flux_interval,
       Real                  a_scale)
{
  MayDay::Error("We didn't write this form of the flux register yet");
}

void 
NewMultiBlockFluxRegister::
setToZero()
{
  LevelFluxRegister::setToZero();
  for (int k = 0; k < m_coarRegisters.size();k++)
    {
      m_coarRegisters[k]->setToZero(1);
      m_cofiRegisters[k]->setToZero(1);
    }
}

void  
NewMultiBlockFluxRegister:: 
incrementCoarse(
                const FArrayBox& a_coarseFlux,
                Real a_scale,
                const DataIndex& a_coarseDataIndex,
                const Interval& a_srcInterval,
                const Interval& a_dstInterval,
                int a_dir,
                Side::LoHiSide a_sd)
{
  CH_TIME("NewMultiBlockFluxRegister::incrementCoarse");
  CH_assert(a_srcInterval.size() == a_dstInterval.size());
  LevelFluxRegister::incrementCoarse(
                                     a_coarseFlux,
                                     a_scale,
                                     a_coarseDataIndex,
                                     a_srcInterval,
                                     a_dstInterval,
                                     a_dir,
                                     a_sd);

  for (SideIterator sd ; sd.ok() ; sd.next())
    {
      int dir = a_dir;
      if (m_coarRegisters[0]->hasInterface(a_coarseDataIndex,dir,sd()))
        {
          FArrayBox fluxComp(a_coarseFlux.box(),1);
          int kSrc=a_srcInterval.begin();
          int kDest=a_dstInterval.begin();
          for (; kSrc <= a_srcInterval.end(); kSrc++, kDest++)  
            { 
              fluxComp.copy(a_coarseFlux,kSrc,0,1);
              m_coarRegisters[kDest]->storeFlux(fluxComp,a_coarseDataIndex,
                                                dir,flip(sd()));
            }
        }
    }
}

void 
NewMultiBlockFluxRegister::
incrementFine(
              const FArrayBox& a_fineFlux,
              Real a_scale,
              const DataIndex& a_fineDataIndex,
              const Interval& a_srcInterval,
              const Interval& a_dstInterval,
              int a_dir,
              Side::LoHiSide a_sd)
{
  CH_TIME("NewMultiBlockFluxRegister::incrementFine");
  LevelFluxRegister::incrementFine(
                                   a_fineFlux,
                                   a_scale,
                                   a_fineDataIndex,
                                   a_srcInterval,
                                   a_dstInterval,
                                   a_dir,
                                   a_sd);
  int iindex = index(a_dir, a_sd);
  Box bx = m_fineBlockBoundaries[a_fineDataIndex][iindex];
  if (!bx.isEmpty())
    {
      FArrayBox& coarsenedFlux = m_fluxCoFi[a_fineDataIndex][a_dir];
      const IntVect& iv=bx.smallEnd();
      IntVect civ = coarsen(iv, m_nRefine);
      int srcComp = a_srcInterval.begin();
      int destComp = a_dstInterval.begin();
      int ncomp = a_srcInterval.size();
      Real denom = 1.0;

      if (m_scaleFineFluxes)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir != a_dir)
                {
                  denom *= m_nRefine[idir];
                }
            }
        }

      Real scale = a_scale/denom;

      FORT_INCREMENTFINE(CHF_CONST_FRA_SHIFT(a_fineFlux, iv),
                         CHF_FRA_SHIFT(coarsenedFlux, civ),
                         CHF_BOX_SHIFT(bx, iv),
                         CHF_CONST_INTVECT(m_nRefine),
                         CHF_CONST_REAL(scale),
                         CHF_CONST_INT(srcComp),
                         CHF_CONST_INT(destComp),
                         CHF_CONST_INT(ncomp));
      
    }
}

#include "NamespaceFooter.H"
