#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Include files:

#include "MultiBlockMeshRefine.H"
#include "MayDay.H"
#include "parstream.H"
#include <list>
#include "NamespaceHeader.H"


    /// Full constructor -- leaves object in usable state  
MultiBlockMeshRefine::MultiBlockMeshRefine(
                       /// Level 0 domain.  used to index into factory object for levels
                       const ProblemDomain& a_baseDomain,
                       /// level 0 grid spacing.  needed for mapping information
                       RealVect a_dxCoarse,
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
                       const int a_maxSize)
  :BRMeshRefine(a_baseDomain, a_refRatios, a_fillRatio, a_blockFactor, a_bufferSize, a_maxSize), 
   m_factory(a_factory), m_dxCoarse(a_dxCoarse) { }


MultiBlockMeshRefine::~MultiBlockMeshRefine()
{
  std::map<ProblemDomain, MultiBlockCoordSys*>::iterator p = m_coords.begin();
  while (p!= m_coords.end())
    {
      delete p->second;
      p++;
    }
}
bool MultiBlockMeshRefine::properlyNested(const Box& a_box,
                                const ProblemDomain& a_domain,
                                const IntVectSet& a_pnd,
                                int   a_totalBufferSize) const
{
  if (m_PNDMode == 0)
    {
      MayDay::Error("MultiBlockMeshRefine only implemented for PNDMode != 0");
    }
  else
    {
      Box growBox(a_box);
      growBox.grow(a_totalBufferSize * m_refineDirs);
      const MultiBlockCoordSys* coords = getCoords(a_domain);
      int blockID = coords->whichBlock(a_box);
      if (blockID == -1)
        {
          return false; // not contained within a single Block
        }
      Box block = coords->mappingBlocks()[blockID];
      if (!a_pnd.contains(growBox & block))
        {
          return false; // first simply check if grown box is nice already without block games
        }
 
 
      const Tuple<BlockBoundary, 2*SpaceDim>& boundaries = coords->boundaries()[blockID];
      for (int dd = 0; dd<CH_SPACEDIM; ++dd)
        {
          if (a_box.smallEnd(dd) <= block.smallEnd(dd)+a_totalBufferSize)
            {
              const BlockBoundary& bb = boundaries[dd];
              if (!bb.isDomainBoundary())
                {
                  const Box& otherBlock = coords->mappingBlocks()[bb.neighbor()];
                  IndicesTransformation tr = bb.getTransformation();
                  Box tbox = tr.transformFwd(growBox);
                  tbox &= otherBlock;
                  if (!a_pnd.contains(tbox))
                    {
                      return false;
                    }
                } 
            }
          if (a_box.bigEnd(dd) >= block.bigEnd(dd)-a_totalBufferSize)
            {
              const BlockBoundary& bb = boundaries[dd+CH_SPACEDIM];
              if (!bb.isDomainBoundary())
                {
                  const Box& otherBlock = coords->mappingBlocks()[bb.neighbor()];
                  IndicesTransformation tr = bb.getTransformation();
                  Box tbox = tr.transformFwd(growBox);
                  tbox &= otherBlock;
                  if (!a_pnd.contains(tbox)) 
                    {
                      return false;
                    }
                } 
            }
        }
    }
  
  return true;
}

 /// Mapped Boundaries
void  MultiBlockMeshRefine::buildSupport(const ProblemDomain& lvldomain, 
                                         Vector<Box>& lvlboxes, 
                                         IntVectSet& modifiedTags)
{
  const MultiBlockCoordSys* coords = getCoords(lvldomain);
  std::list<Box> lastBoxes;
  for (int i=0; i<lvlboxes.size(); i++) lastBoxes.push_back(lvlboxes[i]);
  std::list<Box> nextBoxes;
  std::list<Box>::iterator it;
  for (int d = 0; d <= CH_SPACEDIM; d++)
    {
      for (it=lastBoxes.begin(); it!= lastBoxes.end(); ++it)
        {
          const Box& b = *it;
          int blockID = coords->whichBlock(b);
          if (blockID != -1)
            {
              // easy case, whole box is in one block, jump to next box
              modifiedTags |= b;
              continue;
            }
          blockID = coords->whichBlockOverlap(b);
          if (blockID == -1) 
            {
              //MayDay::Error("Box + buffer reaches into more than one block");
              continue;
            }
          Box block = coords->mappingBlocks()[blockID];
          Box clippedBox = b & block;
          modifiedTags |= clippedBox;
          const Tuple<BlockBoundary, 2*SpaceDim>& boundaries = coords->boundaries()[blockID];
          for (int dd = 0; dd<CH_SPACEDIM; ++dd)
            {
              if (clippedBox.smallEnd(dd) == block.smallEnd(dd))
                {
                  const BlockBoundary& bb = boundaries[dd];
                  if (!bb.isDomainBoundary())
                  {
                    IndicesTransformation tr = bb.getTransformation();
                    Box tbox = tr.transformFwd(b);
                    const Box& otherBlock = coords->mappingBlocks()[bb.neighbor()];
                    nextBoxes.push_back(tbox & otherBlock);
                  } 
                }
              if (clippedBox.bigEnd(dd) == block.bigEnd(dd))
                {
                  const BlockBoundary& bb = boundaries[dd+CH_SPACEDIM];
                  if (!bb.isDomainBoundary())
                  {
                    IndicesTransformation tr = bb.getTransformation();
                    Box tbox = tr.transformFwd(b);
                    const Box& otherBlock = coords->mappingBlocks()[bb.neighbor()];
                    nextBoxes.push_back(tbox & otherBlock);
                  } 
                }
            }
        }
      lastBoxes = nextBoxes;
      nextBoxes.clear();
    }
}


void MultiBlockMeshRefine::clipBox(Box& a_box, const ProblemDomain& a_domain) const
{
  //MultiBlockCoordSys* coords = getCoords(lvldomain);
  // cubed sphere never clips a box by the domain.  all is permitted!
  // so, do nothing.
}

const MultiBlockCoordSys* MultiBlockMeshRefine::getCoords(const ProblemDomain& a_lvldomain) const
{
  std::map<ProblemDomain, MultiBlockCoordSys*>::const_iterator p = m_coords.find(a_lvldomain);
  if (p==m_coords.end())
    {
      MultiBlockCoordSys* newCoords = m_factory->getCoordSys(a_lvldomain, m_dxCoarse);
      m_coords[a_lvldomain] = newCoords;
      return newCoords;
    }
  return p->second;
}

#include "NamespaceFooter.H"
