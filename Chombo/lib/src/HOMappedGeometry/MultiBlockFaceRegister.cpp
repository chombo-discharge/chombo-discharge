#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockFaceRegister.H"
#include "BoxIterator.H"
#include "CH_assert.H"
#include <map>
#include "NamespaceHeader.H"


/*--------------------------------------------------------------------*/
//  Constructor
/**
*//*------------------------------------------------------------------*/

MultiBlockFaceRegister::MultiBlockFaceRegister(
  const MultiBlockCoordSys *const a_CrCoordSys)
  :
  LevelFluxRegister(),
  m_coordSys(a_CrCoordSys)
{
  CH_assert(a_CrCoordSys != NULL);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
*//*------------------------------------------------------------------*/

MultiBlockFaceRegister::~MultiBlockFaceRegister()
{
}

/*--------------------------------------------------------------------*/
//  Define everything
/**
*//*------------------------------------------------------------------*/

void
MultiBlockFaceRegister::define(const DisjointBoxLayout& a_FnGrid,
                               const DisjointBoxLayout& a_CrGrid,
                               const ProblemDomain&     a_FnProbDom,
                               int                      a_nRefine,
                               int                      a_nComp,
                               bool                     a_scaleFineFaces)
{
  defineFine(a_FnGrid,
             a_FnProbDom,
             a_nRefine*IntVect::Unit,
             a_nComp,
             a_scaleFineFaces);
  defineCoarse(a_CrGrid);
}


void
MultiBlockFaceRegister::define(const DisjointBoxLayout& a_FnGrid,
                               const DisjointBoxLayout& a_CrGrid,
                               const ProblemDomain&     a_FnProbDom,
                               IntVect                  a_ivRefine,
                               int                      a_nComp,
                               bool                     a_scaleFineFaces)
{
  defineFine(a_FnGrid,
             a_FnProbDom,
             a_ivRefine,
             a_nComp,
             a_scaleFineFaces);
  defineCoarse(a_CrGrid);
}


/*--------------------------------------------------------------------*/
//  Define the fine register
/**
*//*------------------------------------------------------------------*/

void
MultiBlockFaceRegister::defineFine(const DisjointBoxLayout& a_FnGrid,
                                   const ProblemDomain&     a_FnProbDom,
                                   int                      a_nRefine,
                                   int                      a_nComp,
                                   bool                     a_scaleFineFaces)
{
  CH_assert(m_coordSys->mappingBlocks().size() > 1);
  m_isDefined |= FluxRegDefined;

  m_nRefine = a_nRefine * IntVect::Unit;
  m_scaleFineFluxes = a_scaleFineFaces;
  m_domain = a_FnProbDom;

  DisjointBoxLayout crFnGrid;
  coarsen(crFnGrid, a_FnGrid, a_nRefine);

  m_fineFlux.define(crFnGrid, a_nComp, IntVect::Unit);
  m_isDefined |= FluxRegFineDefined;
}

void
MultiBlockFaceRegister::defineFine(const DisjointBoxLayout& a_FnGrid,
                                   const ProblemDomain&     a_FnProbDom,
                                   IntVect                  a_nRefine,
                                   int                      a_nComp,
                                   bool                     a_scaleFineFaces)
{
  CH_assert(a_nRefine > IntVect::Zero);
  CH_assert(m_coordSys->mappingBlocks().size() > 1);
  m_isDefined |= FluxRegDefined;

  m_nRefine = a_nRefine;
  m_scaleFineFluxes = a_scaleFineFaces;
  m_domain = a_FnProbDom;

  DisjointBoxLayout crFnGrid;
  coarsen(crFnGrid, a_FnGrid, a_nRefine);

  m_fineFlux.define(crFnGrid, a_nComp, IntVect::Unit);
  m_isDefined |= FluxRegFineDefined;
}

/*--------------------------------------------------------------------*/
//  Define the coarse register and copiers
/**
*//*------------------------------------------------------------------*/

void
MultiBlockFaceRegister::defineCoarse(const DisjointBoxLayout& a_CrGrid)
{
  CH_TIME("MultiBlockFaceRegister::defineCoarse");
  const Vector<Box>& blocks = m_coordSys->mappingBlocks();
  CH_assert(blocks.size() > 1);
  m_isDefined |= FluxRegDefined;

//--Setup of base class (since we don't actually call it's define function)

  const DisjointBoxLayout& crFnGrid = m_fineFlux.getBoxes();
  // Note that m_coarFlux and the reverse copier are *not* defined
  for (int i = 0; i != CH_SPACEDIM; ++i)
    {
      m_coarseLocations[i].define(a_CrGrid);
      m_coarseLocations[i+CH_SPACEDIM].define(a_CrGrid);
    }

//--Iterate over all boxes and gather information on those whose boundaries
//--coincide with block boundaries.

  const Vector<Tuple<BlockBoundary, 2*SpaceDim> >& blockBoundaries =
    m_coordSys->boundaries();

  // We construct data structures for accumulating fine faces from remote
  // blocks in each direction.
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      //**FIXME is it possible to put the data iterator on the inner loop?
      for (DataIterator ditCr(a_CrGrid.dataIterator()); ditCr.ok(); ++ditCr)
        {
          // Find the block containing this box.
          Box crBox = a_CrGrid[ditCr];
          int blockID = m_coordSys->whichBlock(crBox);
          CH_assert(blockID != -1);
          const Box& blockBox = blocks[blockID];

          // Now find out whether this face abuts the block boundary and set up
          // any block-block face transfers that are needed.
          Box crFaceBox = crBox;
          crFaceBox.grow(dir, 1);

          //--This would be a correction to coarse cells/faces at the low side of a box

          if (crFaceBox.smallEnd(dir) < blockBox.smallEnd(dir))
            {
              // Associate this box with the (low side) block boundary.
              const BlockBoundary& boundary =
                blockBoundaries[blockID][dir];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation =
                boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              Vector<Box>& coarseLocations =
                m_coarseLocations[dir+SpaceDim][ditCr];
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator litCrFn = crFnGrid.layoutIterator();
                   litCrFn.ok(); ++litCrFn)
                {
                  Box crFnBox = crFnGrid[litCrFn];

                  // Only consider boxes from the neighboring block
                  if (m_coordSys->whichBlock(crFnBox) != nbrBlockID) continue;
                  
                  crFnBox = transformation.transformBack(crFnBox);
                  if (crFnBox.intersectsNotEmpty(crFaceBox))
                    {
                      Box intersect = crFnBox & crFaceBox;
                      // Find the cells just inside the original coarse box
                      intersect = adjCellHi(intersect, dir);
                      
                      // Add this intersection to the set of coarse face locations
                      // so that we can increment coarse faces here.
                      coarseLocations.push_back(intersect);
                    }
                }
            }
          
          //--This would be a correction to coarse cells/faces at the high side of a box
          
          if (crFaceBox.bigEnd(dir) > blockBox.bigEnd(dir))
            {
              // Associate this box with the (hi side) block boundary.
              const BlockBoundary& boundary =
                blockBoundaries[blockID][dir + SpaceDim];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation =
                boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              Vector<Box>& coarseLocations =
                m_coarseLocations[dir][ditCr];
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator litCrFn = crFnGrid.layoutIterator();
                   litCrFn.ok(); ++litCrFn)
                {
                  Box crFnBox = crFnGrid[litCrFn];
                  
                  // Only consider boxes from the neighboring block
                  if (m_coordSys->whichBlock(crFnBox) != nbrBlockID) continue;
                  
                  crFnBox = transformation.transformBack(crFnBox);
                  if (crFnBox.intersectsNotEmpty(crFaceBox))
                    {
                      Box intersect = crFnBox & crFaceBox;
                      // Find the cells just inside the original coarse box
                      intersect = adjCellLo(intersect, dir);
                      
                      // Add this intersection to the set of coarse face locations
                      // so that we can increment coarse faces here.
                      coarseLocations.push_back(intersect);
                    }
                }
            }
        }

      Vector<Box> loSrcBoxes, loDstBoxes, hiSrcBoxes, hiDstBoxes;
      Vector<int> loSrcProcs, loDstProcs, hiSrcProcs, hiDstProcs;
      for (LayoutIterator litCr(a_CrGrid.layoutIterator()); litCr.ok(); ++litCr)
        {
          // Find the block containing this box.
          Box crBox = a_CrGrid[litCr];
          int dstProc = a_CrGrid.procID(litCr());
          int blockID = m_coordSys->whichBlock(crBox);
          CH_assert(blockID != -1);
          const Box& blockBox = blocks[blockID];

          // Now find out whether this face abuts the block boundary and set up
          // any block-block face transfers that are needed.
          Box crFaceBox = crBox;
          crFaceBox.grow(dir, 1);

          //--This would be a correction to coarse cells/faces at the low side of a box

          if (crFaceBox.smallEnd(dir) < blockBox.smallEnd(dir))
            {
              // Associate this box with the (low side) block boundary.
              const BlockBoundary& boundary =
                blockBoundaries[blockID][dir];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation =
                boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator litCrFn = crFnGrid.layoutIterator();
                   litCrFn.ok(); ++litCrFn)
                {
                  Box crFnBox = crFnGrid[litCrFn];
                  int srcProc = crFnGrid.procID(litCrFn());
                  
                  // Only consider boxes from the neighboring block
                  if (m_coordSys->whichBlock(crFnBox) != nbrBlockID) continue;
                  
                  crFnBox = transformation.transformBack(crFnBox);
                  if (crFnBox.intersectsNotEmpty(crFaceBox))
                    {
                      Box intersect = crFnBox & crFaceBox;
                      // Find the cells just inside the original coarse box
                      intersect = adjCellHi(intersect, dir);
                      loDstBoxes.push_back(intersect);
                      Box otherSide = transformation.transformFwd(intersect);
                      loSrcBoxes.push_back(otherSide);

                      loSrcProcs.push_back(srcProc);
                      loDstProcs.push_back(dstProc);
                    }
                }
            }
          
          //--This would be a correction to coarse cells/faces at the high side of a box
          
          if (crFaceBox.bigEnd(dir) > blockBox.bigEnd(dir))
            {
              // Associate this box with the (hi side) block boundary.
              const BlockBoundary& boundary =
                blockBoundaries[blockID][dir + SpaceDim];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation =
                boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator litCrFn = crFnGrid.layoutIterator();
                   litCrFn.ok(); ++litCrFn)
                {
                  Box crFnBox = crFnGrid[litCrFn];
                  int srcProc = crFnGrid.procID(litCrFn());
                  
                  // Only consider boxes from the neighboring block
                  if (m_coordSys->whichBlock(crFnBox) != nbrBlockID) continue;
                  
                  crFnBox = transformation.transformBack(crFnBox);
                  if (crFnBox.intersectsNotEmpty(crFaceBox))
                    {
                      Box intersect = crFnBox & crFaceBox;
                      // Find the cells just inside the original coarse box
                      intersect = adjCellLo(intersect, dir);
                      hiDstBoxes.push_back(intersect);
                      Box otherSide = transformation.transformFwd(intersect);
                      hiSrcBoxes.push_back(otherSide);

                      hiSrcProcs.push_back(srcProc);
                      hiDstProcs.push_back(dstProc);
                    }
                }
            }
        }

      // Now that we have source/destination boxes for coarsened fine faces
      // in this direction, we create a copier that will traffick the data
      // from one to the other.

      // Low side
      {
        map<Box, Box> srcToDstMapping;
        const int bEnd = loSrcBoxes.size();
        for (int b = 0; b != bEnd; ++b)
          {
            srcToDstMapping[loSrcBoxes[b]] = loDstBoxes[b];
          }
        DisjointBoxLayout srcGrid(loSrcBoxes, loSrcProcs);
        DisjointBoxLayout dstGrid(loDstBoxes, loDstProcs);
        m_remoteCopiers[dir][0].define(crFnGrid, // a_from
                                       a_CrGrid, // a_to
                                       srcGrid, // src
                                       dstGrid, // dest
                                       srcToDstMapping,
                                       IntVect::Zero);
      }
      // High side
      {
        map<Box, Box> srcToDstMapping;
        const int bEnd = hiSrcBoxes.size();
        for (int b = 0; b != bEnd; ++b)
          {
            srcToDstMapping[hiSrcBoxes[b]] = hiDstBoxes[b];
          }
        DisjointBoxLayout srcGrid(hiSrcBoxes, hiSrcProcs);
        DisjointBoxLayout dstGrid(hiDstBoxes, hiDstProcs);
        m_remoteCopiers[dir][1].define(crFnGrid, // a_from
                                       a_CrGrid, // a_to
                                       srcGrid, // src
                                       dstGrid, // dest
                                       srcToDstMapping,
                                       IntVect::Zero);
      }
    }  // Loop over directions
}

/*--------------------------------------------------------------------*/
//  Operator applied during copyTo
/**
*//*------------------------------------------------------------------*/

class MBAddFaceOp: public LDOperator<FArrayBox>
{
public:

//--Constructors/destructors

  // Construct
  MBAddFaceOp(const int a_dir,
              const Side::LoHiSide a_side,
              const MultiBlockCoordSys* const a_coordSysPtr,
              const Real a_scale = 1.,
              const bool a_directional = true)
    :
    m_dir(a_dir),
    m_side(a_side),
    m_coordSysPtr(a_coordSysPtr),
    m_directional(a_directional),
    m_scale(a_scale)
  {  }

  // Use synthesized destructor, copy, and assignment

//--Member functions

  /// Operation for distributed memory (src is the buffer)
  virtual void linearIn(FArrayBox& arg,
                        void* buf,
                        const Box& regionTo,
                        const Interval& comps) const
    {
      Real* buffer = (Real*)buf;

      IndicesTransformation itBack = backTransform(arg);
      Box regionFrom = itBack.transformFwd(regionTo);

      FArrayBox tmp(regionFrom, comps.size(), buffer);
      op(arg, regionFrom, comps, regionTo, tmp, comps);
    }

  /// Operation for local memory
  void op(FArrayBox&       a_dest,
          const Box&       a_regionFrom,
          const Interval&  a_Cdest,
          const Box&       a_regionTo,
          const FArrayBox& a_src,
          const Interval&  a_Csrc) const
    {
      CH_assert(a_regionFrom.numPts() == a_regionTo.numPts());
      // Shift to align the face with cells.
      // If low side, then shift each m_dir-face down to the cell below.
      // If high side, then shift each m_dir-face up to the cell above.
      // We'll shift back when we're done.
      a_dest.shiftHalf(m_dir, -sign(m_side));

      CH_assert(a_dest.box().contains(a_regionTo));
      CH_assert(a_src.box().contains(a_regionFrom));

      IndicesTransformation itBack = backTransform(a_dest);

      // Get an understanding of sign transformations in the 'to space'.  Signs
      // are not considered if 'm_directional' is set to false
      //**CHECKME Using transformVectorBack works but is weird.  I guess I don't
      //**        understand why the sign changes based on the direction
      const IntVect signTrfm = (m_directional) ?
        itBack.transformVectorBack(IntVect::Unit) : IntVect::Unit;

      CH_assert(a_Csrc.size() == a_Cdest.size());
      int ncomp = a_Csrc.size();
      if (ncomp % SpaceDim == 0)
        {

//--Not only do we have to transform the to/from locations, but the components
//--themselves are assumed to represent vector data and must be transformed.
          const int compPerDir = ncomp/SpaceDim;
          IntVect compIVDest0;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              compIVDest0[dir] = dir*compPerDir;
            }
          IntVect compIVSrc0 = itBack.transformVectorFwd(compIVDest0);

          if (compPerDir == SpaceDim)
            {

//--Tensor, meaning we also have to transform each range of components

              IntVect compIVDest1(D_DECL6(0, 1, 2, 3, 4, 5));
              IntVect compIVSrc1 = itBack.transformVectorFwd(compIVDest1);

              // For each direction of components
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  // Over all components in a given direction
                  for (int comp = 0; comp != compPerDir; ++comp)
                    {
                      const int compCdest =
                        a_Cdest.begin() + compIVDest0[dir] + compIVDest1[comp];
                      // Sign changes taken from the inner vector transform
                      const Real scale = signTrfm[comp]*m_scale;
                      const int compCsrc = a_Csrc.begin() +
                        abs(compIVSrc0[dir]) + abs(compIVSrc1[comp]);
                      for (BoxIterator bit(a_regionTo); bit.ok(); ++bit)
                        {
                          IntVect ivTo = bit();
                          IntVect ivFrom = itBack.transformFwd(ivTo);
                          a_dest(ivTo, compCdest) +=
                            scale * a_src(ivFrom, compCsrc);
                        }
                    }
                }
            }
          else
            {

//--Vector only (no additional vector meaning within compPerDir)

              // For each direction of components
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  // Over all components in a given direction
                  int compCdest = a_Cdest.begin() + compIVDest0[dir];
                  // Sign changes
                  const Real scale = signTrfm[dir]*m_scale;
                  int compCsrc = a_Csrc.begin() + abs(compIVSrc0[dir]);
                  for (int comp = 0; comp != compPerDir; ++comp)
                    {
                      for (BoxIterator bit(a_regionTo); bit.ok(); ++bit)
                        {
                          IntVect ivTo = bit();
                          IntVect ivFrom = itBack.transformFwd(ivTo);
                          a_dest(ivTo, compCdest) +=
                            scale * a_src(ivFrom, compCsrc);
                        }
                      ++compCsrc;
                      ++compCdest;
                    }
                }
            }
        }
      else
        {

//--The components are assumed to be vector data associated with 'm_dir'.
//--We can directly copy from 'from' direction to 'to direction' and only have
//--to worry about the sign change associated with 'm_dir'
//--locations are transformed.

          int compCsrc = a_Csrc.begin();
          int compCdest = a_Cdest.begin();
          // Sign changes
          const Real scale = signTrfm[m_dir]*m_scale;
          for (int comp = 0; comp < ncomp; comp++)
            {
              for (BoxIterator bit(a_regionTo); bit.ok(); ++bit)
                {
                  IntVect ivTo = bit();
                  IntVect ivFrom = itBack.transformFwd(ivTo);
                  a_dest(ivTo, compCdest) += scale * a_src(ivFrom, compCsrc);
                }
              ++compCsrc;
              ++compCdest;
            }
        }

      // Shift back to original
      a_dest.shiftHalf(m_dir, sign(m_side));
    }

private:

  IndicesTransformation backTransform(const FArrayBox&  a_destFab) const
  {
    // Call whichBlockOverlap() on the FABs' boxes,
    // instead of whichBlock() on a_regionFrom and a_regionTo,
    // because we're going to have cells OUTSIDE valid blocks.
    // int blockFrom = m_coordSysPtr->whichBlockOverlap(a_src.box());
    int blockTo = m_coordSysPtr->whichBlockOverlap(a_destFab.box());
    // petermc, 14 Nov 2011: Do not obtain the IndicesTransformation
    // this way, because there may be more than one way to get from
    // source block to dest block!
    // const IndicesTransformation& it =
    // m_coordSysPtr->blockTransformation(blockFrom, blockTo);
    int faceID = m_dir;
    if (m_side == Side::Hi) faceID += SpaceDim;
    const BlockBoundary& bb = m_coordSysPtr->boundaries()[blockTo][faceID];
    IndicesTransformation itBack = bb.getTransformation();
    return itBack;
  }

  //--Data members

  int m_dir;                          ///< The dir this operation affects
  Side::LoHiSide m_side;              ///< The side this operation affects
  const MultiBlockCoordSys* m_coordSysPtr;
  bool m_directional;                 ///< T - Vector data on faces has a
                                      ///<     direction (normal & default)
                                      ///< F - The vector data does not have
                                      ///<     a direction.  The best example of
                                      ///<     this is <N> which has vector
                                      ///<     components but really just scales
                                      ///<     areas.

public:
  Real m_scale;                       ///< Extra scaling beyond what is normally
                                      ///< needed to integrated the fine fluxes
};

/*--------------------------------------------------------------------*/
//  Overwrite solution on the coarse faces
/** \param[out] a_CrF   Coarse information updated with face
 *                      corrections from the finer level
 *  \param[in]  a_CrFIntv
 *                      Interval in a_CrF to update
 *  \param[in]  a_faceIntv
 *                      Interval in interior structures
 *  \param[in]  a_scale Additional scaling
 *  \param[in]  a_directional
 *                      T - Vector data on faces has a direction
 *                          (normal & default)
 *                      F - The vector data does not have a direction.
 *                          The best example of this is <N> which has
 *                          vector components but really just scales
 *                          areas.
*//*------------------------------------------------------------------*/

void
MultiBlockFaceRegister::reface(LevelData<FluxBox>& a_CrF,
                               const Interval&     a_CrFIntv,
                               const Interval&     a_faceIntv,
                               const Real          a_scale,
                               const bool          a_directional)
{
  CH_TIME("MultiBlockFaceRegister::reface");
  CH_assert(m_coordSys->mappingBlocks().size() > 1);
  CH_assert(m_isDefined & FluxRegFineDefined);

  // Initialize coarse locations to zero
  for (DataIterator dit = a_CrF.dataIterator(); dit.ok(); ++dit)
    {
      FluxBox& CrFFlbx = a_CrF[dit];
      int faceID = 0;
      // Here, side refers to the side of the coarsened-fine box, which
      // is opposite the side of the coarse box.  So be careful when
      // shifting.
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          int sideSign = sign(side);
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              FArrayBox& CrFFabDir = CrFFlbx[idir];

              const Vector<Box>& intersect = m_coarseLocations[faceID][dit];
              const int bEnd = intersect.size();
              for (int b = 0; b < bEnd; ++b)
                {
                  Box faceIntersect(intersect[b]);
                  faceIntersect.shiftHalf(idir, -sideSign);
                  CrFFabDir.setVal(0.,
                                   faceIntersect,
                                   a_CrFIntv.begin(),
                                   a_CrFIntv.size());
                }
              faceID++;
            }
        }
    }

  // Copy all the fine faces from their "native" locations to those abutting
  // the coarse faces.
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      // Alias a direction of the flux box with the original inteval
      LevelData<FArrayBox> CrFDirAlias(
        a_CrF.getBoxes(),
        a_CrF.nComp(),
        a_CrF.ghostVect(),
        FABAliasFlBxDataFactory(&a_CrF, a_CrF.interval(), dir));

      int sideInt = 0;
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          if (m_remoteCopiers[dir][sideInt].isDefined())
            {
              // Copy the coarsened fine faces to the corresponding coarse faces
              // on their neighboring blocks.
              MBAddFaceOp op(dir, side, m_coordSys, -a_scale, a_directional);
              m_fineFlux.copyTo(a_faceIntv,
                                CrFDirAlias,
                                a_CrFIntv,
                                m_remoteCopiers[dir][sideInt],
                                op);
            }
          sideInt++;
        } // end iteration over side
    } // end iteration over dir
}

#include "NamespaceFooter.H"
