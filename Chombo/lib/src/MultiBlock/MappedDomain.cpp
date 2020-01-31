#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MappedDomain.H"
#include "BoxIterator.H"
#include "UGIO.H"
#include "CH_HDF5.H"

#include "NamespaceHeader.H"


//=====================================================================
MappedDomain::MappedDomain ()
{
}

// ---------------------------------------------------------
/// Construct MappedDomain with a single \em a_domBox as computational domain
/**
   This constructor defaults to non-periodic domain
*/
MappedDomain::MappedDomain(const Box& a_domBox,
                           int a_defaultBC) // void* a_defaultBC)
{
  CH_assert(!a_domBox.isEmpty());
  CH_assert(a_defaultBC!= -1); // !=NULL
  // change this to Tuple?
  //  Vector<BlockBoundary> a_boundaries(2*CH_SPACEDIM, BlockBoundary(a_defaultBC));
  Tuple<BlockBoundary, 2*SpaceDim> boundaries;
  for (int faceID = 0; faceID < 2*SpaceDim; faceID++)
    {
      boundaries[faceID] = BlockBoundary(a_defaultBC);
    }
  m_blocks.resize(1, MappedBlock(a_domBox, &BlockMap::Identity, boundaries));
  m_domain = a_domBox;
}


// ---------------------------------------------------------
MappedDomain::MappedDomain(const Vector<MappedBlock>& a_blocks,
                           int a_defaultBC) // void* a_defaultBC
{
  define(a_blocks, a_defaultBC);
}


// ---------------------------------------------------------
void MappedDomain::define(const Vector<MappedBlock>& a_blocks,
                          int a_defaultBC) // void* a_defaultBC
{
  // OK, the user has handed us a vector of MappedBlocks,
  // we need to check things
  // out inside, and clean up dangling connections.

  m_blocks = a_blocks;
  int nblocks = m_blocks.size();

  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      if (!m_blocks[iblock].hasMap())
        {
          MayDay::Abort("MappedBlock object used lacking BlockMap. possibly null constructed");
        }

      const Box& blockBox = m_blocks[iblock].box();
      if (iblock == 0)
        {
          m_domain = blockBox;
        }
      else
        {
          // m_domain |= blockBox;
          IntVect domainLo, domainHi;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              domainLo[idir] =
                Min(m_domain.smallEnd(idir), blockBox.smallEnd(idir));
              domainHi[idir] =
                Max(m_domain.bigEnd(idir), blockBox.bigEnd(idir));
            }
          m_domain = Box(domainLo, domainHi);
        }

      Tuple<BlockBoundary, 2*SpaceDim>& boundaries = m_blocks[iblock].m_boundaries;
      for (int faceID = 0; faceID < 2*SpaceDim; faceID++)
        {
          BlockBoundary& bb = boundaries[faceID];
          if (bb.type() == BlockBoundary::BOUNDARY)
            {
              // that's cool, I trust you are doing a boundary condition correctly.
            }
          else if (bb.type() == BlockBoundary::CONFORMAL ||
                   bb.type() == BlockBoundary::MAPPED)
            {
              // need to check that the other block matches up
              int blockOther = bb.m_neighbor;
              int faceOther  = bb.m_face;
              BlockBoundary& mate =
                m_blocks[blockOther].m_boundaries[faceOther];
              if (mate.type() == BlockBoundary::UNKNOWN)
                {
                  mate.m_neighbor = iblock;
                  mate.m_face = faceID;
                  mate.m_bc = -1; // NULL;
                  mate.m_type = bb.m_type;
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      int dirOther = bb.m_permutation[idir];
                      mate.m_permutation[dirOther] = idir;
                      mate.m_sign[dirOther]        = bb.m_sign[idir];
                      mate.m_translation[dirOther] =
                        -bb.m_sign[idir] * bb.m_translation[idir];
                    }
                }

            }
        }
    }
  // now, all implicitly mated surfaces should be matched up, no boundary should be UNKNOWN now.

  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      Tuple<BlockBoundary, 2*SpaceDim>& boundaries = m_blocks[iblock].m_boundaries;
      for (int b=0; b<2*CH_SPACEDIM; ++b)
        {
          BlockBoundary& bb = boundaries[b];
          if (bb.type() == BlockBoundary::UNKNOWN)
            {
              if (a_defaultBC == -1)
              {
                // was NULL
                MayDay::Abort("UNKNOWN BlockBoundary object remaining");
              }
              else
              {
                bb = BlockBoundary(a_defaultBC);
              }
            }
        }
    }

}


// ---------------------------------------------------------
int MappedDomain::findBlock(const Box&   a_bx) const
{
  int nblocks = m_blocks.size();
  int thisBlock = -1;
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      if (m_blocks[iblock].box().contains(a_bx))
        {
          thisBlock = iblock;
          break;
        }
    }
  CH_assert(thisBlock >= 0);
  return thisBlock;
}


// ---------------------------------------------------------
int MappedDomain::findBlock(const RealVect&   a_physPoint) const
{
  int nblocks = m_blocks.size();
  int thisBlock = -1;
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      if (m_blocks[iblock].map().inThisBlock(a_physPoint))
        {
          thisBlock = iblock;
          break;
        }
    }
  CH_assert(thisBlock >= 0);
  return thisBlock;
}


// ---------------------------------------------------------
int MappedDomain::numBlocks() const
{
  return m_blocks.size();
}


// ---------------------------------------------------------
int MappedDomain::blockValid(const IntVect&   a_cell,
                             int              a_block) const
{
  /*
    If it's not in the domain, then return NO_BLOCK.
  */
  if (! m_domain.contains(a_cell))
    {
      return NO_BLOCK;
    }

  /*
    If you find a block with a domain containing a_cell, return that block.
  */
  int nblocks = m_blocks.size();
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      if (m_blocks[iblock].box().contains(a_cell))
        {
          return iblock;
        }
    }

  /*
    a_cell does not belong to the domain of any block.
    So look for a neighbor of a_block that contains it.
  */
  const MappedBlock& mb = m_blocks[a_block];
  for (int faceID = 0; faceID < 2*SpaceDim; faceID++)
    {
      const BlockBoundary& bb = mb.boundary(faceID);
      if (bb.type() == BlockBoundary::CONFORMAL)
        {
          int blockNbr = bb.neighbor();
          // a_cell in blockNbr's index space
          IntVect cellConverted = bb.convertOldToNew(a_cell);
          if (m_blocks[blockNbr].box().contains(cellConverted))
            {
              return blockNbr;
            }
        }
    }

  /*
    a_cell isn't in domain of any block or of any of its neighbors.
    This could happen for example with corner ghost cells of central block.
    In that case, we don't want to use the cell, so NO_BLOCK is correct.
  */
  return NO_BLOCK;
}


// ---------------------------------------------------------
void MappedDomain::getNeighborhood(DenseIntVectSet&   a_ivs,
                                   IntVect&           a_centerCell,
                                   const IntVect&     a_baseCell,
                                   int                a_block,
                                   int                a_radius)
{
  CH_TIME("MappedDomain::neighborhood");
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // indicesCenterWithin[idir] in range minBox:maxBox
      int minBox = m_domain.smallEnd(idir) + a_radius;
      int maxBox = m_domain.bigEnd(idir) - a_radius;
      int orig = a_baseCell[idir];
      if (orig < minBox)
        a_centerCell[idir] = minBox;
      else if (orig > maxBox)
        a_centerCell[idir] = maxBox;
      else // minBox <= orig <= maxBox
        a_centerCell[idir] = orig;
    }
  Box bxNbrs(-a_radius*IntVect::Unit, a_radius*IntVect::Unit);
  bxNbrs.shift(a_centerCell);
  a_ivs = DenseIntVectSet(bxNbrs, false);
  for (BoxIterator bitNbr(bxNbrs); bitNbr.ok(); ++bitNbr)
    {
      IntVect ivNbr = bitNbr();
      int blockNbr = blockValid(ivNbr, a_block);
      if (blockNbr != NO_BLOCK)
        a_ivs |= ivNbr;
    }
}


// ---------------------------------------------------------
IntVect MappedDomain::convertBetweenBlocks(int              a_blockOld,
                                           int              a_blockNew,
                                           const IntVect&   a_indicesOld)
{
  if (a_blockNew == a_blockOld)
    return a_indicesOld;

  const MappedBlock& mbOld = m_blocks[a_blockOld];
  for (int faceID = 0; faceID < 2*SpaceDim; faceID++)
    {
      const BlockBoundary& bb = mbOld.boundary(faceID);
      if (bb.neighbor() == a_blockNew)
        {
          IntVect indicesNew = bb.convertOldToNew(a_indicesOld);
          return indicesNew;
        }
    }
  // If you reach this point then none of the faces matches.
  MayDay::Abort("MappedDomain::convertBetweenBlocks:  blocks do not abut");
  return IntVect::Zero;
}


// ---------------------------------------------------------
bool MappedDomain::convertGhostBlocks(IntVect&         a_indicesNew,
                                      int              a_blockNew,
                                      const IntVect&   a_indicesOld,
                                      int              a_blockOld)
{
  // a_indicesNew in a_blockNew is the same cell as a_indicesOld in a_blockOld.
  a_indicesNew = convertBetweenBlocks(a_blockOld, a_blockNew, a_indicesOld);

  if (blockValid(a_indicesNew, a_blockNew) == NO_BLOCK)
    { // This means that a_indicesNew is not a valid cell of any block,
      // nor is it a valid cell of any neighbor of a_blockNew after
      // conversion to index space of the neighbor from that of a_blockNew.

      // Among neighbors of a_blockOld, find block where a_indicesOld is valid.
      int blockGoodNbr = NO_BLOCK;
      IntVect cellNbr;
      const MappedBlock& mb = m_blocks[a_blockOld];
      for (int faceID = 0; faceID < 2*SpaceDim; faceID++)
        {
          const BlockBoundary& bb = mb.boundary(faceID);
          if (bb.type() == BlockBoundary::CONFORMAL)
            {
              int blockNbr = bb.neighbor();
              // a_indicesOld in a_blockOld's index space;
              // cellNbr in blockNbr's index space
              cellNbr = bb.convertOldToNew(a_indicesOld);
              if (m_blocks[blockNbr].box().contains(cellNbr))
                {
                  blockGoodNbr = blockNbr;
                  break;
                }
            }
        }
      if (blockGoodNbr == NO_BLOCK)
        {
          return false;
        }
      else
        {
          // convert from cellNbr in blockGoodNbr
          // to a_indicesNew in a_blockNew
          a_indicesNew = convertBetweenBlocks(blockGoodNbr,
                                              a_blockNew,
                                              cellNbr);
          return true;
        }
    }
  else
    {
      return true;
    }
}


// ---------------------------------------------------------
MappedDomain refine(const MappedDomain& a_probdomain,
                    int   a_refinement_ratio)
{
  const Vector<MappedBlock>& vmb = a_probdomain.m_blocks;
  int nblocks = vmb.size();
  Vector<MappedBlock> blocksRefined(nblocks);
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      blocksRefined[iblock] =
        refine(vmb[iblock], a_refinement_ratio);
    }
  MappedDomain mdRefined = MappedDomain(blocksRefined);
  return mdRefined;
}

#include "NamespaceFooter.H"
