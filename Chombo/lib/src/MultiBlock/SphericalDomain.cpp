#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphericalDomain.H"
#include "RectangularUniformMap.H"
#include "SphericalSectorMap.H"
// #include "UGIO.H"
// #include "CH_HDF5.H"

#include "NamespaceHeader.H"


//=====================================================================

SphericalDomain::SphericalDomain()
{
  m_isDefined = false;
}


// ---------------------------------------------------------
/// Construct SphericalDomain with a single \em a_domBox as computational domain
SphericalDomain::SphericalDomain(const Box&        a_centralBox,
                                 const RealVect&   a_center,
                                 Real              a_bxWidth,
                                 Real              a_outerRadius)
{
  define(a_centralBox, a_center, a_bxWidth, a_outerRadius);
}


// ---------------------------------------------------------
SphericalDomain::~SphericalDomain()
{
}


// ---------------------------------------------------------
void SphericalDomain::define(const Box&        a_centralBox,
                             const RealVect&   a_center,
                             Real              a_bxWidth,
                             Real              a_outerRadius)
{
#if CH_SPACEDIM == 3
  m_nblocks = 7;
  m_bxWidth = a_bxWidth;
  m_outerRadius = a_outerRadius;
  m_center = a_center;
  m_centralCornerLo = m_center - (m_bxWidth/2.) * RealVect::Unit;
  m_centralCornerHi = m_center + (m_bxWidth/2.) * RealVect::Unit;

  /*
    Define block domains.  Each is a shifted copy of the central domain.
   */
  Tuple<IntVect, SpaceDim> shifts;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      shifts[idir] = a_centralBox.size(idir) * BASISV(idir);
    }

  // The full domain.
  m_domainBox = grow(a_centralBox, a_centralBox.size());

  /*
    Set m_axis and m_direction for the sectoral blocks.
  */
  m_axis[0] = 0; m_direction[0] = Side::Hi; // +x
  m_axis[1] = 1; m_direction[1] = Side::Hi; // +y
  m_axis[2] = 0; m_direction[2] = Side::Lo; // -x
  m_axis[3] = 1; m_direction[3] = Side::Lo; // -y
  m_axis[4] = 2; m_direction[4] = Side::Hi; // +z
  m_axis[5] = 2; m_direction[5] = Side::Lo; // -z

  Tuple<int, SpaceDim> blockLo, blockHi;
  for (int iblock = 0; iblock < 2*SpaceDim; iblock++)
    {
      int axis = m_axis[iblock];
      switch (m_direction[iblock])
        {
        case Side::Lo :
          {
            blockLo[axis] = iblock;
            break;
          }
        case Side::Hi :
          {
            blockHi[axis] = iblock;
            break;
          }
        default:
          {
          }
        }
    }

  /*
    Define the domain of each block.
   */
  m_boxes.resize(m_nblocks);
  m_boxes[CUBE] = a_centralBox;
  for (int iblock = 0; iblock < 2*SpaceDim; iblock++)
    {
      m_boxes[iblock] = a_centralBox +
        sign(m_direction[iblock]) * shifts[m_axis[iblock]];
    }

  /*
    Define Tuple<BlockBoundary, 2*SpaceDim> m_blockBoundaries[iblock]
    for iblock in 0:m_nblocks-1
   */
  m_blockBoundaries.resize(m_nblocks);

  // Set boundaries of central cube and outer faces of sectoral blocks.
  for (int iblock = 0; iblock < 2*SpaceDim; iblock++)
    {
      // which face of central cube abuts sectoral block iblock
      int faceID = m_axis[iblock];
      if (m_direction[iblock] == Side::Hi) faceID += SpaceDim;

      // which face of sectoral block iblock abuts central cube
      int faceOpposite = (faceID < SpaceDim) ?
        faceID + SpaceDim : faceID - SpaceDim;

      // from CUBE (on faceID) to iblock (on faceOpposite)
      m_blockBoundaries[CUBE][faceID].defineConformal(iblock, faceOpposite);
      // from iblock (on faceOpposite) to CUBE (on faceID)
      m_blockBoundaries[iblock][faceOpposite].defineConformal(CUBE, faceID);

      // from iblock (on faceID) to nothing
      m_blockBoundaries[iblock][faceID].define(0);
    }

  // Set boundaries of other faces of sectoral blocks.
  for (int iblock = 0; iblock < 2*SpaceDim; iblock++)
    {
      int axis = m_axis[iblock];
      Side::LoHiSide direction = m_direction[iblock];
      // faceOther:  dest face on iblockOther
      int faceOther = axis;
      if (direction == Side::Hi) faceOther += SpaceDim;
      for (int axisOther = 0; axisOther < SpaceDim; axisOther++)
        if (axisOther != axis)
          {
            for (SideIterator sit; sit.ok(); sit.next())
              {
                Side::LoHiSide directionOther = sit();
                // iblockOther:  dest block
                int iblockOther = (directionOther == Side::Lo) ?
                  blockLo[axisOther] : blockHi[axisOther];
                // face:  source face on iblock
                int face = axisOther;
                if (directionOther == Side::Hi) face += SpaceDim;
                // FROM face on iblock
                // TO faceOther on iblockOther

                // Recall, for each direction idir:
                // pnew[idir] = m_sign[idir]*pold[m_permutation[idir]]
                // + m_translation[idir]

                // Follow each direction:
                // [iblock] (pold)             -> [iblockOther] (pnew)
                // axis, direction             -> axisOther, directionOther
                // axisOther, directionOther   -> axis, -direction
                // axisRemaining, +            -> axisRemaining, +

                IntVect perm(D_DECL(0, 1, 2)); // keep third dimension fixed
                perm[axisOther] = axis;
                perm[axis] = axisOther;

                IntVect sgn = IntVect::Unit; // +1 in third dimension
                sgn[axisOther] = sign(directionOther) * sign(direction);
                sgn[axis] = sign(directionOther) * -sign(direction);

                // The common edge between the two blocks remains fixed:
                // p[axis] = corner
                // p[axisOther] = cornerOther
                // added by petermc, 27 May 2008
                Box centralBoxNodes = surroundingNodes(a_centralBox);
                int cornerNodes =
                  centralBoxNodes.sideEnd(direction)[axis];
                int cornerNodesOther =
                  centralBoxNodes.sideEnd(directionOther)[axisOther];
                // What happens when we move away from the common edge?
                // On the shared face:
                // direction * (pold[axis] - corner) =
                // directionOther * (pnew[axisOther] - cornerOther)
                // Perpendicular to the shared face:
                // directionOther * (pold[axisOther] - cornerOther) =
                // -direction * (pnew[axis] - corner)

                // Hence
                // pnew[axis] = corner - direction*directionOther *
                //              (pold[axisOther] - cornerOther)
                //            = (-direction*directionOther) * pold[paxisOther]
                //              + corner + direction*directionOther*cornerOther
                // pnew[axisOther] = cornerOther + directionOther*direction *
                //                   (pold[axis] - corner)
                //                 = (directionOther*direction) * pold[axis]
                //                   + cornerOther
                //                   - directionOther*direction*corner
                // Hence trans[axis] =
                //       corner + direction*directionOther*cornerOther
                // and trans[axisOther] =
                //       cornerOther - directionOther*direction*corner
                IntVect trans(IntVect::Zero); // 0 in third dimension
                // correction for signCombined added by petermc, 27 May 2008
                // If adding, need to subtract 1 from sum.
                int signCombined = sign(direction)*sign(directionOther);
                trans[axis] = cornerNodes + signCombined*cornerNodesOther;
                if (signCombined == +1) trans[axis] -= 1;
                trans[axisOther] = cornerNodesOther - signCombined*cornerNodes;
                if (-signCombined == +1) trans[axisOther] -= 1;
                m_blockBoundaries[iblock][face].define(perm, sgn, trans,
                                                       iblockOther, faceOther);
              }
          }
    }

  /*
    Define m_blockMaps[iblock] for iblock in 0:m_nblocks-1.
   */
  m_blockMaps.resize(m_nblocks);
  m_blockMaps[CUBE] = new RectangularUniformMap(m_centralCornerLo,
                                                m_centralCornerHi,
                                                m_boxes[CUBE]);
  for (int iblock = 0; iblock < 2*SpaceDim; iblock++) // sector blocks only
    {
      m_blockMaps[iblock] = new SphericalSectorMap(m_center,
                                                   m_bxWidth,
                                                   m_outerRadius,
                                                   iblock,
                                                   m_boxes[iblock]);
    }

  /*
    Define MappedBlock m_mappedBlock[iblock] for iblock in 0:m_nblocks-1
   */
  m_mappedBlocks.resize(m_nblocks);
  for (int iblock = 0; iblock < m_nblocks; iblock++) // ALL blocks
    {
      m_mappedBlocks[iblock].define(m_boxes[iblock],
                                    m_blockMaps[iblock],
                                    m_blockBoundaries[iblock]);
    }

  defineMappedDomain();
#endif
}

#include "NamespaceFooter.H"
