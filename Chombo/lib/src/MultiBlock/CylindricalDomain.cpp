#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CylindricalDomain.H"
// #include "UGIO.H"
// #include "CH_HDF5.H"

#include "NamespaceHeader.H"


//=====================================================================

CylindricalDomain::CylindricalDomain()
{
  m_isDefined = false;
}


// ---------------------------------------------------------
/// Construct CylindricalDomain with a single \em a_domBox as computational domain
CylindricalDomain::CylindricalDomain(const Box&        a_centralBox,
                                     const RealVect&   a_center,
                                     Real              a_bxWidth,
                                     Real              a_outerRadius)
{
  define(a_centralBox, a_center, a_bxWidth, a_outerRadius);
}


// ---------------------------------------------------------
CylindricalDomain::~CylindricalDomain ()
{
}


// ---------------------------------------------------------
void CylindricalDomain::define(const Box&        a_centralBox,
                               const RealVect&   a_center,
                               Real              a_bxWidth,
                               Real              a_outerRadius)
{
  m_nblocks = 5;
  m_bxWidth = a_bxWidth;
  m_outerRadius = a_outerRadius;
  m_center = a_center;
  m_centralCornerLo = m_center - (m_bxWidth/2.) * RealVect::Unit;
  m_centralCornerHi = m_center + (m_bxWidth/2.) * RealVect::Unit;

  /*
    Define block domains.  Each is a shifted copy of the central domain.
   */
  IntVect sh0 = a_centralBox.size(0) * BASISV(0);
  IntVect sh1 = a_centralBox.size(1) * BASISV(1);
  // The full domain.
  m_domainBox = a_centralBox;
  m_domainBox.grow(sh0 + sh1);

  m_boxes.resize(m_nblocks);
  m_boxes[CUBE] = a_centralBox;
  m_boxes[1] = a_centralBox + sh0;
  m_boxes[2] = a_centralBox + sh1;
  m_boxes[3] = a_centralBox - sh0;
  m_boxes[4] = a_centralBox - sh1;

  /*
    Define Tuple<BlockBoundary, 2*SpaceDim> m_blockBoundaries[iblock]
    for iblock in 0:m_nblocks-1
   */
  m_blockBoundaries.resize(m_nblocks);
  // This is long and tedious.
#if CH_SPACEDIM == 3
  // Going in direction +z or -z from ANY block hits physical boundary.
  for (int iblock = 0; iblock < m_nblocks; iblock++) // ALL blocks
    {
      m_blockBoundaries[iblock][2].define(0);
      m_blockBoundaries[iblock][2 + SpaceDim].define(0);
    }
#endif
  // Boundaries of central cube and outer faces.
  // Order of arguments:  neighbor, face.
  setCubeAndOuterBoundaries(1, 0 + SpaceDim);
  setCubeAndOuterBoundaries(2, 1 + SpaceDim);
  setCubeAndOuterBoundaries(3, 0);
  setCubeAndOuterBoundaries(4, 1);

  // Order:  srcblock, dstblock, srcface, dstface, corner0side, corner1side
  setInnerBoundary(1, 2, 1+SpaceDim, 0+SpaceDim, Side::Hi, Side::Hi);
  setInnerBoundary(2, 3, 0,          1+SpaceDim, Side::Lo, Side::Hi);
  setInnerBoundary(3, 4, 1,          0,          Side::Lo, Side::Lo);
  setInnerBoundary(4, 1, 0+SpaceDim, 1,          Side::Hi, Side::Lo);

  // MappedDomain will fill in the mates of these last four.

  m_blockMaps.resize(m_nblocks);
  m_mappedBlocks.resize(m_nblocks);
}


// ---------------------------------------------------------
void CylindricalDomain::setCubeAndOuterBoundaries(int a_neighbor,
                                                  int a_face)
{
  // defineConformal(neighbor, face)
  int faceOpposite = (a_face < SpaceDim) ?
    a_face + SpaceDim : a_face - SpaceDim;
  // from CUBE (on a_face) to a_neighbor (on faceOpposite)
  m_blockBoundaries[CUBE][a_face].defineConformal(a_neighbor, faceOpposite);
  // from a_neighbor (on faceOpposite) to CUBE (on a_face)
  m_blockBoundaries[a_neighbor][faceOpposite].defineConformal(CUBE, a_face);
  // from a_neighbor (on a_face) to nothing
  m_blockBoundaries[a_neighbor][a_face].define(0);
}


// ---------------------------------------------------------
void CylindricalDomain::setInnerBoundary(int a_srcBlock,
                                         int a_dstBlock,
                                         int a_srcFace,
                                         int a_dstFace,
                                         Side::LoHiSide a_corner0side,
                                         Side::LoHiSide a_corner1side)
{
  Box centralNodes = surroundingNodes(m_boxes[CUBE]);

  // In every case, we swap dimensions.
  IntVect swapDims(D_DECL(1, 0, 2));

  // first two coordinates of corner NODE around which we pivot
  int corn0 = centralNodes.sideEnd(a_corner0side)[0];
  int corn1 = centralNodes.sideEnd(a_corner1side)[1];

  IntVect reflector;
  IntVect trans;
  if ((a_dstBlock == a_srcBlock + 1) || (a_dstBlock == a_srcBlock - 3))
    {
      // Going counterclockwise:
      // new coordinate 0 = opposite-sign old coordinate 1
      // new coordinate 1 = same-sign old coordinate 0
      reflector = IntVect(D_DECL(-1, 1, 1));
      trans = IntVect(D_DECL(corn0 + corn1 - 1, corn1 - corn0, 0));
    }
  else if ((a_dstBlock == a_srcBlock - 1) || (a_dstBlock == a_srcBlock + 3))
    {
      // Going clockwise:
      // new coordinate 0 = same-sign old coordinate 1
      // new coordinate 1 = opposite-sign old coordinate 0
      reflector = IntVect(D_DECL(1, -1, 1));
      trans = IntVect(D_DECL(corn0 - corn1, corn1 + corn0 - 1, 0));
    }
  else
    {
      MayDay::Error("setInnerBoundary:  blocks are not adjacent");
    }

  m_blockBoundaries[a_srcBlock][a_srcFace].define(swapDims, reflector, trans,
                                                  a_dstBlock, a_dstFace);
}

#include "NamespaceFooter.H"
