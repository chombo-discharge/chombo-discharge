#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "XPointCoordSys.H"
#include "XPointBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "BlockBoundary.H"
#include "SpreadingCopier.H"
#include "ReductionOps.H"
#include "CONSTANTS.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"


XPointCoordSys::XPointCoordSys()
{
  //  ParmParse pp("xpsystem");
  //
  //  Vector<Real> decomp(SpaceDim);
  //  pp.getarr("block_decomp", decomp, 0, SpaceDim);
  //  for (int i=0; i<SpaceDim; ++i) {
  //    m_decomp[i] = decomp[i];
  //  }
}


XPointCoordSys::XPointCoordSys(const ProblemDomain& a_levelDomain,
                               const RealVect& a_dx)
{
  define(a_levelDomain, a_dx);
}


void
XPointCoordSys::define(const ProblemDomain& a_levelDomain,
                       const RealVect& a_dx)
{
  // a_levelDomain is the whole domain of indices used by all blocks.
  // Just take the height.
  const Box& levelDomainBox = a_levelDomain.domainBox();
  int baseLength = levelDomainBox.size(1);
  CH_assert(levelDomainBox.size(0) == 15 * baseLength);
  CH_assert(levelDomainBox.size(1) == baseLength);
  CH_assert(levelDomainBox.smallEnd() == IntVect::Zero);

  // In mapped space, the blocks will be, in order:
  // [0:1, 0:1], [2:3, 0:1], [4:5, 0:1], [6:7, 0:1],
  // [8:9, 0:1], [10:11, 0:1], [12:13, 0:1], [14:15, 0:1].
  m_mappingBlocks.resize(8);
  m_coordSysVect.resize(8, NULL);

  // We don't even use the input a_dx, because the mesh spacing is
  // determined by the height of a_levelDomain.  Each block is 1 x 1
  // in mapped space, so dx == 1/baseLength .

  IntVect baseLo = IntVect::Zero;
  IntVect baseHi = (baseLength - 1) * IntVect::Unit;
  Box baseBox(baseLo, baseHi);
  IntVect shiftNextBox = (2*baseLength) * BASISV(0);

  RealVect originBlock = RealVect::Zero;
  ProblemDomain domainBlock(baseBox);
  for (int iblock = 0; iblock < 8; iblock++)
    {
      m_mappingBlocks[iblock] = domainBlock.domainBox();
      m_coordSysVect[iblock] =
        new XPointBlockCoordSys(originBlock, domainBlock, iblock);
      // originBlock += 2.0 * RealVect(BASISV(0));
      domainBlock.shift(shiftNextBox);
    }

  m_gotMappingBlocks = true;
  m_gotCoordSysVect = true;

  /*
    Define block boundaries.
   */
  defineBoundaries();

  // Define and fill in m_validNum, m_validBlocks, m_validTransformations.
  initializeBlockTransformations();
}


XPointCoordSys::~XPointCoordSys()
{
  if (m_gotCoordSysVect)
    {
      for (int i=0; i<8; i++)
        {
          delete m_coordSysVect[i];
        }
    }
}


void 
XPointCoordSys::defineBoundaries()
{
   CH_assert(SpaceDim==2);
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   IntVect permDiff = IntVect(D_DECL6(1, 0, 2,0,0,0));
   IntVect clockwise = IntVect(D_DECL6(-1, 1, 1,0,0,0));
   IntVect anticlockwise = IntVect(D_DECL6(1, -1, 1,0,0,0));

   for (int iblock = 0; iblock < numBlocks(); iblock++) {

      IndicesTransformation it;
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];
      const XPointBlockCoordSys* coord_sys_src = (const XPointBlockCoordSys*)m_coordSysVect[iblock];
      // int block_type = coord_sys_src->blockType();
      int neighbor_block[SpaceDim];

      neighbor_block[0] = (iblock+numBlocks()-1) % numBlocks();
      neighbor_block[1] = (iblock+1) % numBlocks();
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         const XPointBlockCoordSys* coord_sys_dst = (const XPointBlockCoordSys*)m_coordSysVect[neighbor_block[dir]];
         IntVect rotation;
         if (dir==0) {
            rotation = clockwise;
         }
         else {
            rotation = anticlockwise;
         }
         IntVect pold = coord_sys_src->domainSmallEnd();
         IntVect pnew = coord_sys_dst->domainSmallEnd();
         it.defineFromPivot(pold, pnew, permDiff, rotation);
         blockBoundaries[dir].define(it,neighbor_block[dir]);
         blockBoundaries[dir + SpaceDim].define(0);
      }
   }
   m_gotBoundaries = true;
}



void 
XPointCoordSys::blockRemapping(RealVect&       a_xi_valid,
                               int&            a_n_valid,
                               const RealVect& a_xiSrc,
                               int             a_nSrc) const
{
   const XPointBlockCoordSys* coord_sys_Src = (const XPointBlockCoordSys*)m_coordSysVect[a_nSrc];
   RealVect X = coord_sys_Src->realCoord(a_xiSrc);

   Real R = X[0] - coord_sys_Src->centralMajorRadius();
   Real Z = X[1];

   if (R >= 0. && Z > 0.) {     // upper right quadrant
      if (Z > R) {
         a_n_valid = 0;
      }
      else {
         a_n_valid = 1;
      }
   }
   else if (R > 0. && Z <= 0.) { // lower right quadrant
      if (Z > -R) {
         a_n_valid = 2;
      }
      else {
         a_n_valid = 3;
      }
   }
   else if (R <= 0. && Z < 0.) { // lower left quadrant
      if (Z < R) {
         a_n_valid = 4;
      }
      else {
         a_n_valid = 5;
      }
   }
   else if (R < 0. && Z >= 0.) { // upper left quadrant
      if (Z < -R) {
         a_n_valid = 6;
      }
      else {
         a_n_valid = 7;
      }
   }

   const XPointBlockCoordSys* coord_sys_valid = (const XPointBlockCoordSys*)m_coordSysVect[a_n_valid];
   a_xi_valid = coord_sys_valid->mappedCoord(X);

   // Perform a final check that a_xi_valid is contained in the valid region of the block whose
   // index is being returned
   if ( coord_sys_valid->lowerMappedCoordinate(0) > a_xi_valid[0] ||
        coord_sys_valid->upperMappedCoordinate(0) < a_xi_valid[0] ||
        coord_sys_valid->lowerMappedCoordinate(1) > a_xi_valid[1] ||
        coord_sys_valid->upperMappedCoordinate(1) < a_xi_valid[1] ) {
      a_n_valid = -1;
      a_xi_valid = a_xiSrc;
   }
}


MultiBlockCoordSys*
XPointCoordSysFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                   const RealVect& a_dx) const
{
  XPointCoordSys* coordSysPtr = new XPointCoordSys();
  coordSysPtr->define(a_levelDomain, a_dx);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
