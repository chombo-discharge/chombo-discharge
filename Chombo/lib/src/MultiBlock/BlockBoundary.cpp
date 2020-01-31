#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BlockBoundary.H"
#include "UGIO.H"
#include "CH_HDF5.H"

#include "NamespaceHeader.H"

//=============================================================

/// null constructor leaves object in type=UNKNOWN state. It can be overridden by MappedDomain
BlockBoundary::BlockBoundary()
{
  m_bc=-1;
  m_type=UNKNOWN;
  m_neighbor=-1;
}


// ---------------------------------------------------------
/// the boundary is an external boundary condition.
/**
   Boundary condition label.  Application is responsible for interpreting
   the pointer.
*/
// BlockBoundary::BlockBoundary(void*  a_boundaryClass)
BlockBoundary::BlockBoundary(int a_boundaryClass)
{
  define(a_boundaryClass);
}


// ---------------------------------------------------------
// void BlockBoundary::define(void*  a_boundaryClass)
void BlockBoundary::define(int a_boundaryClass)
{
  CH_assert(a_boundaryClass != -1); // !=NULL
  m_type=BOUNDARY;
  m_bc=a_boundaryClass;
  m_neighbor=-1;
}


// ---------------------------------------------------------
/// the boundary is a conformal connection between two blocks.
/**
   The BlockMap between these two MappedBlocks are conformal, hence no spatial
   interpolation needs to be performed. However, the index spaces may require a
   transform.

   these topological connections are redundant between two blocks. the user
   can choose to specify all connections, and have MappedDomain verify the
   topology, or the user can provide just the sufficient non-redundant set
   and have MappedDomain construct the mirrors.
*/
BlockBoundary::BlockBoundary(const IntVect& a_permutation,
                             const IntVect& a_sign,
                             const IntVect& a_translation,
                             int a_neighbor,
                             int a_face)
{
  define(a_permutation, a_sign, a_translation, a_neighbor, a_face);
}


// ---------------------------------------------------------
void BlockBoundary::define(const IntVect& a_permutation,
                           const IntVect& a_sign,
                           const IntVect& a_translation,
                           int a_neighbor,
                           int a_face)
{
  m_bc=-1;
  m_permutation=a_permutation;
  m_sign= a_sign;
  m_translation= a_translation;
  m_type= CONFORMAL;
  m_neighbor = a_neighbor;
  m_face = a_face;
}


// ---------------------------------------------------------
void BlockBoundary::defineConformal(int a_neighbor,
                                    int a_face)
{
  IntVect perm(D_DECL(0, 1, 2));
  // sign is all 1; translation is all 0
  define(perm, IntVect::Unit, IntVect::Zero, a_neighbor, a_face);
}


// ---------------------------------------------------------
/// non-conformal block mating constructor
BlockBoundary::BlockBoundary(int a_neighbor,
                             int a_face)
{
  define(a_neighbor, a_face);
}


// ---------------------------------------------------------
void BlockBoundary::define(int a_neighbor,
                           int a_face)
{
  m_bc=-1;
  m_type=MAPPED;
  m_neighbor = a_neighbor;
  m_face = a_face;
}


// ---------------------------------------------------------
IntVect BlockBoundary::convertOldToNew(const IntVect&   a_ivOld) const
{
  IntVect ivNew;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivNew[idir] =
        m_sign[idir]*a_ivOld[m_permutation[idir]] + m_translation[idir];
    }
  return ivNew;
}


// ---------------------------------------------------------
IntVect BlockBoundary::convertNewToOld(const IntVect&   a_ivNew) const
{
  IntVect ivOld;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivOld[m_permutation[idir]] = m_sign[idir] *
        (a_ivNew[idir] - m_translation[idir]);
    }
  return ivOld;
}


// ---------------------------------------------------------
Box BlockBoundary::convertOldToNew(const Box&   a_bxOld) const
{
  // Box bxOldNodes = surroundingNodes(a_bxOld);
  // IntVect ivNewLo = convertOldToNew(bxOldNodes.smallEnd());
  // IntVect ivNewHi = convertOldToNew(bxOldNodes.bigEnd());
  IntVect ivNewLo = convertOldToNew(a_bxOld.smallEnd());
  IntVect ivNewHi = convertOldToNew(a_bxOld.bigEnd());
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int iLo = Min(ivNewLo[idir], ivNewHi[idir]);
      int iHi = Max(ivNewLo[idir], ivNewHi[idir]);
      ivNewLo[idir] = iLo;
      ivNewHi[idir] = iHi;
    }
  //  Box bxNewNodes(ivNewLo, ivNewHi, IndexType::TheNodeType());
  //  Box bxNewCells = enclosedCells(bxNewNodes);
  Box bxNewCells = Box(ivNewLo, ivNewHi);
  return bxNewCells;
}


// ---------------------------------------------------------
Box BlockBoundary::convertNewToOld(const Box&   a_bxNew) const
{
  //  Box bxNewNodes = surroundingNodes(a_bxNew);
  //  IntVect ivOldLo = convertNewToOld(bxNewNodes.smallEnd());
  //  IntVect ivOldHi = convertNewToOld(bxNewNodes.bigEnd());
  IntVect ivOldLo = convertNewToOld(a_bxNew.smallEnd());
  IntVect ivOldHi = convertNewToOld(a_bxNew.bigEnd());
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int iLo = Min(ivOldLo[idir], ivOldHi[idir]);
      int iHi = Max(ivOldLo[idir], ivOldHi[idir]);
      ivOldLo[idir] = iLo;
      ivOldHi[idir] = iHi;
    }
  //  Box bxOldNodes(ivOldLo, ivOldHi, IndexType::TheNodeType());
  //  Box bxOldCells = enclosedCells(bxOldNodes);
  Box bxOldCells = Box(ivOldLo, ivOldHi);
  return bxOldCells;
}


// ---------------------------------------------------------
BlockBoundary refine(const BlockBoundary&   a_bb,
                     int                    a_refinement_ratio)
{
  //   IntVect translationRefined = a_bb.m_translation * a_refinement_ratio;
  IntVect translationRefined;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_bb.m_sign[idir] == 1)
        translationRefined[idir] =
          a_bb.m_translation[idir]*a_refinement_ratio;
      else // (a_bb.m_sign[idir] == -1)
        translationRefined[idir] =
          (a_bb.m_translation[idir] + 1)*a_refinement_ratio - 1;
    }
  //  BlockBoundary* bbRefined = new BlockBoundary(a_bb.m_permutation,
  //                                               a_bb.m_sign,
  //                                               translationRefined,
  //                                               a_bb.m_neighbor,
  //                                               a_bb.m_face);
  //  bbRefined->m_type = a_bb.m_type;
  //  return *bbRefined;
  BlockBoundary bbRefined(a_bb.m_permutation,
                          a_bb.m_sign,
                          translationRefined,
                          a_bb.m_neighbor,
                          a_bb.m_face);
  bbRefined.m_type = a_bb.m_type;
  return bbRefined;
}

#include "NamespaceFooter.H"
