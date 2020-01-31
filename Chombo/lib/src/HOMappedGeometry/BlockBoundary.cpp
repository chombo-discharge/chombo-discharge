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
BlockBoundary::BlockBoundary(const IndicesTransformation& a_transformation,
                             int a_neighbor)
{
  define(a_transformation, a_neighbor);
}


// ---------------------------------------------------------
void BlockBoundary::define(const IndicesTransformation& a_transformation,
                           int a_neighbor)
{
  m_bc=-1;
  m_transformation = a_transformation;
  m_type= MAPPED;
  m_neighbor = a_neighbor;
}


// ---------------------------------------------------------
void BlockBoundary::defineConformal(int a_neighbor)
{
  define(IndicesTransformation::Identity, a_neighbor);
}


// ---------------------------------------------------------
int BlockBoundary::reorientFace(int a_dir) const
{
  // Recall:
  // vecNew[idir] == m_sign[idir] * vecOld[m_permutation[idir]];
  // So if vecOld == BASISV(a_dir):
  // vecNew[idir] == m_sign[idir] iff m_permutation[idir] == a_dir;
  // vecNew[idir] == 0 iff m_permutation[idir] != a_dir.
  // Hence, this function should return m_sign[idir] where
  // m_permutation[idir] == a_dir.
  const IndicesTransformation& transformation = getTransformation();
  const IndicesTransformation& transInverse = transformation.inverse();
  const IntVect& signInverse = transInverse.getSign();
  return signInverse[a_dir];
}


// ---------------------------------------------------------
int BlockBoundary::dirOther(int a_dir) const
{
  // Recall:
  // vecNew[idir] == m_sign[idir] * vecOld[m_permutation[idir]];
  const IndicesTransformation& transformation = getTransformation();
  IntVect perm = transformation.getPermutation();
  return perm[a_dir];
}


// ---------------------------------------------------------
// friend function
BlockBoundary coarsen(const BlockBoundary&   a_bb,
                      int                    a_refinement_ratio)
{
  const IndicesTransformation& transformation = a_bb.getTransformation();

  IndicesTransformation transformationCoarsened =
    transformation.coarsen(a_refinement_ratio);

  BlockBoundary bbRefined(transformationCoarsened,
                          a_bb.m_neighbor);
  bbRefined.m_type = a_bb.m_type;
  return bbRefined;
}

// ---------------------------------------------------------
// friend function
BlockBoundary refine(const BlockBoundary&   a_bb,
                     int                    a_refinement_ratio)
{
  const IndicesTransformation& transformation = a_bb.getTransformation();

  IndicesTransformation transformationRefined =
    transformation.refine(a_refinement_ratio);

  BlockBoundary bbRefined(transformationRefined,
                          a_bb.m_neighbor);
  bbRefined.m_type = a_bb.m_type;
  return bbRefined;
}

#include "NamespaceFooter.H"
