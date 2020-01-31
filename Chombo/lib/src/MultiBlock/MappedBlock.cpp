#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MappedBlock.H"
#include "UGIO.H"
#include "CH_HDF5.H"

#include "NamespaceHeader.H"


//====================================================================

MappedBlock::MappedBlock()
  :m_map(NULL)
{

}


// ---------------------------------------------------------
MappedBlock::MappedBlock(const Box& a_domain,
                         const BlockMap* a_map,
                         const Tuple<BlockBoundary, 2*SpaceDim>& a_boundaries)
  //  : m_domain(a_domain), m_map(a_map), m_boundaries(a_boundaries)
{
  define(a_domain, a_map, a_boundaries);
  //  CH_assert(a_map != NULL);
  //  CH_assert(!a_domain.isEmpty());
  // need to code in a hypercube completeness check here.
  // a n-dimensional hypercube has 2*n facets of dimension n-1 (boundaries)
}

// ---------------------------------------------------------
void MappedBlock::define(const Box& a_domain,
                         const BlockMap* a_map,
                         const Tuple<BlockBoundary, 2*SpaceDim>& a_boundaries)
{
  CH_assert(a_map != NULL);
  CH_assert(!a_domain.isEmpty());
  m_domain = a_domain;
  m_map = a_map;
  m_boundaries = a_boundaries;
  // need to code in a hypercube completeness check here.
  // a n-dimensional hypercube has 2*n facets of dimension n-1 (boundaries)
}


// ---------------------------------------------------------
MappedBlock refine(const MappedBlock& a_mb,
                   int   a_refinement_ratio)
{
  Box domainRefined = refine(a_mb.m_domain, a_refinement_ratio);
  Tuple<BlockBoundary, 2*SpaceDim> boundariesRefined;
  for (int faceID = 0; faceID < 2*SpaceDim; faceID++)
    {
      boundariesRefined[faceID] =
        refine(a_mb.m_boundaries[faceID], a_refinement_ratio);
    }
  const BlockMap* mapRefinedPtr = a_mb.m_map;
  MappedBlock mb(domainRefined,
                 mapRefinedPtr,
                 boundariesRefined);
  return mb;
}

#include "NamespaceFooter.H"
