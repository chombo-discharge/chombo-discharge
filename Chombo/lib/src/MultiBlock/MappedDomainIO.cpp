#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MappedDomainIO.H"
#include "FArrayBox.H"
#include "BoxIterator.H"
#include "UGIO.H"
#include "LoadBalance.H"
#include "CH_HDF5.H"

#include "NamespaceHeader.H"


//============================================================================


void  MappedDomainIO::writeUnmappedDomainPlotfile(const MappedDomain& a_domain,
                                                  int a_coarseningFactor) const
{
#ifdef USE_HDF
  int nblocks = a_domain.m_blocks.size();
  CH_assert(nblocks > 0);

  Vector<Box> boxes;
  Box dom;
  for (int iblock = 0; iblock < nblocks; ++iblock)
    {
      Box b = a_domain.m_blocks[iblock].box();
      b.coarsen(a_coarseningFactor);
      boxes.push_back(b);
      dom.minBox(b);
    }
  Vector<int> procs(boxes.size());
  LoadBalance(procs, boxes);
  DisjointBoxLayout layout(boxes, procs);
  LevelData<FArrayBox> dummy(layout, 1);
  for (DataIterator dit = dummy.dataIterator(); dit.ok(); ++dit)
    {
      dummy[dit].setVal(dit().intCode());
    }
  WriteUGHDF5(filename, layout, dummy, dom);
#endif
}

// ---------------------------------------------------------
void MappedDomainIO::writeDomainPlotfile(const MappedDomain& a_domain,
                                         int a_coarseningFactor) const
{
#ifdef USE_HDF
  int nblocks = a_domain.m_blocks.size();
  CH_assert(nblocks > 0);

  Vector<Box> boxes;
  Box dom;
  for (int iblock = 0; iblock < nblocks; ++iblock)
    {
      Box b = a_domain.m_blocks[iblock].box();
      b.coarsen(a_coarseningFactor);
      boxes.push_back(b);
      dom.minBox(b);
    }
  Vector<int> procs(boxes.size());
  LoadBalance(procs, boxes);
  DisjointBoxLayout layout(boxes, procs);
  LevelData<FArrayBox> mapping(layout, CH_SPACEDIM, IntVect::Unit);

  for (DataIterator dit = mapping.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&  data = mapping[dit];
      const Box& box = mapping.box(dit());
      FArrayBox  nodeFAB;
      const MappedBlock* mapblock=NULL;
      for (int iblock = 0; iblock < nblocks; ++iblock)
        {
          Box b = a_domain.m_blocks[iblock].box();
          b.coarsen(a_coarseningFactor);
          if (b == box)
            {
              mapblock = &(a_domain.m_blocks[iblock]);
              iblock = nblocks;
            }
        }
      CH_assert(mapblock!= NULL);
      CH_assert(mapblock->hasMap());

      mapblock->map().getCornerCoordinates(nodeFAB, box);
      nodeFAB.shiftHalf(IntVect::Unit);
      data.copy(nodeFAB);
    }
  WriteUGHDF5(filename, layout, mapping, dom);
#endif
}


// ---------------------------------------------------------
void MappedDomainIO::read(MappedDomain& a_domain) const
{


}


// ---------------------------------------------------------
void writeMap(HDF5Handle& a_handle,
              int a_level,
              const DisjointBoxLayout& a_layout,
              const MappedDomain& a_map)
{


}



#include "NamespaceFooter.H"
