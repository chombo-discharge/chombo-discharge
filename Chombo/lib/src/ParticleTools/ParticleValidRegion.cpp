#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#include "ParticleValidRegion.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

ParticleValidRegion::ParticleValidRegion(const DisjointBoxLayout&         a_grids,
                                         const LevelData<BaseFab<bool> >* a_crseValid,
                                         const int                        a_refRatio,
                                         const int                        a_buffer)
  :
  m_isDefined(false)
{
  define(a_grids,a_crseValid,a_refRatio,a_buffer);
}

void ParticleValidRegion::define(const DisjointBoxLayout&         a_grids,
                                 const LevelData<BaseFab<bool> >* a_crseValid,
                                 const int                        a_refRatio,
                                 const int                        a_buffer)
{

  CH_TIME("ParticleValidRegion::define");

  m_isDefined = true;

  // one component, no ghosts
  m_valid.define(a_grids, 1);
  DataIterator di(a_grids);

  // temporary valid with ghost cells
  LevelData<BaseFab<bool> > tempValid(a_grids, 1, IntVect::Unit);

  // initialize tempValid so that only ghost cells are invalid (=false)
  for (di.begin(); di.ok(); ++di)
    {
      tempValid[di].setVal(false);
      tempValid[di].setVal(true, a_grids.get(di), 0);
    }

  // make valid aware of requirements on coarser level
  if (a_crseValid != NULL)
    {
      // create fine coarsened grids
      DisjointBoxLayout crseDBL;
      coarsen (crseDBL, a_grids, a_refRatio);
      const DisjointBoxLayout coarsenedGrids = crseDBL;

      LevelData<BaseFab<bool> > crsenedValid(coarsenedGrids, 1);
      a_crseValid->copyTo(Interval(0,0), crsenedValid, Interval(0,0));

      for (di.begin(); di.ok(); ++di)
        {
          for (BoxIterator bi(a_grids.get(di)); bi.ok(); ++bi)
            {
              const IntVect civ  = coarsen(bi(), a_refRatio);
              tempValid[di](bi(),0) = crsenedValid[di](civ, 0);
            }
        }
    }

  // if buffer has zero size, only coarser buffer cells will be counted
  if (a_buffer > 0)
    {
      // fill in ghost cells
      tempValid.exchange();

      // so copy the values of these
      for (di.begin(); di.ok(); ++di)
        {
          const Box& box = a_grids.get(di);
          for (int dir=0; dir<SpaceDim; ++dir)
            {
              bufferGhostJob(tempValid[di], box, box, a_buffer, dir);
            }
        }
    }

  // define a_valid
  for (di.begin(); di.ok(); ++di)
    {
      const Box& box = a_grids.get(di);
      m_valid[di].copy(tempValid[di],box,0,box,0,1);
    }
}

///
void ParticleValidRegion::bufferGhostJob(BaseFab<bool>& a_valid,
                                         const Box&     a_box,
                                         const Box&     a_domain,
                                         const int      a_buffer,
                                         const int      a_dir)
{

  CH_TIME("ParticleValidRegion::bufferGhostJob");

  // First construct them as "adjacent" boxes one cell deep
  const Box loAdjBox = adjCellLo(a_box,a_dir,1);
  const Box hiAdjBox = adjCellHi(a_box,a_dir,1);

  for (int nxt_dir=a_dir+1; nxt_dir<SpaceDim; ++nxt_dir)
    {
      bufferGhostJob(a_valid,loAdjBox,a_domain,a_buffer,nxt_dir);
      bufferGhostJob(a_valid,hiAdjBox,a_domain,a_buffer,nxt_dir);
    }

  // buffer region around boundary cells
  const Box buffer(-a_buffer*IntVect::Unit, a_buffer*IntVect::Unit);

  Box invalidLo;
  for (BoxIterator abi(loAdjBox); abi.ok(); ++abi)
    {
      if (!a_valid(abi(), 0))
        {
          const Box bbox = (buffer + abi());
	  invalidLo.minBox(bbox);
	}
    }

  const Box bboxlo = invalidLo & a_domain;
  a_valid.setVal(false, bboxlo, 0);

  Box invalidHi;
  for (BoxIterator abi(hiAdjBox); abi.ok(); ++abi)
    {
      if (!a_valid(abi(),0))
        {
          const Box bbox = (buffer+abi());
	  invalidHi.minBox(bbox);
	}
    }

  const Box bboxhi = invalidHi & a_domain;
  a_valid.setVal(false, bboxhi, 0);
    
}

#include "NamespaceFooter.H"
