#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFRemapper.H"
#include "MFAliasFactory.H"
#include "EBPWLFineInterp.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

IntVect direction(const RealVect& normal)
{
  IntVect rtn = IntVect::Zero;
  int maxIndex = 0;
  for (int i=1; i<CH_SPACEDIM; ++i)
    {
      if (abs(normal[i]) > abs(normal[maxIndex])) maxIndex=i;
    }
  rtn[maxIndex] = (normal[maxIndex] > 0.0 ? 1:-1);
  return rtn;
}

MFRemapper::MFRemapper()
{
}

MFRemapper::~MFRemapper()
{
}

void MFRemapper::remap(const MFIndexSpace& a_sourceMF,
                       const LevelData<MFCellFAB>& a_source,
                       const MFIndexSpace& a_destMF,
                       LevelData<MFCellFAB>& a_dest) const
{
  int numPhases =  a_sourceMF.numPhases();
  for (int i=0; i<numPhases ; i++)
  {
    LevelData<EBCellFAB> source, dest;
    aliasMF(source, i, a_source);
    aliasMF(dest, i, a_dest);
    int ncomp = source.nComp();
    DataIterator dit = source.dataIterator();
    DisjointBoxLayout dbl = source.disjointBoxLayout();

    for (dit.begin(); dit.ok(); ++dit)
      {
        const EBCellFAB& s = source[dit];
        EBCellFAB&       d = dest[dit];
        {
          // first all single-valued cells that have remained single-valued
          // get filled in with a regular copy.
          const FArrayBox& sfab = s.getFArrayBox();
          FArrayBox&       dfab = d.getFArrayBox();
          dfab.copy(sfab);
        }
        const EBISBox& sbox = s.getEBISBox();
        const EBISBox& dbox = d.getEBISBox();
        const ProblemDomain& domain = sbox.getDomain();

        Box region = dbl.get(dit());
        //region.grow(1);// need to pass in nghost to do this properly

        {
          // now, cells that were multi-valued and remain multi-valued
          // are copied.  this step might break is multi-valued in one
          // phase get really complicated.  bug-in-waiting perhaps.
          IntVectSet source_mcells = sbox.getMultiCells(region);
          IntVectSet mcells   = dbox.getMultiCells(region);
          mcells &= source_mcells;
          IVSIterator it(mcells);
          for (;it.ok(); ++it)
            {
              Vector<VolIndex> mvols = sbox.getVoFs(it());
              for (int m=0; m<mvols.size(); m++)
              {
                for (int i=0; i<ncomp; i++)
                {
                  d(mvols[m], i) = s(mvols[m], i);
                }
              }
            }
        }
        {
          //cells that WERE multi-valued, but are now single valued
          //or covered
          IntVectSet mcells = sbox.getMultiCells(region);
          IntVectSet dmcells   = dbox.getMultiCells(region);
          mcells -= dmcells;
          IVSIterator it(mcells);
          for (;it.ok(); ++it)
            {
              if (!dbox.isCovered(it()))
              { //don't bother setting covered
                VolIndex single = dbox.getVoFs(it())[0];
                Vector<VolIndex> mvols = sbox.getVoFs(it());
                Real oldVolume = 0.0;
                for (int m=0; m<mvols.size(); m++)
                {
                  oldVolume += sbox.volFrac(mvols[m]);
                }
                for (int i=0; i<ncomp; i++)
                {
                  Real& val = d(single, i);
                  val = 0.0;
                  for (int m=0; m<mvols.size(); m++)
                  {
                    val += s(mvols[m], i)*sbox.volFrac(mvols[m]);
                  }
                  val /= oldVolume; // need to check for very small old
                }
              }
            }
        }
        {
          // cells that WERE single-valued and are now multi-valued
          IntVectSet mcells   = dbox.getMultiCells(region);
          IntVectSet smcells  = sbox.getMultiCells(region);
          mcells -= smcells; //new multivalued, previously single or covered
          IVSIterator it(mcells);
          for (;it.ok(); ++it)
            {
              if (!sbox.isCovered(it()))
              { //exclude previously covered
                VolIndex single = sbox.getVoFs(it())[0];
                Vector<VolIndex> mvols = dbox.getVoFs(it());
                for (int m=0; m<mvols.size(); m++)
                {
                  for (int i=0; i<ncomp; i++)
                  {
                    d(mvols[m], i) = s(single, i);
                  }
                }
              }
            }
        }
        {
          //cells that were covered, and are now single or multi-valued
          // assume that surface can't move far enough in one step
          // to make a covered cell into a regular cell
          IntVectSet cells = dbox.getIrregIVS(region);
          IVSIterator it(cells);
          for (it.begin(); it.ok(); ++it)
            {
              if (sbox.isCovered(it()))
                {
                  Vector<VolIndex> mvols = dbox.getVoFs(it());
                  for (int m=0; m<mvols.size(); m++)
                  {
                    const VolIndex& vi = mvols[m];
                    const RealVect& norm = dbox.normal(vi);
                    IntVect nearest(it());
                    nearest += direction(norm);

                    int bb=0;
                    VolIndex near; near.define(nearest,bb);
                    while ((!domain.contains(nearest) ||
                          !dbox.isConnected(near, vi)) && bb < 30)
                      {
                        bb++;
                        near.define(nearest, bb);

                      }
                    if (bb==30)
                      {
                        // dam, inward normal has steered me wrong. find any neighbour
                        Box close(it(), it());
                        close.grow(1);
                        for (BoxIterator bit(close); bit.ok(); ++bit)
                          {
                            for (int i=0; i<4; ++i)
                            {
                              if ((domain.contains(bit()) && !dbox.isCovered(bit()))
                                 && bit() != it())
                                {
                                  near.define(bit(),i);
                                  if (dbox.isConnected(near, vi))
                                    {
                                      bb=0;
                                      break;
                                    }
                                }
                            }
                          }
                      }
                    if (bb==30)
                    {
                      MayDay::Error("I can't seem to find a near neighbour");
                    }
                    for (int i=0; i<ncomp; i++)
                    {
                      d(vi, i) = d(near, i);
                    }
                  }
                }
            }

        }
      }
    dest.exchange(dest.interval());
  }

}

void MFRemapper::remap(const MFIndexSpace& a_MF,
                       const ProblemDomain& a_domainCoar,
                       const ProblemDomain& a_domainFine,
                       const LevelData<MFCellFAB>& a_source,
                       const LevelData<MFCellFAB>& a_coarse,
                       const int& nref,
                       const int& nghost,
                       LevelData<MFCellFAB>&  a_dest)
{
  int numPhases =  a_MF.numPhases();
  Vector<LevelData<EBCellFAB>* > source(numPhases), coarse(numPhases), dest(numPhases);

  for (int i=0; i<numPhases ; i++)
  {
    source[i] = new LevelData<EBCellFAB>;
    coarse[i] = new LevelData<EBCellFAB>;
    dest[i]   = new LevelData<EBCellFAB>;
  }
  aliasMF(source, numPhases, a_source);
  aliasMF(coarse, numPhases, a_coarse);
  aliasMF(dest  , numPhases, a_dest  );
  DisjointBoxLayout fineDBL(a_dest.getBoxes()), coarDBL(a_coarse.getBoxes());

  for (int i=0; i<numPhases; i++)
    {
      EBISLayout fineEBISL, coarEBISL;
      a_MF.fillEBISLayout(fineEBISL, i, fineDBL,  a_domainFine.domainBox(),nghost);
      a_MF.fillEBISLayout(coarEBISL, i, coarDBL,  a_domainCoar.domainBox(),nghost);
      EBPWLFineInterp interpOp;
      interpOp.define(fineDBL,
                       coarDBL,
                       fineEBISL,
                       coarEBISL,
                       a_domainCoar,
                       nref,
                      source[i]->nComp(),
                      a_MF.EBIS(i));

      interpOp.interpolate(*dest[i],
                           *coarse[i],
                           dest[i]->interval());

      source[i]->copyTo(source[i]->interval(), *dest[i], dest[i]->interval());

      delete source[i];  source[i]=NULL;
      delete dest[i] ;   dest[i]=NULL;
      delete coarse[i];  coarse[i] = NULL;
    }
}
#include "NamespaceFooter.H"
