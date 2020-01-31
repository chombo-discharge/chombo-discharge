#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFStencil.H"
#include "MFCellFAB.H"
#include "NamespaceHeader.H"
/**************/
MFStencil::
MFStencil(const Vector<MFStencil::agg_t>& a_stencil,
          const Box&                      a_grid,
          const Vector<EBISBox>&          a_vectEBISBox,
          const IntVect&                  a_ghostVectLph,
          const IntVect&                  a_ghostVectPhi,
          int                             a_varDest)
  : m_grid(a_grid),
    m_ghostVectPhi( a_ghostVectPhi ),
    m_ghostVectLph( a_ghostVectLph ),
    m_destVar(a_varDest)
{
  m_ebisBox[0]=(a_vectEBISBox[0]);
  m_ebisBox[1]=(a_vectEBISBox[1]);
  computeOffsets(a_stencil);
}
/**************/
void
MFStencil::
computeOffsets(const Vector<MFStencil::agg_t>& a_stencil)
{
  //debugging hook
  //m_srcVoFs = a_vofStencil;

  Box boxPhi = grow(m_grid, m_ghostVectPhi);
  Box boxLph = grow(m_grid, m_ghostVectLph);
  m_lphBox = boxLph;
  m_phiBox = boxPhi;

  //2 is the number of fluids
  IntVectSet ivsPhi[2];
  IntVectSet ivsLph[2];
  BaseIVFAB<Real> baseivfabPhi[2];
  BaseIVFAB<Real> baseivfabLph[2];
  for (int ifluid = 0; ifluid < 2; ifluid++)
    {
      ivsPhi[ifluid] = m_ebisBox[ifluid].getMultiCells(boxPhi);
      ivsLph[ifluid] = m_ebisBox[ifluid].getMultiCells(boxLph);
      EBGraph ebgraph= m_ebisBox[ifluid].getEBGraph();

      baseivfabPhi[ifluid].define(ivsPhi[ifluid], ebgraph, m_destVar+1);
      baseivfabLph[ifluid].define(ivsLph[ifluid], ebgraph, m_destVar+1);
    }

  const IntVect& smallendPhi = boxPhi.smallEnd();
  const IntVect& smallendLph = boxLph.smallEnd();
  IntVect ncellsPhi = boxPhi.size();
  IntVect ncellsLph = boxLph.size();

  m_ebstencil.resize(a_stencil.size());
  m_destTerms.resize(a_stencil.size());
  m_cache.resize(    a_stencil.size());

  for (int isrc = 0; isrc < a_stencil.size(); isrc++)
    {
      const VolIndex& dstvof  = a_stencil[isrc].destVoF;
      const int&      dfluid  = a_stencil[isrc].destFluid;

      m_destTerms[isrc].fluidid = dfluid;
      if (m_ebisBox[dfluid].numVoFs(dstvof.gridIndex()) > 1)
        {//multi-valued (the dataPtr(0) is correct--that is where we start from)
          m_destTerms[isrc].offset =
            baseivfabLph[dfluid].getIndex(dstvof, m_destVar) -
            baseivfabLph[dfluid].dataPtr(0);
          m_destTerms[isrc].multiValued = true;
        }
      else
        {//single-valued
          IntVect ivLph = dstvof.gridIndex()  - smallendLph;
          IntVect ivPhi = dstvof.gridIndex()  - smallendPhi;
          m_destTerms[isrc].offset = ivLph[0] + ivLph[1]*ncellsLph[0] ;
#if CH_SPACEDIM==3
          m_destTerms[isrc].offset +=  ivLph[2]*ncellsLph[0]*ncellsLph[1];
#endif

          //add in term due to variable number
#if CH_SPACEDIM==2
          m_destTerms[isrc].offset += m_destVar*ncellsLph[0]*ncellsLph[1];
#elif CH_SPACEDIM==3
          m_destTerms[isrc].offset += m_destVar*ncellsLph[0]*ncellsLph[1]*ncellsLph[2];
#else
          bogus_spacedim();
#endif
          m_destTerms[isrc].multiValued = false;
        }

      for (int ifluid = 0; ifluid < 2; ifluid++)
        {
          const VoFStencil& sten = a_stencil[isrc].stenFluid[ifluid];

          for (int isten = 0; isten < sten.size(); isten++)
            {
              const VolIndex stencilVof = sten.vof(isten);
              int srcVar = sten.variable(isten);
              stencilTerm stenEntry;
              stenEntry.weight   = sten.weight(isten);
              stenEntry.fluidid  = ifluid;
              if (m_ebisBox[ifluid].numVoFs(stencilVof.gridIndex()) > 1)
                {//multi-valued (the dataPtr(0) is correct--that is where we start from)
                  stenEntry.offset =
                    baseivfabPhi[ifluid].getIndex(stencilVof, srcVar) -
                    baseivfabPhi[ifluid].dataPtr(0);
                  //single vs multi encapsulated in vector length
                  m_ebstencil[isrc].multi.push_back(stenEntry);
                }
              else
                {//single-valued
                  IntVect ivPhi = stencilVof.gridIndex()  - smallendPhi;
                  stenEntry.offset = ivPhi[0] + ivPhi[1]*ncellsPhi[0] ;
#if CH_SPACEDIM==3
                  stenEntry.offset +=  ivPhi[2]*ncellsPhi[0]*ncellsPhi[1];
#endif
                  //add in term due to variable number
#if CH_SPACEDIM==2
                  stenEntry.offset += srcVar*ncellsPhi[0]*ncellsPhi[1];
#elif CH_SPACEDIM==3
                  stenEntry.offset += srcVar*ncellsPhi[0]*ncellsPhi[1]*ncellsPhi[2];
#else
                  bogus_spacedim();
#endif
                  //single vs multi encapsulated in vector length
                  m_ebstencil[isrc].single.push_back(stenEntry);
                }
            }
        }
    }
  //  CH_STOP(t1);
}
/**************/
void
MFStencil::
apply(MFCellFAB& a_lph, const MFCellFAB& a_phi, bool a_incrementOnly) const
{
  CH_TIME("MFStencil::apply");

  CH_assert(a_lph.getPhase(0).getSingleValuedFAB().box() == m_lphBox);
  CH_assert(a_phi.getPhase(0).getSingleValuedFAB().box()    == m_phiBox);

  const Real* singleValuedPtrPhi[2];
  Real*       singleValuedPtrLph[2];
  const Real*  multiValuedPtrPhi[2];
  Real*        multiValuedPtrLph[2];
  for (int ifluid = 0; ifluid < 2; ifluid++)
    {
      singleValuedPtrPhi[ifluid] = a_phi.getPhase(ifluid).getSingleValuedFAB().dataPtr(0);
      singleValuedPtrLph[ifluid] = a_lph.getPhase(ifluid).getSingleValuedFAB().dataPtr(0);
      multiValuedPtrPhi[ifluid] =  a_phi.getPhase(ifluid).getMultiValuedFAB().dataPtr( 0);
      multiValuedPtrLph[ifluid] =  a_lph.getPhase(ifluid).getMultiValuedFAB().dataPtr( 0);
    }

  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      const ebstencil_t& ebstencil = m_ebstencil[isrc];
      Real* lphiPtr = NULL;
      const int& dstFlu = m_destTerms[isrc].fluidid;
      if (m_destTerms[isrc].multiValued)
        {
          lphiPtr  = multiValuedPtrLph[dstFlu] + m_destTerms[isrc].offset;
        }
      else
        {
          lphiPtr  = singleValuedPtrLph[dstFlu] + m_destTerms[isrc].offset;
        }

      Real& lphi = *lphiPtr;
      if (!a_incrementOnly)
        {
          lphi =  0.;
        }
      //single-valued
      for (int isingle = 0; isingle < ebstencil.single.size(); isingle++)
        {
          const int & offset = ebstencil.single[isingle].offset;
          const int & srcFlu = ebstencil.single[isingle].fluidid;
          const Real& weight = ebstencil.single[isingle].weight;

          const Real& phiVal = *(singleValuedPtrPhi[srcFlu] + offset);
          lphi += phiVal*weight;
        }
      //multi-valued
      for (int imulti = 0; imulti < ebstencil.multi.size(); imulti++)
        {
          const int & offset = ebstencil.multi[imulti].offset;
          const int & srcFlu = ebstencil.multi[imulti].fluidid;
          const Real& weight = ebstencil.multi[imulti].weight;

          const Real& phiVal = *(multiValuedPtrPhi[srcFlu] + offset);
          lphi += phiVal*weight;
        }
    }

}

/**************/
void
MFStencil::
cache(const MFCellFAB& a_lph) const
{
  CH_assert(a_lph.getPhase(0).getSingleValuedFAB().box() == m_lphBox);

  const Real* singleValuedPtrLph[2];
  const Real*  multiValuedPtrLph[2];
  for (int ifluid = 0; ifluid < 2; ifluid++)
    {
      singleValuedPtrLph[ifluid] = a_lph.getPhase(ifluid).getSingleValuedFAB().dataPtr(0);
      multiValuedPtrLph[ifluid] =  a_lph.getPhase(ifluid).getMultiValuedFAB().dataPtr( 0);
    }
  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      const int& dstFlu = m_destTerms[isrc].fluidid;
      if (m_destTerms[isrc].multiValued)
        {
          m_cache[isrc] = *(multiValuedPtrLph[dstFlu] + m_destTerms[isrc].offset);
        }
      else
        {
          m_cache[isrc] = *(singleValuedPtrLph[dstFlu] + m_destTerms[isrc].offset);
        }
    }
}
/**************/
void
MFStencil::
uncache(MFCellFAB& a_lph) const
{
  CH_assert(a_lph.getPhase(0).getSingleValuedFAB().box() == m_lphBox);

  Real* singleValuedPtrLph[2];
  Real*  multiValuedPtrLph[2];
  for (int ifluid = 0; ifluid < 2; ifluid++)
    {
      singleValuedPtrLph[ifluid] = a_lph.getPhase(ifluid).getSingleValuedFAB().dataPtr(0);
      multiValuedPtrLph[ifluid] =  a_lph.getPhase(ifluid).getMultiValuedFAB().dataPtr( 0);
    }
  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      const int& dstFlu = m_destTerms[isrc].fluidid;
      if (m_destTerms[isrc].multiValued)
        {
          *(multiValuedPtrLph[dstFlu] + m_destTerms[isrc].offset) = m_cache[isrc];
        }
      else
        {
          *(singleValuedPtrLph[dstFlu] + m_destTerms[isrc].offset) = m_cache[isrc];
        }
    }
}
/**************/

#include "NamespaceFooter.H"
