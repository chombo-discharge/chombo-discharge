#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBAggStencil.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"
///
void
MBAggStencil::
computeSlabSizes()
{
  bool m_reducedDimStencil = (m_stencilDim < SpaceDim);
  m_srcSlabLinearSize = -1;
  m_dstSlabLinearSize = -1;
  m_numSlabs = 1;
  if(m_reducedDimStencil)
    {
      m_srcSlabLinearSize = 1;
      m_dstSlabLinearSize = 1;
      for(int idim = 0; idim < m_stencilDim; idim++)
        {
          m_srcSlabLinearSize *= m_srcBoxSize[idim];
          m_dstSlabLinearSize *= m_dstBoxSize[idim];
        } 
      //if reduced dimension, we just step along the next higher dimensions
      for(int idim = m_stencilDim; idim < SpaceDim; idim++)
        {
          m_numSlabs *= m_dstBoxSize[idim];
        }
    }
}
///
size_t 
MBAggStencil::
getOffset(const IntVect& a_iv, const Box& a_grid)
{
  size_t retval = 0;
  IntVect gsize = a_grid.size();
  IntVect ivdiff = a_iv - a_grid.smallEnd();
  retval = ivdiff[0];
  //I could be fancier with some kind of nested loop here.
  //This seems clearer.
  if(m_stencilDim > 1)
    {
      retval += gsize[0]*ivdiff[1];
    }
  if(m_stencilDim > 2)
    {
      retval += gsize[0]*gsize[1]*ivdiff[2];
    }
  CH_assert(m_stencilDim <= 3);

  return retval;
}
///
void
MBAggStencil::
computeOffsets(const RefCountedPtr<IVSFAB<MBStencil> >  & a_stencil,
               const IntVectSet                         & a_region,
               const Box                                & a_srcBox,
               const Box                                & a_dstBox)
{
  const IntVectSet& ivs = a_region;
  m_dstAccess.resize(ivs.numPts());
  m_mbstencil.resize(ivs.numPts());
  
  int icell = 0;
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit, ++icell)
    {
      const IntVect& ivDst = ivsit();

      size_t offset = getOffset(ivDst, a_dstBox);
      stencil_t ptsten;

      const MBStencil& sten = (*a_stencil)(ivDst, 0);
      ptsten.resize(sten.size());
      for (int isten = 0; isten < sten.size(); isten++)
        {
          const MBStencilElement& elem = sten[isten];
          ptsten[isten].m_offset  = getOffset(elem.cell(), a_srcBox);
          ptsten[isten].m_weight  = elem.weight();
          ptsten[isten].m_blockid = elem.block();
        }
      m_dstAccess[icell] = offset;
      m_mbstencil[icell] = ptsten;
    }
}
///
MBAggStencil::
MBAggStencil(const RefCountedPtr<IVSFAB<MBStencil> >  & a_stencil,
             const IntVectSet                         & a_region,
             const Box                                & a_srcBox,
             const Box                                & a_dstBox,
             const int                                & a_stencilDim)
{
  CH_TIME("MMBAggSten.constructor");
  CH_assert(a_stencilDim <= SpaceDim);
  
  m_stencilDim = a_stencilDim;
  m_srcBoxSize = a_srcBox.size();
  m_dstBoxSize = a_dstBox.size();

  computeSlabSizes();

  computeOffsets(a_stencil, a_region, a_srcBox, a_dstBox);
}
///
void
MBAggStencil::
apply(FArrayBox                 & a_lph,
      const Vector<FArrayBox*>  & a_phi,
      const int                 & a_src,
      const int                 & a_dst,
      const int                 & a_nco,
      const bool                & a_incrementOnly) const
{
  CH_TIME("MMBAggSten::apply");

  for (int icomp = 0; icomp < a_nco; icomp++)
    {
      
      int varDst = a_dst + icomp;
      int varSrc = a_src + icomp;
      Real*               dataPtrLph = a_lph.dataPtr(varDst);
      Vector<const Real*> dataPtrsPhi(a_phi.size());

      for(int iblock = 0; iblock < a_phi.size(); iblock++)
        {
          dataPtrsPhi[iblock] = a_phi[iblock]->dataPtr(varSrc);
        }

      for(int islab = 0; islab < m_numSlabs; islab++)
        {
          for (int idst = 0; idst < m_mbstencil.size(); idst++)
            {
              Real* lphiPtr =  dataPtrLph + m_dstAccess[idst];

              Real& lphi = *lphiPtr;
              if (!a_incrementOnly)
                {
                  lphi =  0.;
                }

              const stencil_t& ebstencil = m_mbstencil[idst];
              for (int isten = 0; isten < ebstencil.size(); isten++)
                {
                  const stenelem_t & elem = ebstencil[isten];

                  const Real   & weight  = elem.m_weight;
                  const size_t & offset  = elem.m_offset;
                  const int    & blockid = elem.m_blockid;

                  const Real& phiVal = *(dataPtrsPhi[blockid] + offset);
                  lphi += phiVal*weight;
                }
            }
          //increment pointers to next slab
          dataPtrLph += m_dstSlabLinearSize;
          for(int iblock = 0; iblock < a_phi.size(); iblock++)
            {
              dataPtrsPhi[iblock] += m_srcSlabLinearSize;
            }
          
        } // end loop over slabs
    } // end of loop over components  
}
///

#include "NamespaceFooter.H"
