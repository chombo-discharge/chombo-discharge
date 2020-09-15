#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MINIIVFABI_H_
#define _MINIIVFABI_H_
#include "MayDay.H"
#include "IntVectSet.H"
#include "VoFIterator.H"
#include "parstream.H"
#include "SPMD.H"
#include "DebugOut.H"
#include "EBDebugOut.H"
#include "SPMD.H"
#include "NamespaceHeader.H"


/******************/
template <class T> inline
MiniIVFAB<T>::MiniIVFAB():m_nComp(0),m_data(0)
{
 
}
/******************/
template <class T> inline
MiniIVFAB<T>::~MiniIVFAB()
{
 
}

template <class T> inline
MiniIVFAB<T>::MiniIVFAB(const Interval& a_comps,
                        MiniIVFAB<T>&   a_original)
  :m_vofs(a_original.m_vofs),
   m_Memory(a_original.m_Memory),
   m_nComp(a_comps.size()),
   m_data(a_original.dataPtr(a_comps.begin())){ ; }

/******************/
template <class T> inline
MiniIVFAB<T>::MiniIVFAB(const IntVectSet& a_ivsin,
                        const EBGraph&    a_ebgraph,
                        const int&        a_nvarin)
{
  define(a_ivsin, a_ebgraph, a_nvarin);
}
/******************/
template <class T> inline
void
MiniIVFAB<T>::define(const IntVectSet& a_ivsin,
                     const EBGraph&    a_ebGraph,
                     const int&        a_nvarin)
{
 
  CH_assert(a_nvarin > 0);

  m_nComp=a_nvarin;

  if (!a_ivsin.isEmpty())
    {
      int nVoFs = 0;
      //figure out how long vector has to be
      IVSIterator ivsit(a_ivsin);
      for (ivsit.reset(); ivsit.ok(); ++ivsit)
        {
          nVoFs += a_ebGraph.numVoFs(ivsit());
        }
      //now allocate the vector set the fab to go into
      //the pool of data at the first component of the first
      //vof
      if (nVoFs > 0)
        {
          if(!m_vofs)
            {
              m_vofs=shared_ptr<Vector<VolIndex> >(new Vector<VolIndex>(nVoFs));
              m_Memory= shared_ptr<Vector<T> >(new Vector<T>(nVoFs*a_nvarin));
              m_data = &(m_Memory->operator[](0));
              
              VoFIterator vofit(a_ivsin, a_ebGraph);
              *m_vofs = vofit.getVector();
            }
        }
      else
        {
          m_vofs.reset();
        }
    }
}

template <class T> inline
void
MiniIVFAB<T>::copy(const Box& a_fromBox,
                   const Interval& a_dstInterval,
                   const Box& a_toBox,
                   const MiniIVFAB<T>& a_src,
                   const Interval& a_srcInterval)
{
  T* to=NULL;
  if(m_vofs)
  {
    to = getIndex((*m_vofs)[0],a_dstInterval.begin());
  }
  for(unsigned int i=0; i<m_vofs->size(); i++, to++)
    {
      const VolIndex& vof = (*m_vofs)[i];
      if(a_fromBox.contains(vof.gridIndex()))
        {
          const T* from= a_src.getIndex(vof,a_srcInterval.begin());
          //T* to  = getIndex(vof,a_dstInterval.begin());
          for(unsigned int c=0; c<a_srcInterval.size(); c++)
            {
              to[c*m_vofs->size()]=from[c*a_src.m_vofs->size()];
            }
        }
    }
}

  
  ///
/******************/
template <class T> inline
int MiniIVFAB<T>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  int count = 0;
  T tmp;
  //create set of cells in fab that are also in the input region
  if(m_vofs){
    for (unsigned int i=0; i<m_vofs->size(); i++)
    {
      if (a_region.contains((*m_vofs)[i].gridIndex())) count++;
    }
  }
  if (count > 0)
    {
      return sizeof(int) + count*CH_XD::linearSize((*m_vofs)[0]) + count*a_comps.size()*CH_XD::linearSize(tmp);
    }
  return sizeof(int);
}
/********************/
template <class T> inline
void MiniIVFAB<T>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  char* buffer = (char*)a_buf;
  buffer += sizeof(int);
  int count = 0;
  if (m_vofs)
    {
      const T* ptr = &(this->m_data[0]);
      for (unsigned int i=0; i<m_vofs->size(); i++, ptr++)
        {
          const VolIndex& v = (*m_vofs)[i];
          if (a_region.contains(v.gridIndex()))
            {
              count++;
              CH_XD::linearOut(buffer, v);
              buffer+= CH_XD::linearSize(v);
              for (int c=a_comps.begin(); c<=a_comps.end(); c++)
                {
                  CH_XD::linearOut(buffer, *(ptr+c*(m_vofs->size())));
                  buffer += CH_XD::linearSize(*ptr);
                }
            }
        }
    }
  int* b = (int*)a_buf;
  *b = count;
}
/********************/
template <class T> inline
void MiniIVFAB<T>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  int* b = (int*)a_buf;
  int count = *b;
  char* buffer = (char*)a_buf;
  buffer += sizeof(int);
  for (int i=0; i<count; i++)
    {
      VolIndex v;
      CH_XD::linearIn(v, buffer);
      buffer += linearSize(v);
      for (int c=a_comps.begin(); c<=a_comps.end(); c++)
        {
          T* ptr = getIndex(v, c);
          CH_XD::linearIn(*ptr, buffer);
          buffer+=CH_XD::linearSize(*ptr);
        }
    }
}
/********************/
template <class T> inline
T*
MiniIVFAB<T>::getIndex(const VolIndex& a_vof, const int& a_comp) const
{

  CH_assert((a_comp >= 0) && (a_comp < this->m_nComp));

  T* dataPtr =  (T*)&(this->m_data[0]);
  for (unsigned int i=0; i<m_vofs->size(); ++i)
    {
      if (a_vof == (*m_vofs)[i]) break;
      dataPtr++;
    }
  dataPtr += a_comp*m_vofs->size();
  return dataPtr;
}

/********************/
template <class T> inline
void MiniIVFAB<T>::setVal(const T& a_val)
{
  if(m_vofs)
    for(int i=0; i<m_vofs->size()*m_nComp; i++) m_data[i]=a_val;
}
/********************/
template <class T> inline
void MiniIVFAB<T>::setVal(int a_comp, const T& a_val)
{
  if(m_vofs)
    for(int i=a_comp*m_nComp; i<m_vofs->size()*(a_comp+1)*m_nComp; i++) m_data[i]=a_val;
}

/********************/
template <class T> inline
void MiniIVFAB<T>::setVal(const T& a_val, const Box& a_region, int a_startcomp, int a_ncomp)
{
  if(m_vofs)
    forall(*this, a_region, a_startcomp, a_startcomp, a_ncomp, false, [a_val](T& d, const T& s){d=a_val;});
}



/********************/
template <class T> inline
int MiniIVFAB<T>::numVoFs() const
{
  if(m_vofs)
    return m_vofs->size();
  return 0;
}
/********************/
template <class T> inline
Vector<VolIndex>
MiniIVFAB<T>::getVoFs() const
{
  //if this fails, nothing to get
  if(m_vofs)
    {
      return *m_vofs;
    }
  else
    {
      return Vector<VolIndex>();
    }
}
/***********/
template <class T>
template <typename F>
void MiniIVFAB<T>::forall(const MiniIVFAB<T>& a_src, const Box& a_box, int a_srccomp, int a_destcomp, 
                          int a_numcomp, bool sameRegBox, const F& func)
{
  if(!m_vofs) return;
  if (sameRegBox)
    {
      Real* l = dataPtr(a_destcomp);
      const Real* r = a_src.dataPtr(a_srccomp);
      int nvof = m_vofs->size();
      for (int i=0; i<a_numcomp*nvof; i++)
        {
          func(l[i], r[i]);
        }
    }
  else
    {
      const Vector<VolIndex>& vofs = *m_vofs;
      for(int i=0; i<vofs.size(); i++)
        {
          const VolIndex& vof=vofs[i];
          if(a_box.contains(vof.gridIndex()))
            {
              for (int icomp = 0; icomp < a_numcomp; ++icomp)
                {
                  func(this->operator()(vof, a_destcomp+icomp), a_src(vof, a_srccomp+icomp));
                }
            }
        }
    }
}


#include "NamespaceFooter.H"
#endif

