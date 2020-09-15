#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IVSFABI_H_
#define _IVSFABI_H_

#include "MayDay.H"
#include "IntVectSet.H"
#include "parstream.H"
#include "SPMD.H"
#include "DebugOut.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

template <class T>
bool IVSFAB<T>::s_verbose = false;

template <class T>
void
IVSFAB<T>::setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}
/*************/
template <class T> inline
const IntVectSet&
IVSFAB<T>::getIVS() const
{
  return m_ivs;
}
/******************/
template <class T> inline
IVSFAB<T>::IVSFAB():m_dataPtr(0),m_nComp(0),m_nIvs(0),m_loVect(IntVect::Unit), m_hiVect(IntVect::Zero), m_isDefined(false)
{

}
/******************/
template <class T> inline
IVSFAB<T>::~IVSFAB()
{
  clear();
}
/******************/
template <class T> inline
IVSFAB<T>::IVSFAB(const IntVectSet& a_ivsin,
                  const int&        a_nvarin)
{
  define( a_ivsin, a_nvarin );
}


/******************/
template<class T> inline
void IVSFAB<T>::define(const IntVectSet& a_ivsin,
                       const int& a_nvarin )
{ 
  m_dataPtr = NULL;
  m_nComp = a_nvarin;
  m_nIvs = 0;
  m_loVect = IntVect::Zero;
  m_hiVect = IntVect::Unit;
  m_ivmap.define( a_ivsin.minBox(), a_nvarin );
  m_ivs = a_ivsin;

  CH_assert(a_nvarin > 0);
  if (!a_ivsin.isEmpty())
    {
      IVSIterator ivsit(a_ivsin);
      for (ivsit.reset(); ivsit.ok(); ++ivsit)
        {
          const IntVect& iv = ivsit();
          m_ivmap(iv,0) = m_nIvs;
          m_nIvs++;
        }

      if ( (m_nIvs > 0) && m_nComp > 0)
        {
          m_dataPtr = new T[m_nComp*m_nIvs];
        }
    }
  if (m_nIvs > 0)
    {
      //set up face low and high vectors
      for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
          m_loVect[dir] = 1;
          m_hiVect[dir] = 1;
        }
      m_hiVect[0] = m_nIvs;
      
    }
  else
    {
      m_loVect = IntVect::Unit;
      m_hiVect = IntVect::Zero;
      m_dataPtr = NULL;
    }
  m_isDefined=true;

}
/******************/
template <class T> inline
void
IVSFAB<T>::setVal(const T& a_value)
{
  CH_assert(isDefined());
  for (int ivec = 0; ivec < m_nIvs*m_nComp; ivec++)
    m_dataPtr[ivec] = a_value;
}
/******************/
template <class T> inline
void
IVSFAB<T>::copy(const Box& a_fromBox,
                   const Interval& a_dstInterval,
                   const Box& a_toBox,
                   const IVSFAB<T>& a_src,
                   const Interval& a_srcInterval)
{
  IntVect shift = a_toBox.smallEnd() - a_fromBox.smallEnd();
  IntVect shift2 = a_toBox.bigEnd() - a_fromBox.bigEnd();
  CH_assert(shift==shift2);
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  CH_assert(a_srcInterval.size() == a_dstInterval.size());
  CH_assert(a_dstInterval.begin() >= 0);
  CH_assert(a_srcInterval.begin() >= 0);
  CH_assert(a_dstInterval.end()   < m_nComp);
  CH_assert(a_srcInterval.end()   < a_src.m_nComp);

  if ((!m_ivs.isEmpty()) && (!a_src.m_ivs.isEmpty()))
    {
      // my points in toBox
      IntVectSet ivsIntersect = m_ivs & a_toBox;
      // his points in fromBox
      IntVectSet set2 = a_src. m_ivs & a_fromBox;
      // his points mapped to my Box
      set2.shift(shift);
      // overlap
      ivsIntersect &= set2;
      
      IVSFAB<T>& thisFAB = *this;
      int compSize = a_srcInterval.size();
      for (IVSIterator ivit(ivsIntersect); ivit.ok(); ++ivit)
        {
          const IntVect& iv = ivit();
          IntVect iv2 = iv - shift;
          for (int icomp = 0; icomp < compSize; icomp++)
            {
              int isrccomp = a_srcInterval.begin() + icomp;
              int idstcomp = a_dstInterval.begin() + icomp;
              thisFAB(iv, idstcomp) = a_src(iv2, isrccomp);
            }
        }
    }
}
/******************/
//template specializations for real
template <> inline
int IVSFAB<Real>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  // size of actual data plus size of count ..
  // specialization for T=Real.
  IntVectSet subset = m_ivs & a_region;
  int numberThings = subset.numPts();
  int retval = numberThings*sizeof(Real)*a_comps.size() + sizeof(int);
  retval += SpaceDim*sizeof(int)*numberThings; // for intvects
  retval += 2*SpaceDim*sizeof(Real); // box done as reals for some reason
  return retval;
}

/********************/
template <> inline void IVSFAB<Real>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  // message includes count
  IntVectSet subset = m_ivs & a_region;

  // input the box
  IntVect ibl, ibh;
  Real* rptr = (Real*)a_buf;
  for ( int idir=0; idir<SpaceDim; idir++ ) {
    ibl[idir] = *rptr++;
    ibh[idir] = *rptr++;
  }

  // compute the offset (for periodic case)
  IntVect shift = a_region.smallEnd() - ibl;
  IntVect shift2 = a_region.bigEnd() - ibh;
  CH_assert(shift==shift2);

  int* iptr = (int*)rptr;
  int count_received = *iptr++;
  rptr = (Real*)iptr;
  for ( int i=0; i<count_received; i++ )
  {
    IntVect iv2;
    iptr = (int*)rptr;
    for ( int idir=0; idir<SpaceDim; idir++ ) iv2[idir] = *iptr++;
    rptr = (Real*)iptr;

    // perform shift
    iv2 += shift;

    for ( int c=a_comps.begin(); c<=a_comps.end(); c++ )
    {
      int iloc = getIndex(iv2, c);
      m_dataPtr[iloc] = *rptr++;
    }
  }
}

/********************/
template <> inline void IVSFAB<Real>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  IntVectSet subset = m_ivs & a_region;

  Real* rptr = (Real*)a_buf;
  for ( int idir=0; idir<SpaceDim; idir++ ) {
    *rptr++ = a_region.smallEnd(idir);
    *rptr++ = a_region.bigEnd(idir);
  }

  int* iptr = (int*)rptr;
  int count = subset.numPts();
  *iptr++ = count;
  rptr = (Real*)iptr;

  for ( IVSIterator ivsit(subset); ivsit.ok(); ++ivsit )
  {
    const IntVect& iv = ivsit();
    iptr = (int*)rptr;
    for ( int idir=0; idir<SpaceDim; idir++ ) *iptr++ = iv[idir];
    rptr = (Real*)iptr;
    for ( int c=a_comps.begin(); c<=a_comps.end(); c++ )
    {
      int iloc = getIndex(iv, c);
      *rptr++ = m_dataPtr[iloc];
    }
  }
  CH_assert( ((char*)rptr) - ((char*)a_buf) == size(a_region,a_comps) );
}


/******************/
//template specializations for int
template <> inline
int IVSFAB<int>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  // size of actual data plus size of count ..
  // specialization for T=int
  IntVectSet subset = m_ivs & a_region;
  int numberThings = subset.numPts();
  int retval = numberThings*sizeof(int)*a_comps.size() + sizeof(int);
  retval += SpaceDim*sizeof(int)*numberThings; // for intvects
  retval += 2*SpaceDim*sizeof(int); // box 
  return retval;
}


/********************/
template <> inline void IVSFAB<int>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  IntVectSet subset = m_ivs & a_region;

  int* rptr = (int*)a_buf;
  for ( int idir=0; idir<SpaceDim; idir++ ) 
    {
      *rptr++ = a_region.smallEnd(idir);
      *rptr++ = a_region.bigEnd(idir);
    }

  int* iptr = (int*)rptr;
  int count = subset.numPts();
  *iptr++ = count;
  rptr = (int*)iptr;

  for ( IVSIterator ivsit(subset); ivsit.ok(); ++ivsit )
  {
    const IntVect& iv = ivsit();
    iptr = (int*)rptr;
    for ( int idir=0; idir<SpaceDim; idir++ ) *iptr++ = iv[idir];
    rptr = (int*)iptr;
    for ( int c=a_comps.begin(); c<=a_comps.end(); c++ )
    {
      int iloc = getIndex(iv, c);
      *rptr++ = m_dataPtr[iloc];
    }
  }
  //  CH_assert( ((char*)rptr) - ((char*)a_buf) == size(a_region,a_comps) );
}
/********************/
template <> inline void IVSFAB<int>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  // message includes count
  IntVectSet subset = m_ivs & a_region;

  // input the box
  IntVect ibl, ibh;
  int* rptr = (int*)a_buf;
  for ( int idir=0; idir<SpaceDim; idir++ ) 
    {
      ibl[idir] = *rptr++;
      ibh[idir] = *rptr++;
    }

  // compute the offset (for periodic case)
  IntVect shift = a_region.smallEnd() - ibl;
  IntVect shift2 = a_region.bigEnd() - ibh;
  CH_assert(shift==shift2);

  int* iptr = (int*)rptr;
  int count_received = *iptr++;
  rptr = (int*)iptr;
  for ( int i=0; i<count_received; i++ )
  {
    IntVect iv2;
    iptr = (int*)rptr;
    for ( int idir=0; idir<SpaceDim; idir++ ) iv2[idir] = *iptr++;
    rptr = (int*)iptr;

    // perform shift
    iv2 += shift;

    for ( int c=a_comps.begin(); c<=a_comps.end(); c++ )
    {
      int iloc = getIndex(iv2, c);
      m_dataPtr[iloc] = *rptr++;
    }
  }
}
////////////////////////////////////////
template <class T> inline
int IVSFAB<T>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  IntVectSet subset = m_ivs & a_region;
  int numberThings = subset.numPts();

// can't do this for objects like vectors ...
//  T dummy;
//  int sizeofThing = dummy.linearSize();
//  int retval = numberThings*sizeofThing + sizeof(int);
 
  int retval=sizeof(int); // count
  for( int icomp=a_comps.begin(); icomp<=a_comps.end(); icomp++ ) {
    for( IVSIterator ivsit(subset); ivsit.ok(); ++ivsit ) {
      const IntVect& iv = ivsit();
      int iloc = getIndex(iv,icomp);
      retval += m_dataPtr[iloc].linearSize(); // size of each object
    }
  }
  retval += SpaceDim*sizeof(int)*numberThings; // for intvects
  retval += 2*SpaceDim*sizeof(Real); // box done as reals for some reason
  return retval;
}
////////////////////
template <class T> inline
void IVSFAB<T>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  // Note:  it WAS assumed that there is one component.
  //        and all "T" have equal size.
  // message includes count
  IntVectSet subset = m_ivs & a_region;


  Real* rptr = (Real*)a_buf;
  for ( int idir=0; idir<SpaceDim; idir++ ) {
    *rptr++ = a_region.smallEnd(idir);
    *rptr++ = a_region.bigEnd(idir);
  }

  int* iptr = (int*)rptr;

  int count = subset.numPts();
  *iptr++ = count;
// cannot count on T having constant size
//  T dummy;
//  int Tsz = dummy.linearSize();
  char* cptr = (char*)iptr;
  for ( IVSIterator ivsit(subset); ivsit.ok(); ++ivsit )
  {
    const IntVect& iv = ivsit();

    iptr = (int*)cptr;
    for ( int idir=0; idir<SpaceDim; idir++ ) *iptr++ = iv[idir];
    cptr = (char*)iptr;
    
    for( int c=a_comps.begin(); c<=a_comps.end(); c++ ) {
      int iloc = getIndex(iv, c); // assume 1 component ...
      m_dataPtr[iloc].linearOut( (void*)cptr );
      cptr += m_dataPtr[iloc].linearSize();
    }
  }
if( cptr - ((char*)a_buf) != size(a_region,a_comps) ) {
  pout() << "    region " << a_region << endl;
  pout() << "    comps " << a_comps.begin() << " to " << a_comps.end() << endl;
  pout() << "    size " << size(a_region,a_comps) << endl;
  pout() << "    ptrdiff  " << (cptr - ((char*)a_buf)) << endl;
}
  CH_assert( cptr - ((char*)a_buf) == size(a_region,a_comps) );
}

template <class T> inline
void IVSFAB<T>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  // Note:  it is assumed that there is one component.
  //        and all "T" have same size
  // message includes count
  IntVectSet subset = m_ivs & a_region;

  // input the box
  IntVect ibl, ibh;
  Real* rptr = (Real*)a_buf;
  for ( int idir=0; idir<SpaceDim; idir++ ) {
    ibl[idir] = *rptr++;
    ibh[idir] = *rptr++;
  }

  // compute the offset (for periodic case)
  IntVect shift = a_region.smallEnd() - ibl;
  IntVect shift2 = a_region.bigEnd() - ibh;
  CH_assert(shift==shift2);

  int* iptr = (int*)rptr;

  int count_received = *iptr++;
// cannot count on T having constant size
//  T dummy;
//  int Tsz = dummy.linearSize();
  char* cptr = (char*)iptr;
  for ( int i=0; i<count_received; i++ )
  {
    IntVect iv2;
    iptr = (int*)cptr;
    for ( int idir=0; idir<SpaceDim; idir++ ) iv2[idir] = *iptr++;
    cptr = (char*)iptr;

    // perform shift
    iv2 += shift;

    for( int c=a_comps.begin(); c<=a_comps.end(); c++ ) {
      int iloc = getIndex(iv2, c);
      m_dataPtr[iloc].linearIn( (void*)cptr );
//    cptr += Tsz;
      cptr += m_dataPtr[iloc].linearSize();
    }
  }
}

template <class T> inline
int
IVSFAB<T>::getIndex(const IntVect& a_iv, const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(m_ivs.contains(a_iv));
  CH_assert((a_comp >= 0) && (a_comp < m_nComp));

  int ioffset = m_ivmap(a_iv, 0);
  CH_assert(ioffset >= 0);
  CH_assert(ioffset < m_nIvs);
  //now add offset from componentnitude
  ioffset += m_nIvs*a_comp;
  return ioffset;
}
/********************/
template <class T> inline
void
IVSFAB<T>::clear()
{
  m_nComp = 0;
  m_nIvs = 0;
  m_ivs.makeEmpty();
  m_ivmap.clear();
  if (m_dataPtr != NULL)
    {
      delete[] m_dataPtr;
      m_dataPtr = NULL;
    }
  m_isDefined = false;
}
/*************************/
template <class T> inline
bool
IVSFAB<T>::isDefined() const
{
  return (m_isDefined);
}
/*************************/
template <class T> inline
int
IVSFAB<T>::numIvs() const
{
  return m_nIvs;
}
/*************************/
template <class T> inline
int
IVSFAB<T>::nComp() const
{
  return m_nComp;
}
/*************************/
template <class T> inline
T&
IVSFAB<T>::operator() (const IntVect& a_ndin,
                       const int& a_comp)
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  int iloc = getIndex(a_ndin, a_comp);
  return(m_dataPtr[iloc]);
}
/**************************/
template <class T> inline
const T&
IVSFAB<T>::operator() (const IntVect& a_ndin,
                       const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  int iloc = getIndex(a_ndin, a_comp);
  return(m_dataPtr[iloc]);
}
/******************/
template <class T> inline
const T*
IVSFAB<T>::dataPtr(const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  return m_dataPtr + a_comp*m_nIvs;
}
/******************/
template <class T> inline
T*
IVSFAB<T>::dataPtr(const int& a_comp)
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  return m_dataPtr + a_comp*m_nIvs;
}
/******************/
template <class T> inline
const int*
IVSFAB<T>::loVect() const
{
  return m_loVect.getVect();
}
/******************/
template <class T> inline
const int*
IVSFAB<T>::hiVect() const
{
  return m_hiVect.getVect();
}
/******************/
/*
template <class T> inline
void
IVSFAB<T>::setDefaultValues()
{
  m_isDefined = false;
  m_dataPtr = NULL;
  m_nIvs = 0;
  m_nComp = 0;
  m_loVect = IntVect::Unit;
  m_hiVect = IntVect::Zero;
}
*/
#include "NamespaceFooter.H"
#endif
