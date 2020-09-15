#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MINIIVFAB_H_
#define _MINIIVFAB_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "SPACE.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "VolIndex.H"
#include "BaseFab.H"
#include "EBGraph.H"
#include "BaseIVFAB.H"
#include "NamespaceHeader.H"

///
/**
   MiniIVFAB is a templated
   data holder defined at the VoFs of an irregular domain.

Implemented as just a raw vector of vofs and data, more optimized
for smaller memory footprint and faster linearIn/linearOut.  will
be more brutal for vof-by-vof indexing.
bvs
*/
template <class T>
class MiniIVFAB //: public BaseIVFAB<T>
{
public:
  ///
  /**
     Null constructor.  
  */
  MiniIVFAB();

  ///
  /**
     Defining constructor.  Specifies the irregular domain
     and the number of data components per VoF.  The irregular domain
     must lie completely within the EBGraph.  The
     contents are uninitialized.  Calls full define function.
  */
  MiniIVFAB( const IntVectSet& a_region,
             const EBGraph&    a_ebgraph,
             const int&        a_nvarin);

  ///
  bool hasVoFs() const
  {
    if(m_vofs) return true;
    return false;
  }
  ///
  /**

     This aliased MiniIVFAB will have a_comps.size() components, starting at zero.
  */
  MiniIVFAB(const Interval& a_comps,
            MiniIVFAB<T>&   a_original);


  ///
  virtual ~MiniIVFAB();

  ///
  /**
     Full define function.  Specifies the irregular domain and the
     number of data components per VoF.  The irregular domain
     must lie completely within the EBGraph.   The contents
     are uninitialized.  If it has previously been defined, the old definition
     data is overwritten and lost.
  */
  virtual void define(const IntVectSet& a_region,
              const EBGraph&    a_ebgraph,
              const int&        a_nvarin);


  ///
  void copy(const Box& a_fromBox,
            const Interval& a_dstInterval,
            const Box& a_toBox,
            const MiniIVFAB<T>& a_src,
            const Interval& a_srcInterval);
  ///
  int size(const Box& R, const Interval& comps) const ;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  Vector<VolIndex>  getVoFs() const;

  inline T* dataPtr(int a_ivar)
  {
    if(m_vofs) return m_data+a_ivar*m_vofs->size();
    return m_data;
  }
  
  const T* dataPtr(int a_ivar) const
  {
    if(m_vofs) return m_data+a_ivar*m_vofs->size();
    return m_data;
  }
  
  void setVal(const T& a_val);

  void setVal(int a_comp, const T& a_val);

  void setVal(const T& a_val, const Box& a_box, int a_startcomp, int a_ncomp);


  ///get index into vector
  /**
     needs to be public so that bulk stencils can be constructed
   */
  virtual T* getIndex(const VolIndex& a_vof,const int& a_comp) const;

 ///for AggStencil
  long offset(const BaseIndex& a_vof, const int& a_ivar) const
  {
    const VolIndex* vof = dynamic_cast< const VolIndex* >(&a_vof);
    if (vof == NULL) MayDay::Error("cast failed:MiniIVFAB only takes vofs for indexing");
    unsigned int i = 0;
    for(;i<m_vofs->size(); i++)
      if(*vof == (*m_vofs)[i]) break;
    i+=a_ivar*m_vofs->size();
    return i;
  }

  ///
  T& operator() (const VolIndex& a_vof, int a_comp)
  {  return *(getIndex(a_vof, a_comp));}

  ///
  const T& operator() (const VolIndex& a_vof, int a_comp) const
  {
    return *(getIndex(a_vof, a_comp));
  }

  template <typename F>
  void forall(const MiniIVFAB& a_src, const Box& a_box, int a_srccomp, int a_destcomp, int a_numcomp, bool sameRegBox, const F& func);

  int numVoFs() const;

  int nComp() const {return m_nComp;}
private:

protected:

  shared_ptr<Vector<VolIndex> > m_vofs;
  shared_ptr<Vector<T> > m_Memory;
  int m_nComp = 0;
  T* m_data   = nullptr;

};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "MiniIVFABI.H"
#endif

#endif
