#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFFluxFAB.H"
#include "NamespaceHeader.H"

MFFluxFAB::~MFFluxFAB()
{
  for (int i=0; i<m_phase.size(); ++i)
    {
      delete m_phase[i];
    }
}

MFFluxFAB::MFFluxFAB(const Vector<EBISBox>& a_phaseGraphs,
                     const Box& a_region, const Vector<int>& a_comps)
{
  CH_assert(a_phaseGraphs.size()==  a_comps.size());
  m_phase.resize(a_phaseGraphs.size(), NULL);
  for (int i=0; i<a_phaseGraphs.size(); i++)
    {
      m_phase[i] = new EBFluxFAB(a_phaseGraphs[i], a_region, a_comps[i]);
    }
  m_box = a_region;
}

void MFFluxFAB::setVal(Vector<Real> a_value)
{
  CH_assert(a_value.size()==m_phase.size());

  for (int i=0; i<a_value.size(); i++)
    {
      m_phase[i]->setVal(a_value[i]);
    }
}

void MFFluxFAB::setVal(const Real& a_value)
{

  for (int i=0; i<m_phase.size(); i++)
    {
      m_phase[i]->setVal(a_value);
    }
}

void MFFluxFAB::copy(const Box& RegionFrom,
                     const Interval& destInt,
                     const Box& RegionTo,
                     const MFFluxFAB& source,
                     const Interval& srcInt)
{
  for (int i=0; i< m_phase.size(); ++i)
    {
      m_phase[i]->copy(RegionFrom,destInt,RegionTo,
                       *(source.m_phase[i]), srcInt);
    }
}

 ///
int MFFluxFAB::size(const Box& R, const Interval& comps) const
{
  int size = m_phase.size();
  for (int i=0; i<m_phase.size(); ++i)
    {
      size+=m_phase[i]->size(R, comps);
    }
  return size;
}
 ///
void MFFluxFAB::linearOut(void* buf, const Box& R,
                          const Interval& comps) const
{
  int* buffer = (int*)buf;
  for (int i=0; i<m_phase.size(); ++i)
    {
      *buffer = m_phase[i]->size(R, comps);
      ++buffer;
    }
  int* size = (int*)buf;
  unsigned char* ebbuffer = (unsigned char*)buffer;
  for (int i=0; i<m_phase.size(); ++i)
    {
      m_phase[i]->linearOut(ebbuffer, R, comps);
      ebbuffer+= size[i];
    }
}
  ///
void MFFluxFAB::linearIn(void* buf, const Box& R, const Interval& comps)
{
  int* size = (int*)buf;
  for (int i=0; i<m_phase.size(); ++i)
    {
      CH_assert(size[i] == m_phase[i]->size(R, comps));
    }
  unsigned char* ebbuffer  = (unsigned char*)(size+m_phase.size());
  for (int i=0; i<m_phase.size(); ++i)
    {
      m_phase[i]->linearIn(ebbuffer, R, comps);
      ebbuffer += size[i];
    }
}


MFFluxFAB&
MFFluxFAB::copy(const MFFluxFAB& a_src)
{
  for (int i=0;i<m_phase.size(); ++i)
    {
      EBFluxFAB& srcFAB = *a_src.m_phase[i];
      Interval comps(0, srcFAB.nComp()-1);
      Box thisBox = srcFAB.getRegion();
      m_phase[i]->copy(thisBox,
                       comps,
                       thisBox,
                       srcFAB,
                       comps);
    }
  return *this;
}

MFFluxFactory::MFFluxFactory(Vector<EBISLayout>& a_ebis,
                             const  Vector<int>& a_ncomp)
{
  define(a_ebis, a_ncomp);
}

void MFFluxFactory::define(Vector<EBISLayout>& a_ebis,
                           const  Vector<int>& a_ncomp)
{
  CH_assert(a_ebis.size() == a_ncomp.size());
  m_ebis = a_ebis;
  m_ncomp = a_ncomp;
}

MFFluxFactory::MFFluxFactory(const MFIndexSpace& a_mf,
                             const DisjointBoxLayout& a_dbl,
                             const Box& a_domain,
                             const Vector<int>& a_ncomps,
                             int ghost)
{
  int p = a_mf.numPhases();
  CH_assert(a_ncomps.size() == a_mf.numPhases());

  Vector<EBISLayout> layouts(p);
  a_mf.fillEBISLayout(layouts,
                      a_dbl,
                      a_domain,
                      ghost);
  define(layouts, a_ncomps);
}

MFFluxFactory::~MFFluxFactory()
{
}

MFFluxFAB* MFFluxFactory::create(const Box& a_box, int a_ncompsIgnored,
                                 const DataIndex& a_dit) const
{
  Vector<EBISBox> ebox(m_ebis.size());
  for (int i=0; i<ebox.size(); ++i)
    {
      ebox[i] = m_ebis[i][a_dit];
    }
  MFFluxFAB* mfab = new MFFluxFAB(ebox, a_box, m_ncomp);
  if (mfab == NULL)
  {
    MayDay::Error("Out of memory in MFFluxFactory::create");
  }
  return mfab;
}
#include "NamespaceFooter.H"
