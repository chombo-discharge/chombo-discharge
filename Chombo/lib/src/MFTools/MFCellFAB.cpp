#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFCellFAB.H"
#include "NamespaceHeader.H"

MFCellFAB::~MFCellFAB()
{
  for (int i=0; i<m_phase.size(); ++i)
    {
      delete m_phase[i];
    }
}

MFCellFAB::MFCellFAB(const Vector<EBISBox>& a_phaseGraphs,
                     const Box& a_region, const Vector<int>& a_comps)
{
  CH_assert(a_phaseGraphs.size()==  a_comps.size());
  m_phase.resize(a_phaseGraphs.size(), NULL);
  for (int i=0; i<a_phaseGraphs.size(); i++)
    {
      m_phase[i] = new EBCellFAB(a_phaseGraphs[i], a_region, a_comps[i]);
    }
  m_box = a_region;
}

void MFCellFAB::setVal(Vector<Real> a_value)
{
  CH_assert(a_value.size()==m_phase.size());

  for (int i=0; i<a_value.size(); i++)
    {
      m_phase[i]->setVal(a_value[i]);
    }
}

void MFCellFAB::setVal(const Real& a_value)
{

  for (int i=0; i<m_phase.size(); i++)
    {
      m_phase[i]->setVal(a_value);
    }
}

void MFCellFAB::setVal(Real a_value, const Box& region,
                       int startcomp, int ncomp)
{

  for (int i=0; i<m_phase.size(); i++)
    {
      m_phase[i]->setVal(a_value, region, startcomp, ncomp);
    }
}

void MFCellFAB::copy(const Box& RegionFrom,
          const Interval& destInt,
          const Box& RegionTo,
          const MFCellFAB& source,
          const Interval& srcInt)
{
  for (int i=0; i< m_phase.size(); ++i)
    {
      m_phase[i]->copy(RegionFrom,destInt,RegionTo,
                       *(source.m_phase[i]), srcInt);
    }
}

 ///
int MFCellFAB::size(const Box& R, const Interval& comps) const
{
  int size = m_phase.size()*sizeof(int);
  for (int i=0; i<m_phase.size(); ++i)
    {
      size+=m_phase[i]->size(R, comps);
    }
  return size;
}
 ///
void MFCellFAB::linearOut(void* buf, const Box& R,
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
void MFCellFAB::linearIn(void* buf, const Box& R, const Interval& comps)
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
///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator+=(const MFCellFAB& a_mffab)
{
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] += *a_mffab.m_phase[i];
    }
  return *this;
}


///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB&
MFCellFAB::plus(const MFCellFAB& a_mffab,
                int a_srccomp,
                int a_destcomp,
                int a_numcomp)
{
  for (int i=0;i<m_phase.size(); ++i)
    {
      m_phase[i]->plus(*(a_mffab.m_phase[i]), a_srccomp,
                       a_destcomp, a_numcomp);
    }
  return *this;
}

MFCellFAB&
MFCellFAB::plus(const MFCellFAB& a_src,
                const Box&       a_srcbox,
                const Box&       a_destbox,
                const Real&      a_scale,
                int              a_srccomp,
                int              a_destcomp,
                int              a_numcomp)
{
  CH_assert(a_srcbox == a_destbox);
  CH_assert(a_srcbox == m_box);
  return this->plus(a_src, a_scale);
}

///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator-=(const MFCellFAB& a_mffab)
 {
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] -= *a_mffab.m_phase[i];
    }
  return *this;
}

///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB&
MFCellFAB::minus(const MFCellFAB& a_mffab,
                 int a_srccomp,
                 int a_destcomp,
                 int a_numcomp)
{
  for (int i=0;i<m_phase.size(); ++i)
    {
      m_phase[i]->minus(*(a_mffab.m_phase[i]), a_srccomp,
                        a_destcomp, a_numcomp);
    }
  return *this;
}


///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator*=(const MFCellFAB& a_mffab)
 {
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] *= *a_mffab.m_phase[i];
    }
  return *this;
}

///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB&
MFCellFAB::mult(const MFCellFAB& a_mffab,
                int a_srccomp,
                int a_destcomp,
                int a_numcomp)
{
  for (int i=0;i<m_phase.size(); ++i)
    {
      m_phase[i]->mult(*(a_mffab.m_phase[i]), a_srccomp,
                       a_destcomp, a_numcomp);
    }
  return *this;
}


///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator/=(const MFCellFAB& a_mffab)
 {
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] /= *a_mffab.m_phase[i];
    }
  return *this;
}

///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB&
MFCellFAB::divide(const MFCellFAB& a_mffab,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp)
{
  for (int i=0;i<m_phase.size(); ++i)
    {
      m_phase[i]->divide(*(a_mffab.m_phase[i]), a_srccomp,
                         a_destcomp, a_numcomp);
    }
  return *this;
}



///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator+=(const Real& a_scalar)
 {
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] += a_scalar;
    }
  return *this;
}
///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator-=(const Real& a_scalar)
  {
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] -= a_scalar;
    }
  return *this;
}
///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator*=(const Real& a_scalar)
   {
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] *= a_scalar;
    }
  return *this;
}
///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::mult(Real a_scalar)
   {
  for (int i=0;i<m_phase.size(); ++i)
    {
      m_phase[i]->mult(a_scalar);
    }
  return *this;
}
///
/**
   Both fabs need the same ebisBox and region.
*/
MFCellFAB& MFCellFAB::operator/=(const Real& a_scalar)
    {
  for (int i=0;i<m_phase.size(); ++i)
    {
      *m_phase[i] /= a_scalar;
    }
  return *this;
}
///
/**
   Current FAB += a_src FAB * a_scalar.  Both fabs need the same ebisBox
   and region.
*/
MFCellFAB& MFCellFAB::plus(const MFCellFAB& a_src,
                           Real             a_scalar)
 {
  for (int i=0;i<m_phase.size(); ++i)
    {
      m_phase[i]->plus(*a_src.m_phase[i], a_scalar);
    }
  return *this;
}

MFCellFAB& MFCellFAB::copy(const MFCellFAB& a_src)
 {
  for (int i=0;i<m_phase.size(); ++i)
    {
      m_phase[i]->copy(*a_src.m_phase[i]);
    }
  return *this;
}

/// (Not implemented) Returns the Lp-norm of this MFCellFAB
/**
   (Not implemented) Returns the Lp-norm of this MFCellFAB using components
   (a_comp : a_comp + a_numcomp - 1).  a_power < 0 -> ERROR.
   a_power = 0  -> infinity norm (max norm).
   a_power = 1  -> L1-norm
   a_power > 1  -> Lp-norm
*/
 Real MFCellFAB::norm(int a_power,
                  int a_comp,
                  int a_numComp) const
{
  MayDay::Error("Not implemented\n");
  return -1;
}

/// (Not implemented) Returns the Lp-norm of this MFCellFAB within a region
  /**
     (Not implemented) Returns the Lp-norm of this MFCellFAB using components
     (a_comp : a_comp + a_numcomp - 1) and within the a_subbox.  a_power < 0
     -> ERROR.
     a_power = 0 -> infinity norm (max norm).
     a_power = 1 -> L1-norm
     a_power > 1 -> Lp-norm
  */
 Real MFCellFAB::norm(const Box& a_subbox,
                  int        a_power,
                  int        a_comp,
                    int        a_numComp) const
{
   MayDay::Error("Not implemented\n");
   return -1;
}

/// (Not implemented) Returns a sum of powers of a subset of this MFCellFAB
/**
     (Not implemented) Returns a sum of powers of a subset of this MFCellFAB,
     specifically components a_comp to a_comp+a_numcomp-1 within a_subbox.
     a_power >= 2 only.

*/
 Real MFCellFAB::sumPow(const Box& a_subbox,
                    int        a_power,
                    int        a_comp,
                      int        a_numComp) const
{
   MayDay::Error("Not implemented\n");
   return -1;
}

/// (Not implemented) Return the dot product of this MFCellFAB with another
/**
   (Not implemented) Return the dot product of this MFCellFAB and "ebfab2"
   over their overlap region and all components.
*/
Real MFCellFAB::dotProduct(const MFCellFAB& ebfab2) const

{
   MayDay::Error("Not implemented\n");
   return -1;
}

Real MFCellFAB::dotProduct(const MFCellFAB& ebfab2, const Box& a_box) const

{
   MayDay::Error("Not implemented\n");
   return -1;
}

MFCellFactory::MFCellFactory(Vector<EBISLayout>& a_ebis,
                             const  Vector<int>& a_ncomp)
{
  define(a_ebis, a_ncomp);
}

void MFCellFactory::define(Vector<EBISLayout>& a_ebis,
                           const  Vector<int>& a_ncomp)
{
  CH_assert(a_ebis.size() == a_ncomp.size());
  m_ebis = a_ebis;
  m_ncomp = a_ncomp;
}

MFCellFactory::MFCellFactory(const MFIndexSpace& a_mf,
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

MFCellFactory::~MFCellFactory()
{
}

MFCellFAB* MFCellFactory::create(const Box& a_box, int a_ncompsIgnored,
                                 const DataIndex& a_dit) const
{
  Vector<EBISBox> ebox(m_ebis.size());
  for (int i=0; i<ebox.size(); ++i)
    {
      ebox[i] = m_ebis[i][a_dit];
    }
  MFCellFAB* mfab = new MFCellFAB(ebox, a_box, m_ncomp);
  if (mfab == NULL)
  {
    MayDay::Error("Out of memory in MFCellFactory::create");
  }
  return mfab;
}
#include "NamespaceFooter.H"
