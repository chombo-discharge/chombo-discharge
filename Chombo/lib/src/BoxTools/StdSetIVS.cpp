#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Box.H"
#include "SPMD.H"
#include "BoxIterator.H"
#include "StdSetIVS.H"
#include "NamespaceHeader.H"
///
StdSetIVS::
StdSetIVS(const Box& a_box)
{
  for(BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      *this |= bit();
    }
}
///
StdSetIVS::
StdSetIVS(const StdSetIVS& a_sivs)
{
  m_stdSet = a_sivs.m_stdSet; 
}
///
void 
StdSetIVS::
define(const Box& a_box)
{
  *this = StdSetIVS(a_box);
}
///
void 
StdSetIVS::
define(const StdSetIVS& a_sivs)
{
  m_stdSet = a_sivs.m_stdSet; 
}
///
StdSetIVS& 
StdSetIVS::
operator=(const StdSetIVS& a_sivs)
{
  m_stdSet = a_sivs.m_stdSet; 
  return *this;
}
///
StdSetIVS& 
StdSetIVS::
operator|=(const StdSetIVS& a_sivs)
{
  const std::set<IntVect, lex_compare_iv> inputset = a_sivs.m_stdSet;
  std::set<IntVect,lex_compare_iv>::iterator it;
  for(it = inputset.begin(); it!=  inputset.end(); ++it)
    {
      m_stdSet.insert(*it);
    }
  return *this;
}
///
StdSetIVS & 
StdSetIVS::
operator|=(const IntVect& a_iv)
{
  m_stdSet.insert(a_iv);
  return *this;
}
///
StdSetIVS& 
StdSetIVS::
operator|=(const Box& a_box)
{
  for(BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      m_stdSet.insert(bit());
    }
  return *this;
}
///
StdSetIVS& 
StdSetIVS::
operator&=(const StdSetIVS& a_sivs)
{
  if(&a_sivs != this)
    {
      std::set<IntVect, lex_compare_iv> newSet;
      std::set<IntVect, lex_compare_iv>::iterator it;
      for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
        {
          const IntVect& iv = *it;
          if(contains(iv) && a_sivs.contains(iv))
            {
              newSet.insert(iv);
            }
        }
      m_stdSet = newSet;
    }
  return *this;
}

///and
StdSetIVS& 
StdSetIVS::
operator&=(const Box& a_box)
{
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      if(!a_box.contains(iv))
        {
          m_stdSet.erase(it);
        }
    }
  return *this;
}
///and
StdSetIVS& 
StdSetIVS::
operator&=(const ProblemDomain& a_domain)
{
  return (*this &= a_domain.domainBox());
}
///not
StdSetIVS& 
StdSetIVS::
operator-=(const StdSetIVS& a_sivs)
{
  std::set<IntVect, lex_compare_iv>::iterator it;
  //leaving out the ++it because  erase 
  for(it = m_stdSet.begin(); it!=  m_stdSet.end(); )
    {
      if(a_sivs.contains(*it))
        {
          m_stdSet.erase(it++);
        }
      else
        {
          ++it;
        }
    }
  return *this;
}
///not
StdSetIVS& 
StdSetIVS::
operator-=(const IntVect& a_iv)
{
  if(contains(a_iv))
    {
      m_stdSet.erase(m_stdSet.find(a_iv));
    }
  return *this;
}
///not
StdSetIVS& 
StdSetIVS::
operator-=(const Box& a_box)
{
  for(BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      *this -= bit();
    }
  return *this;
}
///
bool 
StdSetIVS::
operator==(const StdSetIVS& a_lhs) const
{
  if(a_lhs.m_stdSet.size() != m_stdSet.size())
    {
      return false;
    }

  bool retval = true;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
    {
      if((!contains(*it)) || (!a_lhs.contains(*it)))
        {
          retval = false;
          break;
        }
    }
  return retval;
}

///
bool 
StdSetIVS::
contains(const IntVect& a_iv) const
{
  std::set<IntVect, lex_compare_iv>::iterator it = m_stdSet.find(a_iv);
  return (it != m_stdSet.end());
}

///
bool 
StdSetIVS::
contains(const Box& a_box) const
{
  bool retval = true;
  for(BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      if(!contains(bit()))
        {
          retval = false;
          break;
        }
    }
  return retval;
}

///
void 
StdSetIVS::
grow(int igrow)
{
  StdSetIVS newSet;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.grow(igrow);
      newSet |= grid;
    }
  *this = newSet;
}

///
void 
StdSetIVS::
grow(int idir, int igrow)
{
  StdSetIVS newSet;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.grow(idir, igrow);
      newSet |= grid;
    }
  *this = newSet;
}

///
void 
StdSetIVS::
growHi()
{
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      growHi(idir);
    }
}

///
void 
StdSetIVS::
growHi(int a_dir)
{
  StdSetIVS newSet;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.growHi(a_dir);
      newSet |= grid;
    }
  *this = newSet;
}

///
void 
StdSetIVS::
refine(int iref)
{
  StdSetIVS newSet;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.refine(iref);
      newSet |= grid;
    }
  *this = newSet;
}

///
void 
StdSetIVS::
coarsen(int iref)
{
  std::set<IntVect, lex_compare_iv> newSet;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.coarsen(iref);
      for(BoxIterator bit(grid); bit.ok(); ++bit)
        {
          newSet.insert(bit());
        }
    }
  m_stdSet = newSet;
}

///
void 
StdSetIVS::
shift(const IntVect& a_iv)
{
  std::set<IntVect, lex_compare_iv> newSet;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      IntVect iv = *it;
      iv.shift(a_iv);
      newSet.insert(iv);
    }
  m_stdSet = newSet;
}

///
void 
StdSetIVS::
clear()
{
  std::set<IntVect, lex_compare_iv> newSet;
  m_stdSet = newSet;
}

///
Box 
StdSetIVS::
minBox() const
{
  int bignum = 100000;
  IntVect lo = bignum*IntVect::Unit;
  IntVect hi =-bignum*IntVect::Unit;
  std::set<IntVect, lex_compare_iv>::iterator it;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          lo[idir] = Min(lo[idir], iv[idir]);
          hi[idir] = Max(hi[idir], iv[idir]);
        }
    }
  
  Box retval(lo, hi);
  return retval;
}

///
bool 
StdSetIVS::
isEmpty() const
{
  return (m_stdSet.size() == 0);
}
///
void 
StdSetIVS::
getVectorIV(Vector<IntVect>& a_vect) const
{
  a_vect.resize(m_stdSet.size());
  int bignum = 100000;
  IntVect lo = bignum*IntVect::Unit;
  IntVect hi =-bignum*IntVect::Unit;
  std::set<IntVect, lex_compare_iv>::iterator it;
  int ivec = 0;
  for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      a_vect[ivec] = *it;
      ivec++;
    }
}
///
void 
StdSetIVS::
makeEmpty() 
{
  clear();
}

///
int
StdSetIVS::
numPts() const
{
  return m_stdSet.size();
}
///
int 
StdSetIVS::
linearSize() const
{
  Vector<IntVect> vect;
  getVectorIV(vect);
  int retval = linearListSize(vect);
  return retval;
}

///
void 
StdSetIVS::
define(const Vector<IntVect>& a_vect)
{
  makeEmpty();
  for(int ivec = 0; ivec  < a_vect.size(); ivec++)
    {
      (*this ) |= a_vect[ivec];
    }
}

///
void 
StdSetIVS::
linearIn(const void* const inBuf)
{
  Vector<IntVect> vect;
  linearListIn(vect, inBuf);
  define(vect);
}

///
void 
StdSetIVS::
linearOut(void* const a_outBuf) const
{
  Vector<IntVect> vect;
  getVectorIV(vect);
  linearListOut(a_outBuf, vect);
}


StdSetIVSIterator::
StdSetIVSIterator()
{
  m_ivs = NULL;
}

///
StdSetIVSIterator::
StdSetIVSIterator(const StdSetIVS& ivs)
{
  m_ivs = &ivs;
  m_iter = m_ivs->m_stdSet.begin();
}

///
void 
StdSetIVSIterator::
define(const StdSetIVS& a_ivs)
{
  m_ivs = &a_ivs;
  m_iter = m_ivs->m_stdSet.begin();
}

///
const IntVect& 
StdSetIVSIterator::
operator()() const 
{
  return *m_iter;
}

///
bool 
StdSetIVSIterator::
ok() const
{
  CH_assert(m_ivs != NULL);
  return (m_iter != m_ivs->m_stdSet.end());
}

///
void 
StdSetIVSIterator::
operator++()
{
  m_iter++;
}

///
void 
StdSetIVSIterator::
begin()
{
  CH_assert(m_ivs != NULL);
  m_iter = m_ivs->m_stdSet.begin();
}

///
void 
StdSetIVSIterator::
end()
{
  CH_assert(m_ivs != NULL);
  m_iter = m_ivs->m_stdSet.end();
}

///
void 
StdSetIVSIterator::
clear()
{
  m_ivs = NULL;
}

#include "NamespaceFooter.H"
