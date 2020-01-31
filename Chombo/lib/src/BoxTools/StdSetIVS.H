#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STDSETIVS_H_
#define _STDSETIVS_H_

#include "Box.H"
#include "ProblemDomain.H"
#include <set>
#include "NamespaceHeader.H"
using namespace std;
class StdSetIVSIterator;

struct lex_compare_iv
{
  bool operator() (const IntVect& lhs, const IntVect& rhs) const
  {
    return lhs.lexLT(rhs);
  }
};
/// IntVectSet implementation based std::set
/**
 */
class StdSetIVS
{
public:
  friend class StdSetIVSIterator; 

  ///
  StdSetIVS(){;}

  ///
  StdSetIVS(const Box& a_box);

  ///
  StdSetIVS(const StdSetIVS& a_sivs);

  ///
  ~StdSetIVS(){;}

  ///
  void define(const Box&);

  ///
  void define(const StdSetIVS& a_sivs);

  ///
  StdSetIVS& operator=(const StdSetIVS& a_sivs);

  ///or
  StdSetIVS& operator|=(const StdSetIVS& a_sivs);

  ///
  StdSetIVS& operator|=(const IntVect& a_iv);

  ///
  StdSetIVS& operator|=(const Box& a_box);

  ///and
  StdSetIVS& operator&=(const StdSetIVS& s_sivs);

  ///and
  StdSetIVS& operator&=(const Box& a_box);

  ///and
  StdSetIVS& operator&=(const ProblemDomain& a_domain);

  ///not
  StdSetIVS& operator-=(const StdSetIVS& a_sivs);

  ///not
  StdSetIVS& operator-=(const IntVect& a_iv);

  ///not
  StdSetIVS& operator-=(const Box& a_box);

  ///returns true if
  bool operator==(const StdSetIVS& lhs) const;

  ///
  bool contains(const IntVect& iv) const;

  ///
  bool contains(const Box& box) const;

  ///
  void grow(int igrow);

  ///
  void grow(int idir, int igrow);

  ///
  void growHi();

  ///
  void growHi(int a_dir);

  ///
  void refine(int iref = 2);

  ///
  void coarsen(int iref = 2);

  ///
  void shift(const IntVect& iv);

  ///
  void clear();

  ///
  Box minBox() const;

  ///
  bool isEmpty() const;

  ///
  void makeEmpty() ;

  ///
  int numPts() const;

  ///
  int linearSize() const;

  ///
  void linearIn(const void* const inBuf);

  ///
  void linearOut(void* const a_outBuf) const;

  void define(const Vector<IntVect>& a_vect);
private:

  void getVectorIV(Vector<IntVect>& a_vect) const;
  std::set<IntVect, lex_compare_iv> m_stdSet;
};

///
class StdSetIVSIterator
{
public:
  ///
  StdSetIVSIterator();

  ///
  StdSetIVSIterator(const StdSetIVS& ivs);

  ///
  void define(const StdSetIVS& ivs);

  ///
  const IntVect& operator()() const ;

  ///
  bool ok() const;

  ///
  void operator++();

  ///
  void begin();

  ///
  void end();

  ///
  void clear();
private:
  const StdSetIVS* m_ivs;
  std::set<IntVect, lex_compare_iv>::iterator m_iter;
};

#include "NamespaceFooter.H"
#endif //  STDSETIVS_H
