#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BinItem.H"
#include "NamespaceHeader.H"

BinItem::BinItem()
  : m_position(D_DECL6(1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8))
{
}

BinItem::BinItem(const RealVect& a_position)
  : m_position(a_position)
{
}

BinItem::~BinItem()
{
}

void BinItem::define(const RealVect& a_position)
{
  m_position = a_position;
}

void BinItem::setPosition(const RealVect& a_position)
{
  m_position = a_position;
}

void BinItem::setPosition(const Real a_position, const int a_dimension)
{
  m_position[a_dimension] = a_position;
}

const RealVect& BinItem::position() const
{
  return m_position;
}

RealVect& BinItem::position()
{
  return m_position;
}

Real BinItem::position(const int a_dir) const
{
  return m_position[a_dir];
}

bool BinItem::operator == (const BinItem& a_p) const
{
  return ( m_position  == a_p.m_position );
}

bool BinItem::operator == (const BinItem* a_p) const
{
  return (*this == *a_p);
}

bool BinItem::operator != (const BinItem& a_p) const
{
  return !(*this == a_p);
}

int BinItem::size() const
{
  return SpaceDim*sizeof(Real);
}

void BinItem::linearOut(void* buf) const
{
  Real* buffer = (Real*) buf;
  for (int dim=0; dim<SpaceDim; dim++)
    {
      *buffer = m_position[dim];
      ++buffer;
    }
}

void BinItem::linearIn(void* buf)
{
  Real* buffer = (Real*) buf;
  for (int dir=0; dir<SpaceDim; ++dir)
    {
      m_position[dir] = *buffer;
      ++buffer;
    }
}

// Write a human-readable version to stdout.  This is for debugging.
void BinItem::print(void) const
{
  std::cout << *this;
}

std::ostream& operator<<( std::ostream& ostr, const BinItem& p)
{
  ostr << "(BinItem (position " << p.position(0);
  for (int i=1; i<SpaceDim; ++i)
    {
      ostr << " " << p.position(i);
    }
  ostr << "))";

  return ostr;
}
#include "NamespaceFooter.H"
