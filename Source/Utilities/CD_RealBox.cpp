/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RealBox.cpp
  @brief  Implementation of CD_RealBox.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_RealBox.H>
#include <CD_NamespaceHeader.H>

RealBox::RealBox()
{
  m_lo = RealVect::Zero;
  m_hi = RealVect::Zero;
}

RealBox::RealBox(const RealVect a_lo, const RealVect a_hi)
{
  m_lo = a_lo;
  m_hi = a_hi;
}

RealBox::RealBox(const Box a_box, const RealVect a_origin, const Real a_dx)
{

  const IntVect lo = a_box.smallEnd();
  const IntVect hi = a_box.bigEnd();

  m_lo = a_origin + a_dx * RealVect(lo);
  m_hi = a_origin + a_dx * RealVect(hi);
}

RealBox::~RealBox()
{}

RealVect
RealBox::getLo() const
{
  return m_lo;
}

RealVect
RealBox::getHi() const
{
  return m_hi;
}

Vector<RealVect>
RealBox::getCorners() const
{

  Vector<RealVect> corners(pow(2, SpaceDim));

  const RealVect DD = m_hi - m_lo;
  const RealVect DX = RealVect(D_DECL(DD[0], 0.0, 0.0));
  const RealVect DY = RealVect(D_DECL(0.0, DD[1], 0.0));
#if CH_SPACEDIM == 3
  const RealVect DZ = RealVect(D_DECL(0.0, 0.0, DD[2]));
#endif

  corners[0] = m_lo;
  corners[1] = m_lo + DX;
  corners[2] = m_lo + DY;
  corners[3] = m_lo + DX + DY; // = m_hi
#if CH_SPACEDIM == 3
  corners[4] = m_lo + DZ;
  corners[5] = m_lo + DZ + DX;
  corners[6] = m_lo + DZ + DY;
  corners[7] = m_lo + DX + DX + DY; // = m_hi
#endif

  return corners;
}

bool
RealBox::intersect(const RealBox& a_box) const
{

  //  bool ret = false;

  const RealVect LO = a_box.getLo();
  const RealVect HI = a_box.getHi();

  if (LO[0] < m_hi[0] && HI[0] > m_lo[0] //  Input x-edges either to left or righ
      && LO[1] < m_hi[1] && HI[1] > m_lo[1]
#if CH_SPACEDIM == 3
      && LO[2] < m_hi[2] && HI[2] > m_lo[2]
#endif
  ) {
    return true;
  }
  else {
    return false;
  }

  // const Vector<RealVect> corners = a_box.getCorners();

  // for (int i = 0; i < corners.size(); i++){
  //   if(isPointInside(corners[i])){
  //     ret = true;
  //   }
  // }

  // return ret;
}

bool
RealBox::isPointInside(const RealVect a_point) const
{

  bool ret = false;

  if (a_point[0] >= m_lo[0] && a_point[0] <= m_hi[0] && a_point[1] >= m_lo[1] && a_point[1] <= m_hi[1]
#if CH_SPACEDIM == 3
      && a_point[2] >= m_lo[2] && a_point[2] <= m_hi[2]
#endif
  ) {
    ret = true;
  }

  return ret;
}

bool
RealBox::isBoxInside(const RealBox& a_box) const
{

  Vector<RealVect> corners = a_box.getCorners();

  bool inside = true;
  return true;
  for (int i = 0; i < corners.size(); i++) {

    // Check if any of the corners of the input box lies inside this box. If one of the
    // corners of the input box is outside, a_box is NOT completely contained by this box.
    const bool is_corner_inside = this->isPointInside(corners[i]);

    if (!is_corner_inside) {
      inside = false;
      break;
    }
  }

  return inside;
}

#include <CD_NamespaceFooter.H>
