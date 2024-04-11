/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GeoCoarsener.cpp
  @brief  Implementation of CD_GeoCoarsener.H
  @author Robert Marskar
*/

// Std includes
#include <sstream>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_GeoCoarsener.H>
#include <CD_NamespaceHeader.H>

GeoCoarsener::GeoCoarsener()
{
  CH_TIME("GeoCoarsener::GeoCoarsener()");

  // Reset boxes
  m_coarsenBoxes.resize(0);
  m_coarsenLevels.resize(0);

  // Parse information from input script
  ParmParse pp("GeoCoarsener");

  int num_boxes = 0;
  pp.query("num_boxes", num_boxes);

  if (num_boxes > 0) {
    m_coarsenBoxes.resize(num_boxes);
    m_coarsenLevels.resize(num_boxes);
    m_inverse.resize(num_boxes, 0);

    // Read string in type format GeoCoarsener.box1_lo and GeoCoarsener.box2_hi etc.
    for (int ibox = 0; ibox < num_boxes; ibox++) {
      const std::string boxNum = std::to_string(1 + ibox);
      const std::string str1   = std::string("box") + boxNum + std::string("_lo");
      const std::string str2   = std::string("box") + boxNum + std::string("_hi");
      const std::string str3   = std::string("box") + boxNum + std::string("_lvl");
      const std::string str4   = std::string("box") + boxNum + std::string("_inv");

      Vector<Real> loCorner(SpaceDim);
      Vector<Real> hiCorner(SpaceDim);
      int          finestBoxLvl;
      bool         inverse;

      pp.getarr(str1.c_str(), loCorner, 0, SpaceDim);
      pp.getarr(str2.c_str(), hiCorner, 0, SpaceDim);
      pp.get(str3.c_str(), finestBoxLvl);
      pp.get(str4.c_str(), inverse);

      const RealVect c1 = RealVect(D_DECL(loCorner[0], loCorner[1], loCorner[2]));
      const RealVect c2 = RealVect(D_DECL(hiCorner[0], hiCorner[1], hiCorner[2]));

      m_coarsenBoxes[ibox]  = RealBox(c1, c2);
      m_coarsenLevels[ibox] = finestBoxLvl;
      m_inverse[ibox]       = inverse;
    }
  }
}

GeoCoarsener::~GeoCoarsener()
{
  CH_TIME("GeoCoarsener::~GeoCoarsener");
}

void
GeoCoarsener::coarsenTags(Vector<IntVectSet>& a_tags, const Vector<Real>& a_dx, const RealVect& a_probLo) const
{
  CH_TIME("GeoCoarsener::coarsenTags");

  if (!(m_coarsenBoxes.size() == m_coarsenLevels.size())) {
    pout() << "GeoCoarsener::coarsenTags - m_geoCoarsen is not well defined. Skipping the coarsening step" << endl;
  }
  else {
    if (m_coarsenBoxes.size() > 0) {
      for (int lvl = 0; lvl < a_tags.size(); lvl++) {
        const int        numCoarsen = m_coarsenBoxes.size();
        const IntVectSet tmp        = a_tags[lvl];

        for (IVSIterator it(tmp); it.ok(); ++it) {
          const IntVect  iv  = it();
          const RealVect pos = a_probLo + RealVect(iv) * a_dx[lvl] + 0.5 * a_dx[lvl] * RealVect::Unit;

          bool insideCoarsenBox  = false;
          bool insideInverseBox  = false;
          bool outsideInverseBox = false;

          // Check the regular box
          for (int ibox = 0; ibox < numCoarsen; ibox++) {
            const RealVect lo = m_coarsenBoxes[ibox].getLo();
            const RealVect hi = m_coarsenBoxes[ibox].getHi();

            const bool inverse    = (m_inverse[ibox] == 1);
            const bool insideBox  = (pos > lo) && (pos < hi);
            const bool coarsenLvl = lvl >= m_coarsenLevels[ibox];

            if (insideBox && !inverse && coarsenLvl) {
              insideCoarsenBox = true;
            }
          }

          // Check the inverse boxes
          for (int ibox = 0; ibox < numCoarsen; ibox++) {
            const RealVect lo = m_coarsenBoxes[ibox].getLo();
            const RealVect hi = m_coarsenBoxes[ibox].getHi();

            const bool inverse    = (m_inverse[ibox] == 1);
            const bool insideBox  = pos > lo && pos < hi;
            const bool coarsenLvl = lvl >= m_coarsenLevels[ibox];

            if (insideBox && inverse && coarsenLvl) { // Protect tag
              insideInverseBox = true;
            }
            else if (!insideBox && inverse && coarsenLvl) { // Remove tag
              outsideInverseBox = true;
            }
          }

          if (insideCoarsenBox && !insideInverseBox) {
            a_tags[lvl] -= iv;
          }
          if (!insideInverseBox && outsideInverseBox) {
            a_tags[lvl] -= iv;
          }
        }
      }
    }
  }
}

Vector<RealBox>
GeoCoarsener::getCoarsenBoxes() const
{
  CH_TIME("GeoCoarsener::getCoarsenBoxes()");

  return m_coarsenBoxes;
}

Vector<int>
GeoCoarsener::getCoarsenLevels() const
{
  CH_TIME("GeoCoarsener::getCoarsenLevels()");

  return m_coarsenLevels;
}

#include <CD_NamespaceFooter.H>
