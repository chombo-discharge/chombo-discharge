/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CellTagger.cpp
  @brief  Implementation of CD_CellTagger.H
  @author Robert marskar
*/

// Std includes
#include <string>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CellTagger.H>
#include <CD_NamespaceHeader.H>

CellTagger::CellTagger()
{
  CH_TIME("CellTagger::CellTagger()");

  // Defaeult settings
  m_verbosity = -1;
  m_buffer    = 0;
  m_name      = "CellTagger";
}

CellTagger::~CellTagger() { CH_TIME("CellTagger::~CellTagger()"); }

void
CellTagger::preRegrid() noexcept
{
  CH_TIME("CellTagger::preRegrid()");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid()" << endl;
  }
}

void
CellTagger::prePlot() const noexcept
{
  CH_TIME("CellTagger::prePlot()");
  if (m_verbosity > 5) {
    pout() << m_name + "::prePlot()" << endl;
  }
}

int
CellTagger::getNumberOfPlotVariables() const
{
  CH_TIME("CellTagger::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  return 0;
}

Vector<std::string>
CellTagger::getPlotVariableNames() const
{
  CH_TIME("CellTagger::getPlotVariableNames)");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames" << endl;
  }

  return Vector<std::string>();
}

void
CellTagger::writePlotData(LevelData<EBCellFAB>& a_output,
                          int&                  a_icomp,
                          const std::string     a_outputRealm,
                          const int             a_level) const
{
  CH_TIME("CellTagger::writePlotData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }
}

int
CellTagger::getBuffer() const
{
  CH_TIME("CellTagger::getBuffer()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getBuffer()" << endl;
  }

  return (m_buffer);
}

void
CellTagger::parseRuntimeOptions()
{
  CH_TIME("CellTagger::parseRunTimeOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }
}

void
CellTagger::parseTagBoxes()
{
  CH_TIME("CellTagger::parseTagBoxes()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseTagBoxes()" << endl;
  }

  // TLDR: This code parses strings from the command line in the form
  //
  //    m_name.num_boxes = 0
  //
  // If num_boxes > 0 we process lines
  //
  //    m_name.tag_box1_lo = 0 0 0
  //    m_name.tag_box1_hi = 1 1 1
  //    m_name.tag_box2_lo = 0 0 0
  //    m_name.tag_box2_hi = 1 1
  //
  // We then populate m_tagBoxes with RealBox (an axis-aligned box with physical coordinates). These
  // boxes can be used in tagCells(...) to prune regions in space where tagging is unnecessary. Users
  // are free to leave num_boxes=0 if they don't want to that.

  ParmParse pp(m_name.c_str());

  int numBoxes = 0;
  pp.query("num_tag_boxes", numBoxes);

  m_tagBoxes.resize(0);

  if (numBoxes > 0) {
    m_tagBoxes.resize(numBoxes);

    for (int ibox = 0; ibox < numBoxes; ibox++) {
      const std::string str1 = "tag_box" + std::to_string(1 + ibox) + "_lo";
      const std::string str2 = "tag_box" + std::to_string(1 + ibox) + "_hi";

      Vector<Real> cornerLo(SpaceDim);
      Vector<Real> cornerHi(SpaceDim);

      pp.getarr(str1.c_str(), cornerLo, 0, SpaceDim);
      pp.getarr(str2.c_str(), cornerHi, 0, SpaceDim);

      const RealVect c1 = RealVect(D_DECL(cornerLo[0], cornerLo[1], cornerLo[2]));
      const RealVect c2 = RealVect(D_DECL(cornerHi[0], cornerHi[1], cornerHi[2]));

      m_tagBoxes[ibox] = RealBox(c1, c2);
    }
  }
}

void
CellTagger::parseRefinementBoxes()
{
  CH_TIME("CellTagger::parseRefinementBoxes()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRefinementBoxes()" << endl;
  }

  // TLDR: This code parses strings from the command line in the form
  //
  //    m_name.num_boxes = 0
  //
  // If num_boxes > 0 we process lines
  //
  //    m_name.ref_box1_lo  = 0 0 0
  //    m_name.ref_box1_hi  = 1 1 1
  //    m_name.ref_box1_lvl = 1
  //    m_name.ref_box2_lo  = 0 0 0
  //    m_name.ref_box2_hi  = 1 1 1
  //    m_name.ref_box2_lvl = 3
  //
  // and so on.
  //
  // We then populate m_tagBoxes with RealBox (an axis-aligned box with physical coordinates). These
  // boxes can be used in tagCells(...) to prune regions in space where tagging is unnecessary. Users
  // are free to leave num_boxes=0 if they don't want to that.

  int numBoxes = 0;
  m_refBoxes.resize(0);

  // Get from input script.
  ParmParse pp(m_name.c_str());
  pp.query("num_ref_boxes", numBoxes);

  // Construct the refinement boxes.
  if (numBoxes > 0) {
    m_refBoxes.resize(numBoxes);

    for (int ibox = 0; ibox < numBoxes; ibox++) {
      const std::string str1 = "ref_box" + std::to_string(1 + ibox) + "_lo";
      const std::string str2 = "ref_box" + std::to_string(1 + ibox) + "_hi";
      const std::string str3 = "ref_box" + std::to_string(1 + ibox) + "_lvl";

      int          refLevel = 0;
      Vector<Real> cornerLo(SpaceDim);
      Vector<Real> cornerHi(SpaceDim);

      pp.getarr(str1.c_str(), cornerLo, 0, SpaceDim);
      pp.getarr(str2.c_str(), cornerHi, 0, SpaceDim);
      pp.get(str3.c_str(), refLevel);

      const RealVect c1 = RealVect(D_DECL(cornerLo[0], cornerLo[1], cornerLo[2]));
      const RealVect c2 = RealVect(D_DECL(cornerHi[0], cornerHi[1], cornerHi[2]));

      if (refLevel > 0) {
        m_refBoxes[ibox] = std::pair<RealBox, int>(RealBox(c1, c2), refLevel);
      }
    }
  }
}

void
CellTagger::parseBuffer()
{
  CH_TIME("CellTagger::parseBuffer()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseBuffer()" << endl;
  }

  // Just parse the refinement buffer. This
  ParmParse pp(m_name.c_str());

  pp.get("buffer", m_buffer);

  m_buffer = std::max(0, m_buffer);
}

void
CellTagger::parseVerbosity()
{
  CH_TIME("CellTagger::parseVerbosity()");

  ParmParse pp(m_name.c_str());

  pp.get("verbosity", m_verbosity);

  if (m_verbosity > 5) {
    pout() << m_name + "::parseVerbosity()" << endl;
  }
}

bool
CellTagger::insideTagBox(const RealVect a_pos) const
{
  CH_TIME("CellTagger::insideTagBox(RealVect)");
  if (m_verbosity > 5) {
    pout() << m_name + "::insideTagBox(RealVect)" << endl;
  }

  bool doThisRefine = (m_tagBoxes.size() > 0) ? false : true; // If we don't have any boxes, everything goes.

  for (int ibox = 0; ibox < m_tagBoxes.size(); ibox++) {
    const RealVect lo = m_tagBoxes[ibox].getLo();
    const RealVect hi = m_tagBoxes[ibox].getHi();

    // In this case a_pos is inside the tag box and we permit cell tagging.
    if (a_pos >= lo && a_pos <= hi) {
      doThisRefine = true;
    }
  }

  return doThisRefine;
}

int
CellTagger::getManualRefinementLevel(const RealVect a_pos) const
{
  CH_TIME("CellTagger::insideTagBox(RealVect)");
  if (m_verbosity > 5) {
    pout() << m_name + "::insideTagBox(RealVect)" << endl;
  }

  int refToThisLevel = -1;

  for (int i = 0; i < m_refBoxes.size(); i++) {
    const RealBox& refBox = m_refBoxes[i].first;
    const int&     refLvl = m_refBoxes[i].second;

    const RealVect lo = refBox.getLo();
    const RealVect hi = refBox.getHi();

    if (a_pos >= lo && a_pos <= hi) {
      refToThisLevel = std::max(refToThisLevel, refLvl);
    }
  }

  return refToThisLevel;
}

#include <CD_NamespaceFooter.H>
