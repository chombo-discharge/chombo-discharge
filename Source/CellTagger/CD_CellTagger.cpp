/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CellTagger.cpp
  @brief  Implementation of CD_CellTagger.H
  @author Robert marskar
*/

// Chombo includes
#include <EBArith.H>
#include <ParmParse.H>

// Our includes
#include <CD_DataOps.H>
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

int
CellTagger::getNumberOfPlotVariables() const
{
  CH_TIME("CellTagger::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  return 0;
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
CellTagger::parseBoxes()
{
  CH_TIME("CellTagger::parseBoxes()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseBoxes()" << endl;
  }

  // TLDR: This code parses strings from the command line in the form
  //
  //    m_name.num_boxes = 0
  //
  // If num_boxes > 0 we process lines
  //
  //    m_name.box1_lo = 0 0 0
  //    m_name.box1_hi = 1 1 1
  //    m_name.box2_lo = 0 0 0
  //    m_name.box2_hi = 1 1
  //
  // We then populate m_tagBoxes with RealBox (an axis-aligned box with physical coordinates). These
  // boxes can be used in tagCells(...) to prune regions in space where tagging is unnecessary. Users
  // are free to leave num_boxes=0 if they don't want to that.

  ParmParse pp(m_name.c_str());

  int numBoxes;
  pp.get("num_boxes", numBoxes);

  m_tagBoxes.resize(0);

  if (numBoxes > 0) {
    m_tagBoxes.resize(numBoxes);

    const int ndigits = (int)log10(1.0 * numBoxes) + 1;

    for (int ibox = 0; ibox < numBoxes; ibox++) {
      char* cstr = new char[ndigits];
      sprintf(cstr, "%d", 1 + ibox);

      const std::string str1 = "box" + std::string(cstr) + "_lo";
      const std::string str2 = "box" + std::string(cstr) + "_hi";

      Vector<Real> cornerLo(SpaceDim);
      Vector<Real> cornerHi(SpaceDim);

      pp.getarr(str1.c_str(), cornerLo, 0, SpaceDim);
      pp.getarr(str2.c_str(), cornerHi, 0, SpaceDim);

      const RealVect c1 = RealVect(D_DECL(cornerLo[0], cornerLo[1], cornerLo[2]));
      const RealVect c2 = RealVect(D_DECL(cornerHi[0], cornerHi[1], cornerHi[2]));

      m_tagBoxes[ibox] = RealBox(c1, c2);

      delete[] cstr;
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

void
CellTagger::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const
{
  CH_TIME("CellTagger::writePlotData(EBAMRCellData, Vector<std::string>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData(EBAMRCellData, Vector<std::string>, int)" << endl;
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

#include <CD_NamespaceFooter.H>
