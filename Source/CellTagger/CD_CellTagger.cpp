/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
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

CellTagger::CellTagger(){
  CH_TIME("CellTagger::CellTagger");
  m_verbosity = 10;
  if(m_verbosity > 5){
    pout() << "CellTagger::CellTagger" << endl;
  }

  m_name  = "CellTagger";
}

CellTagger::~CellTagger(){

}

int CellTagger::getNumberOfPlotVariables(){
  return 0;
}

int CellTagger::getBuffer(){
  return m_buffer;
}

void CellTagger::parseRuntimeOptions(){

}

void CellTagger::parseBoxes(){
  
  ParmParse pp(m_name.c_str());

  int num_boxes = 0;
  pp.get("num_boxes", num_boxes);

  m_tagBoxes.resize(0);
  if(num_boxes > 0){
    m_tagBoxes.resize(num_boxes);

    const int ndigits = (int) log10((double) num_boxes) + 1;
      
    for (int ibox = 0; ibox < num_boxes; ibox++){
      char* cstr = new char[ndigits];
      sprintf(cstr, "%d", 1+ibox);

      std::string str1 = "box" + std::string(cstr) + "_lo";
      std::string str2 = "box" + std::string(cstr) + "_hi";

      Vector<Real> corner_lo(SpaceDim);
      Vector<Real> corner_hi(SpaceDim);

      pp.getarr(str1.c_str(), corner_lo, 0, SpaceDim);
      pp.getarr(str2.c_str(), corner_hi, 0, SpaceDim);

      const RealVect c1 = RealVect(D_DECL(corner_lo[0], corner_lo[1], corner_lo[2]));
      const RealVect c2 = RealVect(D_DECL(corner_hi[0], corner_hi[1], corner_hi[2]));

      m_tagBoxes[ibox] = RealBox(c1,c2);
	
      delete cstr;
    }
  }
}

void CellTagger::parseBuffer(){

  ParmParse pp(m_name.c_str());
  pp.get("buffer", m_buffer);
  m_buffer = Max(0, m_buffer);
}

void CellTagger::parseVerbosity(){
  CH_TIME("CellTagger::parseVerbosity");

  ParmParse pp(m_name.c_str());
  pp.get("verbosity", m_verbosity);
}

void CellTagger::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp){
  CH_TIME("CellTagger::writePlotData");
  if(m_verbosity > 3){
    pout() << "CellTagger::writePlotData" << endl;
  }
}

bool CellTagger::insideTagBox(const RealVect a_pos){
  bool do_this_refine = (m_tagBoxes.size() > 0) ? false : true;
  for (int ibox = 0; ibox < m_tagBoxes.size(); ibox++){
    const RealVect lo = m_tagBoxes[ibox].getLo();
    const RealVect hi = m_tagBoxes[ibox].getHi();

    if(a_pos >= lo && a_pos <= hi){
      do_this_refine = true;
    }
  }

  return do_this_refine;
}

#include <CD_NamespaceFooter.H>
