/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD-BrownianWalkerTagger.cpp
  @brief  Implementation of CD_BrownianWalkerTagger.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_BrownianWalkerTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>


using namespace Physics::BrownianWalker;

BrownianWalkerTagger::BrownianWalkerTagger(RefCountedPtr<ItoSolver>& a_solver,
					       RefCountedPtr<AmrMesh>&   a_amr){
  m_solver    = a_solver;
  m_amr       = a_amr;
  m_name      = "BrownianWalker";
  m_verbosity = -1;
}

BrownianWalkerTagger::~BrownianWalkerTagger(){

}

void BrownianWalkerTagger::regrid(){

}

void BrownianWalkerTagger::parseOptions(){
  ParmParse pp(m_name.c_str());
  pp.get("refine_magn", m_refi_magn);
  
  parseBuffer(); // Derived from CellTagger
}

bool BrownianWalkerTagger::tagCells(EBAMRTags& a_tags){
  return true;
}

#include <CD_NamespaceFooter.H>
