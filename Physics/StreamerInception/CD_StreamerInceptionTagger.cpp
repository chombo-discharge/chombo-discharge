/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionTagger.cpp
  @brief  Implementation of CD_StreamerInceptionTagger.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_StreamerInceptionTagger.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::StreamerInception;

StreamerInceptionTagger::StreamerInceptionTagger(const RefCountedPtr<AmrMesh>& a_amrMesh,
						 const EBAMRCellData* const a_electricField,
						 const phase::which_phase a_phase) {

  m_amr = a_amrMesh;
  m_electricField = a_electricField;
  m_phase = a_phase;
}


StreamerInceptionTagger::~StreamerInceptionTagger() {

}

void
StreamerInceptionTagger::parseOptions() {
}

void
StreamerInceptionTagger::regrid() {

}

bool
StreamerInceptionTagger::tagCells(EBAMRTags& a_tags) {

}

#include <CD_NamespaceFooter.H>
