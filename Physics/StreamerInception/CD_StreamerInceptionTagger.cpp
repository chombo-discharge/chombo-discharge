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
#include <ParmParse.H>

// Our includes
#include <CD_StreamerInceptionTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::StreamerInception;

StreamerInceptionTagger::StreamerInceptionTagger(const RefCountedPtr<AmrMesh>& a_amrMesh,
                                                 const EBAMRCellData* const    a_electricField,
                                                 const phase::which_phase      a_phase)
{
  CH_TIME("StreamerInceptionTagger::StreamerInceptionTagger()");

  m_verbosity = -1;
  m_buffer    = 0;
  m_name = "StreamerInceptionTagger";

  m_amr           = a_amrMesh;
  m_electricField = a_electricField;
  m_phase         = a_phase;
}

StreamerInceptionTagger::~StreamerInceptionTagger() {
  CH_TIME("StreamerInceptionTagger::~StreamerInceptionTagger()");
  if(m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::~StreamerInceptionTagger" << endl;
  }  
}

void
StreamerInceptionTagger::parseOptions() {
  CH_TIME("StreamerInceptionTagger::parseOptions()");
  if(m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::parseOptions()" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("verbosity", m_verbosity);  
  pp.get("buffer", m_buffer);
  pp.get("ref_curv", m_refCurv);
}

void
StreamerInceptionTagger::regrid()
{
  CH_TIME("StreamerInceptionTagger::regrid()");
  if(m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::regrid()" << endl;
  }
  
  CH_assert(m_electricField != nullptr);

  m_realm = m_electricField->getRealm();
  
  m_amr->allocate(m_tracerField, m_realm, m_phase, 1);
  DataOps::setValue(m_tracerField, 0.0);
}

bool
StreamerInceptionTagger::tagCells(EBAMRTags& a_tags)
{
  CH_TIME("StreamerInceptionTagger::tagCells()");
  if(m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::tagCells()" << endl;
  }

  CH_assert(m_electricField != nullptr);

  EBAMRCellData magnE;
  EBAMRCellData gradE;

  m_amr->allocate(magnE, m_realm, m_phase, 1);
  m_amr->allocate(gradE, m_realm, m_phase, SpaceDim);

  // Compute |E|
  DataOps::vectorLength(magnE, *m_electricField);

  m_amr->averageDown(magnE, m_realm, m_phase);
  m_amr->interpGhost(magnE, m_realm, m_phase);

  // Compute |grad(|E|)| onto m_tracerField
  m_amr->computeGradient(gradE, magnE, m_realm, m_phase);

  DataOps::vectorLength(m_tracerField, gradE);

  // Compute |grad(|E|)| * dx / E
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real dx = m_amr->getDx()[lvl];

    DataOps::scale(*m_tracerField[lvl], dx);
  }

  DataOps::divideFallback(m_tracerField, magnE, 0.0);

  // Flag cells for refinement.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++ ) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      DenseIntVectSet& tags = (*a_tags[lvl])[dit()];

      const EBCellFAB& tracerField = (*m_tracerField[lvl])[dit()];
      const FArrayBox& tracerFieldReg = tracerField.getFArrayBox();



      auto regularKernel = [&](const IntVect& iv) -> void {
	if(tracerFieldReg(iv, 0) >= m_refCurv) {
	  tags |= iv;
	}
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
	if(tracerField(vof, 0) >= m_refCurv) {
	  tags |= vof.gridIndex();
	}
      };

      
      // Run the kernels.
      Box cellBox = dbl[dit()];
      VoFIterator vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      // Execute the kernels.
      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);	
    }
  }
}

#include <CD_NamespaceFooter.H>
