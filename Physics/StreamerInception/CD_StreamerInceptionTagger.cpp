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
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::StreamerInception;

StreamerInceptionTagger::StreamerInceptionTagger(const RefCountedPtr<AmrMesh>& a_amrMesh,
                                                 const EBAMRCellData* const    a_electricField,
                                                 const phase::which_phase      a_phase)
{
  CH_TIME("StreamerInceptionTagger::StreamerInceptionTagger()");

  m_verbosity = -1;
  m_buffer    = 0;
  m_name      = "StreamerInceptionTagger";

  m_amr           = a_amrMesh;
  m_electricField = a_electricField;
  m_phase         = a_phase;
  m_plotTracer    = true;
}

StreamerInceptionTagger::~StreamerInceptionTagger()
{
  CH_TIME("StreamerInceptionTagger::~StreamerInceptionTagger()");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::~StreamerInceptionTagger" << endl;
  }
}

void
StreamerInceptionTagger::parseOptions()
{
  CH_TIME("StreamerInceptionTagger::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::parseOptions()" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("verbosity", m_verbosity);
  pp.get("buffer", m_buffer);
  pp.get("ref_curv", m_refCurv);
  pp.get("plot_curv", m_plotTracer);
}

void
StreamerInceptionTagger::regrid()
{
  CH_TIME("StreamerInceptionTagger::regrid()");
  if (m_verbosity > 5) {
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
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::tagCells()" << endl;
  }

  this->computeTracerField();

  int foundTags = -1;

  const int finestLevel    = m_amr->getFinestLevel();
  const int maxLevel       = m_amr->getMaxAmrDepth();
  const int finestTagLevel = (finestLevel == maxLevel) ? maxLevel - 1 : finestLevel;

  for (int lvl = 0; lvl <= finestTagLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      DenseIntVectSet& tags = (*a_tags[lvl])[dit()];
      tags.makeEmptyBits();

      const EBCellFAB& tracerField    = (*m_tracerField[lvl])[dit()];
      const FArrayBox& tracerFieldReg = tracerField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv) && tracerFieldReg(iv, 0) >= m_refCurv) {
          foundTags = 1;

          tags |= iv;
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (tracerField(vof, 0) >= m_refCurv) {
          foundTags = 1;

          tags |= vof.gridIndex();
        }
      };

      // Run the kernels.
      Box         cellBox = dbl[dit()];
      VoFIterator vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      // Execute the kernels.
      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  return (ParallelOps::max(foundTags) > 0) ? true : false;
}

int
StreamerInceptionTagger::getNumberOfPlotVariables() const
{
  CH_TIME("StreamerInceptionTagger::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::getNumberOfPlotVariables()" << endl;
  }

  return m_plotTracer ? 1 : 0;
}

void
StreamerInceptionTagger::writePlotData(EBAMRCellData&       a_output,
                                       Vector<std::string>& a_plotVariableNames,
                                       int&                 a_icomp) const
{
  CH_TIME("StreamerInceptionTagger::writePlotData()");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::writePlotData()" << endl;
  }

  if (m_plotTracer) {
    this->computeTracerField();

    DataOps::copy(a_output, m_tracerField, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_plotVariableNames.push_back("Tagging curvature field");

    a_icomp++;
  }
}

void
StreamerInceptionTagger::computeTracerField() const noexcept
{
  CH_TIME("StreamerInceptionTagger::computeTracerField()");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionTagger::computeTracerField()" << endl;
  }

  CH_assert(m_electricField != nullptr);

  EBAMRCellData magnE;
  EBAMRCellData gradE;

  m_amr->allocate(magnE, m_realm, m_phase, 1);
  m_amr->allocate(gradE, m_realm, m_phase, SpaceDim);

  // Compute |E|
  DataOps::vectorLength(magnE, *m_electricField);

  m_amr->averageDown(magnE, m_realm, m_phase);
  m_amr->interpGhostMG(magnE, m_realm, m_phase);

  // Compute |grad(|E|)| onto m_tracerField
  m_amr->computeGradient(gradE, magnE, m_realm, m_phase);

  DataOps::vectorLength(m_tracerField, gradE);

  // Compute |grad(|E|)| * dx / E
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real dx = m_amr->getDx()[lvl];

    DataOps::scale(*m_tracerField[lvl], dx);
  }

  DataOps::divideFallback(m_tracerField, magnE, 0.0);
}

#include <CD_NamespaceFooter.H>
