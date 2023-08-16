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
                                                 const std::function<Real(const Real E, const RealVect x)>& a_alphaEff,
                                                 const phase::which_phase                                   a_phase)
{
  CH_TIME("StreamerInceptionTagger::StreamerInceptionTagger()");

  m_verbosity = -1;
  m_buffer    = 0;
  m_name      = "StreamerInceptionTagger";

  m_amr           = a_amrMesh;
  m_electricField = a_electricField;
  m_alphaEff      = a_alphaEff;
  m_phase         = a_phase;
  m_plot          = true;
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
  pp.get("ref_alpha", m_refAlpha);
  pp.get("max_voltage", m_maxVoltage);
  pp.get("plot", m_plot);
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
        if (ebisbox.isRegular(iv) && tracerFieldReg(iv, 0) >= m_refAlpha) {
          foundTags = 1;

          tags |= iv;
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (tracerField(vof, 0) >= m_refAlpha) {
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

  return m_plot ? 1 : 0;
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

  if (m_plot) {
    this->computeTracerField();

    DataOps::copy(a_output, m_tracerField, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_plotVariableNames.push_back(m_name + " alpha-criterion");

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

  const RealVect probLo = m_amr->getProbLo();

  // Compute alphaEff * dx on all grid levels.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real               dx  = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBCellFAB& electricField    = (*(*m_electricField)[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      EBCellFAB& tracerField    = (*m_tracerField[lvl])[dit()];
      FArrayBox& tracerFieldReg = tracerField.getFArrayBox();

      const EBISBox& ebisbox = electricField.getEBISBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect x  = probLo + dx * (0.5 * RealVect::Unit + RealVect(iv));
        const RealVect EE = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const Real     E  = EE.vectorLength();

        tracerFieldReg(iv, 0) = m_alphaEff(m_maxVoltage * E, x) * dx;
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect x  = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
        const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E  = EE.vectorLength();

        tracerField(vof, 0) = m_alphaEff(m_maxVoltage * E, x) * dx;
      };

      Box          cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  m_amr->conservativeAverage(m_tracerField, m_realm, m_phase);
}

#include <CD_NamespaceFooter.H>
