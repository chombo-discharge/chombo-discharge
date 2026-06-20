/*
 * SPDX-FileCopyrightText: 2022-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_DischargeInceptionTagger.cpp
   @brief  Implementation of CD_DischargeInceptionTagger.H
   @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_DischargeInceptionTagger.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::DischargeInception;

DischargeInceptionTagger::DischargeInceptionTagger(
  const RefCountedPtr<AmrMesh>&                              a_amrMesh,
  const EBAMRCellData* const                                 a_electricField,
  const std::function<Real(const Real E, const RealVect x)>& a_alphaEff,
  const phase::which_phase                                   a_phase)
  : m_amr(a_amrMesh), m_electricField(a_electricField), m_alphaEff(a_alphaEff), m_phase(a_phase), m_plot(true)
{
  CH_TIME("DischargeInceptionTagger::DischargeInceptionTagger()");

  m_verbosity = -1;
  m_buffer    = 0;
  m_name      = "DischargeInceptionTagger";
}

DischargeInceptionTagger::~DischargeInceptionTagger()
{
  CH_TIME("DischargeInceptionTagger::~DischargeInceptionTagger()");
  if (m_verbosity > 5) {
    pout() << "DischargeInceptionTagger::~DischargeInceptionTagger" << endl;
  }
}

void
DischargeInceptionTagger::parseOptions()
{
  CH_TIME("DischargeInceptionTagger::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "DischargeInceptionTagger::parseOptions()" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("verbosity", m_verbosity);
  pp.get("buffer", m_buffer);
  pp.get("ref_alpha", m_refAlpha);
  pp.get("max_voltage", m_maxVoltage);
  pp.get("plot", m_plot);

  this->parseRefinementBoxes();
}

void
DischargeInceptionTagger::regrid()
{
  CH_TIME("DischargeInceptionTagger::regrid()");
  if (m_verbosity > 5) {
    pout() << "DischargeInceptionTagger::regrid()" << endl;
  }

  CH_assert(m_electricField != nullptr);

  m_realm = m_electricField->getRealm();

  m_amr->allocate(m_tracerField, m_realm, m_phase, 1);

  DataOps::setValue(m_tracerField, 0.0);
}

bool
DischargeInceptionTagger::tagCells(EBAMRTags& a_tags)
{
  CH_TIME("DischargeInceptionTagger::tagCells()");
  if (m_verbosity > 5) {
    pout() << "DischargeInceptionTagger::tagCells()" << endl;
  }

  this->computeTracerField();

  int foundTags = 0;

  const int finestLevel    = m_amr->getFinestLevel();
  const int maxLevel       = m_amr->getMaxAmrDepth();
  const int finestTagLevel = (finestLevel == maxLevel) ? maxLevel - 1 : finestLevel;

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= finestTagLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : foundTags)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const EBISBox& ebisbox = ebisl[din];

      DenseIntVectSet& tags = (*a_tags[lvl])[din];

      const EBCellFAB& tracerField    = (*m_tracerField[lvl])[din];
      const FArrayBox& tracerFieldReg = tracerField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv) && tracerFieldReg(iv, 0) >= m_refAlpha) {
          foundTags = 1;

          tags |= iv;
        }

        // Check for manual refinement
        const RealVect physPos = probLo + (iv + 0.5 * RealVect::Unit) * dx;

        if (this->getManualRefinementLevel(physPos) > lvl) {
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
      Box          cellBox = dbl[din];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

      // Execute the kernels. Not vectorizable: data-dependent DenseIntVectSet insertion (tags |= iv) +
      // out-of-line getManualRefinementLevel call + control flow. One-time tagging (per regrid).
      BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  return (ParallelOps::max(foundTags) > 0) ? true : false;
}

int
DischargeInceptionTagger::getNumberOfPlotVariables() const
{
  CH_TIME("DischargeInceptionTagger::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << "DischargeInceptionTagger::getNumberOfPlotVariables()" << endl;
  }

  return m_plot ? 1 : 0;
}

void
DischargeInceptionTagger::writePlotData(EBAMRCellData&       a_output,
                                        Vector<std::string>& a_plotVariableNames,
                                        int&                 a_icomp) const
{
  CH_TIME("DischargeInceptionTagger::writePlotData()");
  if (m_verbosity > 5) {
    pout() << "DischargeInceptionTagger::writePlotData()" << endl;
  }

  if (m_plot) {
    this->computeTracerField();

    DataOps::copy(a_output, m_tracerField, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_plotVariableNames.push_back(m_name + " alpha-criterion");

    a_icomp++;
  }
}

void
DischargeInceptionTagger::computeTracerField() const noexcept
{
  CH_TIME("DischargeInceptionTagger::computeTracerField()");
  if (m_verbosity > 5) {
    pout() << "DischargeInceptionTagger::computeTracerField()" << endl;
  }

  CH_assert(m_electricField != nullptr);

  const RealVect probLo = m_amr->getProbLo();

  // Compute alphaEff * dx on all grid levels.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = m_amr->getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const EBCellFAB& electricField    = (*(*m_electricField)[lvl])[din];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      EBCellFAB& tracerField    = (*m_tracerField[lvl])[din];
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

      Box          cellBox = dbl[din];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

      // Not vectorizable: m_alphaEff (effective Townsend coefficient) is a std::function call per cell.
      BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  m_amr->conservativeAverage(m_tracerField, m_realm, m_phase);
}

#include <CD_NamespaceFooter.H>
