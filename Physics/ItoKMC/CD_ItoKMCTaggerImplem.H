/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCTaggerImplem.H
  @brief  Implementation of CD_ItoKMCTagger.H
  @author Robert marskar
*/

#ifndef CD_ItoKMCTaggerImplem_H
#define CD_ItoKMCTaggerImplem_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoKMCTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

template <typename S>
ItoKMCTagger<S>::ItoKMCTagger()
{
  CH_TIME("ItoKMCTagger::ItoKMCTagger");
  m_verbosity = -1;
  if (m_verbosity > 5) {
    pout() << "ItoKMCTagger::ItoKMCTagger" << endl;
  }

  m_name      = "ItoKMCTagger";
  m_phase     = phase::gas;
  m_realm     = Realm::Primal;
  m_isDefined = false;
  m_plotData  = false;
}

template <typename S>
ItoKMCTagger<S>::ItoKMCTagger(const RefCountedPtr<ItoKMCPhysics>& a_physics,
                              const RefCountedPtr<S>&             a_timeStepper,
                              const RefCountedPtr<AmrMesh>&       a_amr) noexcept
  : ItoKMCTagger<S>()
{
  this->define(a_physics, a_timeStepper, a_amr);
}

template <typename S>
ItoKMCTagger<S>::~ItoKMCTagger() noexcept
{}

template <typename S>
void
ItoKMCTagger<S>::define(const RefCountedPtr<ItoKMCPhysics>& a_physics,
                        const RefCountedPtr<S>&             a_timeStepper,
                        const RefCountedPtr<AmrMesh>&       a_amr) noexcept

{
  CH_TIME("ItoKMCTagger::define");
  if (m_verbosity > 5) {
    pout() << m_name + "::define" << endl;
  }

  m_physics     = a_physics;
  m_timeStepper = a_timeStepper;
  m_amr         = a_amr;
  m_isDefined   = true;
}

template <typename S>
void
ItoKMCTagger<S>::regrid() noexcept
{
  CH_TIME("ItoKMCTagger::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  CH_assert(m_isDefined);

  if (m_numTagFields > 0) {

    m_tagFields.resize(m_numTagFields);
    m_gradTagFields.resize(m_numTagFields);

    for (int i = 0; i < m_numTagFields; i++) {
      m_amr->allocate(m_tagFields[i], m_realm, m_phase, 1);
      m_amr->allocate(m_gradTagFields[i], m_realm, m_phase, SpaceDim);
    }
  }
}

template <typename S>
void
ItoKMCTagger<S>::setPhase(const phase::which_phase a_phase) noexcept
{
  CH_TIME("ItoKMCTagger::setPhase");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

template <typename S>
int
ItoKMCTagger<S>::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("ItoKMCTagger::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  return m_plotData ? m_numTagFields : 0;
}

template <typename S>
Vector<EBAMRCellData>&
ItoKMCTagger<S>::getTagFields() noexcept
{
  CH_TIME("ItoKMCTagger::getTagFields");
  if (m_verbosity > 5) {
    pout() << m_name + "::getTagFields" << endl;
  }

  CH_assert(m_isDefined);

  return m_tagFields;
}

template <typename S>
void
ItoKMCTagger<S>::writePlotData(EBAMRCellData&       a_output,
                               Vector<std::string>& a_plotVariableNames,
                               int&                 a_icomp) const noexcept
{
  CH_TIME("ItoKMCTagger::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  CH_assert(m_isDefined);

  if (m_plotData) {
    this->computeTagFields();

    for (int i = 0; i < m_numTagFields; i++) {
      const EBAMRCellData& tagField = m_tagFields[i];

      const Interval srcInterv(0, 0);
      const Interval dstInterv(a_icomp, a_icomp);

      a_output.copy(srcInterv, tagField, dstInterv);

      DataOps::setCoveredValue(a_output, a_icomp, 0.0);

      a_plotVariableNames.push_back("Tag field-" + std::to_string(i));

      a_icomp++;
    }
  }
}

template <typename S>
bool
ItoKMCTagger<S>::tagCells(EBAMRTags& a_tags) noexcept
{
  CH_TIME("ItoKMCTagger::tagCells");
  if (m_verbosity > 5) {
    pout() << m_name + "::tagCells" << endl;
  }

  CH_assert(m_isDefined);

  bool gotNewTags = false;

  const RealVect probLo         = m_amr->getProbLo();
  const Real     curTime        = m_timeStepper->getTime();
  const int      finestLevel    = m_amr->getFinestLevel();
  const int      maxAmrDepth    = m_amr->getMaxAmrDepth();
  const int      finestTagLevel = (finestLevel == maxAmrDepth) ? maxAmrDepth - 1 : finestLevel;

  if (m_numTagFields > 0) {
    this->computeTagFields();

    for (int lvl = 0; lvl <= finestTagLevel; lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box      box     = dbl.get(dit());
        const EBISBox& ebisbox = ebisl[dit()];

        const IntVectSet irregCells = ebisbox.getIrregIVS(box);
        const IntVectSet backupTags = IntVectSet((*a_tags[lvl])[dit()]);

        DenseIntVectSet coarsenCells(box, false);
        DenseIntVectSet refinedCells(box, false);

        Vector<EBCellFAB*> tagFields;
        Vector<EBCellFAB*> gradTagFields;

        for (int i = 0; i < m_numTagFields; i++) {
          tagFields.push_back(&((*m_tagFields[i][lvl])[dit()]));
          gradTagFields.push_back(&((*m_gradTagFields[i][lvl])[dit()]));
        }

        this->tagCellsBox(refinedCells,
                          coarsenCells,
                          tagFields,
                          gradTagFields,
                          lvl,
                          dit(),
                          box,
                          ebisbox,
                          curTime,
                          dx,
                          probLo);

        // Check if we got any new tags, or we are just recycling old tags.
        // Basically we will check if (current_tags + refined_tags - coarsenCells) == current_tags
        DenseIntVectSet& tags = (*a_tags[lvl])[dit()];

        DenseIntVectSet cpy1 = tags;
        tags -= coarsenCells;
        tags |= refinedCells;
        DenseIntVectSet cpy2 = tags;

        cpy2 -= cpy1; // = new tags minus old tags. If nonzero, we got some new tags.
        cpy1 -= tags; // = old_tags minus new tags. If nonzero, we got some new tags
        if (cpy1.numPts() != 0 || cpy2.numPts() != 0) {
          gotNewTags = true;
        }
      }
    }
  }

  // Some ranks may have gotten new tags while others have not. This little code snippet
  // sets gotNewTags = true for all ranks if any rank probLoally had gotNewTags = true
  const int foundTags = gotNewTags ? 1 : 0;

  return (ParallelOps::max(foundTags) > 0);
}

template <typename S>
void
ItoKMCTagger<S>::tagCellsBox(DenseIntVectSet&          a_refinedCells,
                             DenseIntVectSet&          a_coarsenedCells,
                             const Vector<EBCellFAB*>& a_tagFields,
                             const Vector<EBCellFAB*>& a_gradTagFields,
                             const int                 a_lvl,
                             const DataIndex           a_dit,
                             const Box                 a_box,
                             const EBISBox&            a_ebisbox,
                             const Real                a_time,
                             const Real                a_dx,
                             const RealVect            a_probLo) const noexcept
{
  CH_TIME("ItoKMCTagger::refineCellsBox");
  if (m_verbosity > 5) {
    pout() << "ItoKMCTagger::refineCellsBox" << endl;
  }

  CH_assert(m_isDefined);

  Vector<FArrayBox*> tagFieldsReg;
  Vector<FArrayBox*> gradTagFieldsReg;

  for (int i = 0; i < m_numTagFields; i++) {
    tagFieldsReg.push_back(&(a_tagFields[i]->getFArrayBox()));
    gradTagFieldsReg.push_back(&(a_gradTagFields[i]->getFArrayBox()));
  }

  // Regular kernel.
  auto regularKernel = [&](const IntVect& iv) -> void {
    const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;

    if (this->insideTagBox(pos) && a_ebisbox.isRegular(iv)) {

      Vector<Real>     tr(m_numTagFields);
      Vector<RealVect> gt(m_numTagFields);

      for (int i = 0; i < m_numTagFields; i++) {
        tr[i] = (*tagFieldsReg[i])(iv, 0);
        gt[i] = RealVect(
          D_DECL((*gradTagFieldsReg[i])(iv, 0), (*gradTagFieldsReg[i])(iv, 1), (*gradTagFieldsReg[i])(iv, 2)));
      }

      if (this->refineCell(pos, a_time, a_dx, a_lvl, tr, gt)) {
        a_refinedCells |= iv;
      }

      if (this->coarsenCell(pos, a_time, a_dx, a_lvl, tr, gt)) {
        a_coarsenedCells |= iv;
      }
    }
  };

  // Irregular box loop
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect pos = a_probLo + Location::position(Location::Cell::Center, vof, a_ebisbox, a_dx);

    if (this->insideTagBox(pos)) {

      Vector<Real>     tr(m_numTagFields);
      Vector<RealVect> gt(m_numTagFields);

      for (int i = 0; i < m_numTagFields; i++) {
        tr[i] = (*a_tagFields[i])(vof, 0);
        gt[i] = RealVect(
          D_DECL((*a_gradTagFields[i])(vof, 0), (*a_gradTagFields[i])(vof, 1), (*a_gradTagFields[i])(vof, 2)));
      }

      if (this->refineCell(pos, a_time, a_dx, a_lvl, tr, gt)) {
        a_refinedCells |= vof.gridIndex();
      }

      if (this->coarsenCell(pos, a_time, a_dx, a_lvl, tr, gt)) {
        a_coarsenedCells |= vof.gridIndex();
      }
    }
  };

  // Run the kernels.
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];

  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

#include <CD_NamespaceFooter.H>

#endif
