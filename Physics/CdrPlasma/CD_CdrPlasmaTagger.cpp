/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaTagger.cpp
  @brief  Implementation of CD_CdrPlasmaTagger.H
  @author Robert marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_CdrPlasmaTagger.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include <CD_BoxLoops.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaTagger::CdrPlasmaTagger()
{
  CH_TIME("CdrPlasmaTagger::CdrPlasmaTagger");
  if (m_verbosity > 5) {
    pout() << m_name + "::define" << endl;
  }

  // Default settings
  m_name  = "CdrPlasmaTagger";
  m_phase = phase::gas;
}

CdrPlasmaTagger::CdrPlasmaTagger(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
                                 const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
                                 const RefCountedPtr<AmrMesh>&               a_amr,
                                 const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
  : CdrPlasmaTagger()
{
  CH_TIME("CdrPlasmaTagger::CdrPlasmaTagger(RefCountedPtr<...> x 4)");
  if (m_verbosity > 5) {
    pout() << m_name + "::define(RefCountedPtr<...>x4)" << endl;
  }

  // Call the define function
  this->define(a_physics, a_timeStepper, a_amr, a_computationalGeometry);
}

CdrPlasmaTagger::~CdrPlasmaTagger()
{
  CH_TIME("CdrPlasmaTagger::~CdrPlasmaTagger");
  if (m_verbosity > 5) {
    pout() << m_name + "::~CdrPlasmaTagger" << endl;
  }
}

void
CdrPlasmaTagger::define(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
                        const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
                        const RefCountedPtr<AmrMesh>&               a_amr,
                        const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("CdrPlasmaTagger::define(RefCountedPtr<...>x4)");
  if (m_verbosity > 5) {
    pout() << m_name + "::define(RefCountedPtr<...>x4)" << endl;
  }

  m_physics               = a_physics;
  m_timeStepper           = a_timeStepper;
  m_amr                   = a_amr;
  m_realm                 = Realm::Primal;
  m_computationalGeometry = a_computationalGeometry;
  m_phase                 = phase::gas;
}

void
CdrPlasmaTagger::prePlot() const noexcept
{
  CH_TIME("CdrPlasmaTagger::prePlot()");
  if (m_verbosity > 5) {
    pout() << m_name + "::prePlot()" << endl;
  }

  this->computeTracers();
}

void
CdrPlasmaTagger::preRegrid() noexcept
{
  CH_TIME("CdrPlasmaTagger::preRegrid()");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid()" << endl;
  }

  for (int i = 0; i < m_numTracers; i++) {
    m_tracers[i].clear();
    m_gradTracers[i].clear();
  }
}

void
CdrPlasmaTagger::regrid()
{
  CH_TIME("CdrPlasmaTagger::regrid()");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid()" << endl;
  }

  // TLDR: This reallocates storage for the tracer fields.

  if (m_numTracers > 0) {
    m_tracers.resize(m_numTracers);
    m_gradTracers.resize(m_numTracers);

    // Allocate storage for scalar tracer and it's gradient.
    for (int i = 0; i < m_numTracers; i++) {
      m_amr->allocate(m_tracers[i], m_realm, m_phase, 1);
      m_amr->allocate(m_gradTracers[i], m_realm, m_phase, SpaceDim);
    }
  }
}

int
CdrPlasmaTagger::getNumberOfPlotVariables() const
{
  CH_TIME("CdrPlasmaTagger::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  return m_numTracers;
}

Vector<std::string>
CdrPlasmaTagger::getPlotVariableNames() const
{
  CH_TIME("CdrPlasmaTagger::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames" << endl;
  }

  Vector<std::string> plotVars;

  for (int i = 0; i < m_numTracers; i++) {
    plotVars.push_back("Tracer field-" + std::to_string(i));
  }

  return plotVars;
}

void
CdrPlasmaTagger::writePlotData(LevelData<EBCellFAB>& a_output,
                               int&                  a_icomp,
                               const std::string     a_outputRealm,
                               const int             a_level) const
{
  CH_TIME("CdrPlasmaTagger::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  // Go through the fields and add them to file.
  for (int i = 0; i < m_numTracers; i++) {
    const Interval srcInterv(0, 0);
    const Interval dstInterv(a_icomp, a_icomp);

    const EBAMRCellData& tagField = m_tracers[i];

    // Copy data to the ouput data holder. Covered data is bogus.
    m_amr->copyData(a_output, *tagField[a_level], a_level, a_outputRealm, tagField.getRealm(), dstInterv, srcInterv);

    DataOps::setCoveredValue(a_output, a_icomp, 0.0);

    a_icomp++;
  }
}

bool
CdrPlasmaTagger::tagCells(EBAMRTags& a_tags)
{
  CH_TIME("CdrPlasmaTagger::tagCells(EBAMRTags)");
  if (m_verbosity > 5) {
    pout() << m_name + "::tagCells(EBAMRTags)" << endl;
  }

  bool gotNewTags = false;

  // Lower left corner and current time.
  const RealVect probLo = m_amr->getProbLo();
  const Real     time   = m_timeStepper->getTime();

  // Determine how deep we should flag cells for refinement. We will never flag cells on AmrMesh's maximum AMR depth, so we restrict to that.
  const int finestLevel    = m_amr->getFinestLevel();
  const int maxDepth       = m_amr->getMaxAmrDepth();
  const int finestTagLevel = (finestLevel == maxDepth) ? maxDepth - 1 : finestLevel;

  if (m_numTracers > 0) {

    // Compute tracer fields -- note that this is implemented by subclasses.
    this->computeTracers();

    // Go through the AMR levels, being careful not to add cell tags on the maximum possible AMR depth (because grid level l+1 is generated from tags on level l).
    for (int lvl = 0; lvl <= finestTagLevel; lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];

      // Go through the grid patches.
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box      box     = dbl[dit()];
        const EBISBox& ebisbox = ebisl[dit()];

        // Create data holders that hold which cells were coarsened and which cells were refined
        DenseIntVectSet coarsenTags(box, false); // Cells that will be coarsened
        DenseIntVectSet refineTags(box, false);  // Cells that will be refined

        // Handle to tracer fields and their gradients on this grid patch.
        Vector<EBCellFAB*> tracers;
        Vector<EBCellFAB*> gtracers;

        for (int i = 0; i < m_numTracers; i++) {
          tracers.push_back(&((*m_tracers[i][lvl])[dit()]));
          gtracers.push_back(&((*m_gradTracers[i][lvl])[dit()]));
        }

        // Current cell flags. If we add a tag, the cell will be refined. Remove one, and it will be coarsened.
        DenseIntVectSet& tags = (*a_tags[lvl])[dit()];

        // Calls the patch versions which figure out which tags will be refined in each grid patch.
        this->refineCellsBox(refineTags, tracers, gtracers, lvl, dit(), box, ebisbox, time, dx, probLo);
        this->coarsenCellsBox(coarsenTags, tracers, gtracers, lvl, dit(), box, ebisbox, time, dx, probLo);

        // Check if we got any new tags, or we are just recycling old tags. If we did not get new tags then
        // we will ask the Driver to skip the regrid completely. Basically we will check if (current_tags + refined_tags - coarsenTags) == current_tags
        DenseIntVectSet cpy1 = tags;
        tags -= coarsenTags;
        tags |= refineTags;
        DenseIntVectSet cpy2 = tags;

        cpy2 -= cpy1; // = new tags minus old tags. If nonzero, we got some new tags.
        cpy1 -= tags; // = old_tags minus new tags. If nonzero, we got some new tags
        if (cpy1.numPts() != 0 || cpy2.numPts() != 0) {
          gotNewTags = true;
        }

        // No tags allowed outside the current grid patch.
        tags &= box;
      }
    }
  }

  // Some ranks may have gotten new tags while others have not. This little code snippet
  // sets gotNewTags = true for all ranks if any rank probLoally had gotNewTags = true
  int hasTags = gotNewTags ? 1 : 0;
  hasTags     = ParallelOps::max(hasTags);

  if (hasTags > 0) {
    gotNewTags = true;
  }

  return gotNewTags;
}

void
CdrPlasmaTagger::refineCellsBox(DenseIntVectSet&          a_refinedCells,
                                const Vector<EBCellFAB*>& a_tracers,
                                const Vector<EBCellFAB*>& a_gradTracers,
                                const int                 a_lvl,
                                const DataIndex           a_dit,
                                const Box                 a_box,
                                const EBISBox&            a_ebisbox,
                                const Real                a_time,
                                const Real                a_dx,
                                const RealVect            a_probLo)
{
  CH_TIME("CdrPlasmaTagger::refineCellsBox(...)");
  if (m_verbosity > 5) {
    pout() << m_name + "::refineCellsBox(...)" << endl;
  }

  // Get a handle to the single-valued data.
  Vector<FArrayBox*> tracersReg;
  Vector<FArrayBox*> gradientsReg;

  for (int i = 0; i < m_numTracers; i++) {
    tracersReg.push_back(&(a_tracers[i]->getFArrayBox()));
    gradientsReg.push_back(&(a_gradTracers[i]->getFArrayBox()));
  }

  // Cell-wise tracer fields. Used in the kernels.
  Vector<Real>     tr(m_numTracers);
  Vector<RealVect> gt(m_numTracers);

  // Regular kernel.
  auto regularKernel = [&](const IntVect& iv) -> void {
    const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;

    // If position is inside any of the tagging boxes, we can refine
    if (this->insideTagBox(pos) && a_ebisbox.isRegular(iv)) {

      // Get the tracer fields and the gradients of them.
      for (int i = 0; i < m_numTracers; i++) {
        tr[i] = (*tracersReg[i])(iv, 0);
        gt[i] = RealVect(D_DECL((*gradientsReg[i])(iv, 0), (*gradientsReg[i])(iv, 1), (*gradientsReg[i])(iv, 2)));
      }

      // Call the per-cell-refinement method.
      const bool refine = this->refineCell(pos, a_time, a_dx, a_lvl, tr, gt);

      // If we refine, grow with buffer and increment to a_refinedCells
      if (refine) {
        a_refinedCells |= iv;
      }
    }
  };

  // Irregular kernel
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect pos = a_probLo + Location::position(Location::Cell::Center, vof, a_ebisbox, a_dx);

    // If position is inside any of the tagging boxes, we can refine
    if (this->insideTagBox(pos)) {

      // Get the tracer fields and the gradients of them.
      for (int i = 0; i < m_numTracers; i++) {
        tr[i] = (*a_tracers[i])(vof, 0);
        gt[i] = RealVect(D_DECL((*a_gradTracers[i])(vof, 0), (*a_gradTracers[i])(vof, 1), (*a_gradTracers[i])(vof, 2)));
      }

      // Call the per-cell-refinement method.
      const bool refine = this->refineCell(pos, a_time, a_dx, a_lvl, tr, gt);

      if (refine) {
        a_refinedCells |= vof.gridIndex();
      }
    }
  };

  // Irregular kernel region. Should probably be stored in its own VoFIterator but I'm lazy so let's define it right here.
  VoFIterator vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];

  // Execute the kernels.
  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

void
CdrPlasmaTagger::coarsenCellsBox(DenseIntVectSet&          a_coarsenedCells,
                                 const Vector<EBCellFAB*>& a_tracers,
                                 const Vector<EBCellFAB*>& a_gradTracers,
                                 const int                 a_lvl,
                                 const DataIndex           a_dit,
                                 const Box                 a_box,
                                 const EBISBox&            a_ebisbox,
                                 const Real                a_time,
                                 const Real                a_dx,
                                 const RealVect            a_probLo)
{
  CH_TIME("CdrPlasmaTagger::coarsenCellsBox(...)");
  if (m_verbosity > 5) {
    pout() << m_name + "::coarsenCellsBox(...)" << endl;
  }

  // Get a handle to the single-valued data.
  Vector<FArrayBox*> tracersReg;
  Vector<FArrayBox*> gradientsReg;

  for (int i = 0; i < m_numTracers; i++) {
    tracersReg.push_back(&(a_tracers[i]->getFArrayBox()));
    gradientsReg.push_back(&(a_gradTracers[i]->getFArrayBox()));
  }

  // Tracer field and gradient of tracer fields. Used in the kernels.
  Vector<Real>     tr(m_numTracers);
  Vector<RealVect> gt(m_numTracers);

  // Regular box loop
  auto regularKernel = [&](const IntVect& iv) -> void {
    const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;

    // If position is inside any of the tagging boxes, we can refine
    const bool isPointInside = this->insideTagBox(pos);
    if (!isPointInside) {
      a_coarsenedCells |= iv; // Always coarsen outside the refinement boxes.
    }
    else if (isPointInside && a_ebisbox.isRegular(iv)) {

      // Reconstruct the tracer fields and tracer field gradients
      for (int i = 0; i < m_numTracers; i++) {
        tr[i] = (*tracersReg[i])(iv, 0);
        gt[i] = RealVect(D_DECL((*gradientsReg[i])(iv, 0), (*gradientsReg[i])(iv, 1), (*gradientsReg[i])(iv, 2)));
      }

      // Call the per-cell-coarsening method and see if we should coarsen this cell.
      const bool coarsen = this->coarsenCell(pos, a_time, a_dx, a_lvl, tr, gt);

      // If we refine, grow with buffer and increment to a_refinedCells
      if (coarsen) {
        a_coarsenedCells |= iv;
      }
    }
  };

  // Irregular box loop
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect pos = a_probLo + Location::position(Location::Cell::Center, vof, a_ebisbox, a_dx);

    // If position is inside any of the tagging boxes, we can refine
    const bool isPointInside = this->insideTagBox(pos);
    if (!isPointInside) {
      a_coarsenedCells |= vof.gridIndex();
    }
    else
      // Reconstruct the tracer fields and gradients again.
      for (int i = 0; i < m_numTracers; i++) {
        tr[i] = (*a_tracers[i])(vof, 0);
        gt[i] = RealVect(D_DECL((*a_gradTracers[i])(vof, 0), (*a_gradTracers[i])(vof, 1), (*a_gradTracers[i])(vof, 2)));
      }

    // Call the per-cell refinement method.
    const bool coarsen = this->coarsenCell(pos, a_time, a_dx, a_lvl, tr, gt);

    if (coarsen) {
      a_coarsenedCells |= vof.gridIndex();
    }
  };

  // Kernel region for the irregular kernel.
  VoFIterator vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];

  // Execute the kernels
  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

#include <CD_NamespaceFooter.H>
