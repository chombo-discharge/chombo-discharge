/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaTagger.cpp
  @brief  Implementation of CD_ItoPlasmaTagger.H
  @author Robert marskar
*/

// Chombo includes
#include <EBArith.H>
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaTagger.H>
#include <CD_ItoIterator.H>
#include <CD_RtIterator.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaTagger::ItoPlasmaTagger()
{
  CH_TIME("ItoPlasmaTagger::ItoPlasmaTagger");
  m_verbosity = 10;
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaTagger::ItoPlasmaTagger" << endl;
  }

  m_name  = "ItoPlasmaTagger";
  m_phase = phase::gas;
  m_realm = Realm::Primal;
}

ItoPlasmaTagger::ItoPlasmaTagger(const RefCountedPtr<ItoPlasmaPhysics>&      a_physics,
                                 const RefCountedPtr<ItoPlasmaStepper>&      a_timeStepper,
                                 const RefCountedPtr<AmrMesh>&               a_amr,
                                 const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
  : ItoPlasmaTagger()
{
  this->define(a_physics, a_timeStepper, a_amr, a_computationalGeometry);
}

ItoPlasmaTagger::~ItoPlasmaTagger() {}

void
ItoPlasmaTagger::define(const RefCountedPtr<ItoPlasmaPhysics>&      a_physics,
                        const RefCountedPtr<ItoPlasmaStepper>&      a_timeStepper,
                        const RefCountedPtr<AmrMesh>&               a_amr,
                        const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("ItoPlasmaTagger::define");
  if (m_verbosity > 5) {
    pout() << m_name + "::define" << endl;
  }

  m_physics               = a_physics;
  m_timeStepper           = a_timeStepper;
  m_amr                   = a_amr;
  m_computationalGeometry = a_computationalGeometry;
}

void
ItoPlasmaTagger::regrid()
{
  CH_TIME("ItoPlasmaTagger::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  if (m_num_tracers > 0) {
    m_tracer.resize(m_num_tracers);
    m_grad_tracer.resize(m_num_tracers);
    for (int i = 0; i < m_num_tracers; i++) {
      m_amr->allocate(m_tracer[i], m_realm, m_phase, 1);
      m_amr->allocate(m_grad_tracer[i], m_realm, m_phase, SpaceDim);
    }
  }
}

void
ItoPlasmaTagger::setPhase(const phase::which_phase a_phase)
{
  CH_TIME("ItoPlasmaTagger::setPhase");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

int
ItoPlasmaTagger::getNumberOfPlotVariables()
{
  CH_TIME("ItoPlasmaTagger::getNumberOfPlotVariables_cells");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  return m_num_tracers;
}

Vector<EBAMRCellData>&
ItoPlasmaTagger::getTracerFields()
{
  return m_tracer;
}

void
ItoPlasmaTagger::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp)
{
  CH_TIME("ItoPlasmaTagger::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  this->computeTracers();
  for (int i = 0; i < m_num_tracers; i++) {
    const EBAMRCellData& tracer = m_tracer[i];

    const Interval src_interv(0, 0);
    const Interval dst_interv(a_icomp, a_icomp);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      tracer[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);

      DataOps::setCoveredValue(*a_output[lvl], a_icomp, 0.0);
    }

    // Add component and name
    a_plotVariableNames.push_back("Tracer field-" + std::to_string(i));

    a_icomp++;
  }
}

bool
ItoPlasmaTagger::tagCells(EBAMRTags& a_tags)
{
  CH_TIME("ItoPlasmaTagger::tagCells");
  if (m_verbosity > 5) {
    pout() << m_name + "::tagCells" << endl;
  }

  bool got_new_tags = false;

  const RealVect origin       = m_amr->getProbLo();
  const Real     time         = m_timeStepper->getTime();
  const int      finest_level = m_amr->getFinestLevel();
  const int      max_depth    = m_amr->getMaxAmrDepth();
  const int finest_tag_level = (finest_level == max_depth) ? max_depth - 1 : finest_level; // Never tag on max_amr_depth

  if (m_num_tracers > 0) {

    // Compute tracer fields. This is an overriden pure function
    computeTracers();

    for (int lvl = 0; lvl <= finest_tag_level; lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const Box      box     = dbl.get(dit());
        const EBISBox& ebisbox = ebisl[dit()];

        const IntVectSet irreg_ivs = ebisbox.getIrregIVS(box);
        const IntVectSet prev_tags = IntVectSet((*a_tags[lvl])[dit()]);

        DenseIntVectSet coarsenTags(box, false); // Cells that will be coarsened
        DenseIntVectSet refine_tags(box, false); // Cells that will be refined

        Vector<EBCellFAB*> tracers;
        Vector<EBCellFAB*> gtracers;

        for (int i = 0; i < m_num_tracers; i++) {
          tracers.push_back(&((*m_tracer[i][lvl])[dit()]));
          gtracers.push_back(&((*m_grad_tracer[i][lvl])[dit()]));
        }

        DenseIntVectSet& tags = (*a_tags[lvl])[dit()];

        // Refinement and coarsening
        refineCellsBox(refine_tags, tracers, gtracers, lvl, box, ebisbox, time, dx, origin);
        coarsenCellsBox(coarsenTags, tracers, gtracers, lvl, box, ebisbox, time, dx, origin);

        // Check if we got any new tags, or we are just recycling old tags.
        // Basically we will check if (current_tags + refined_tags - coarsenTags) == current_tags
        DenseIntVectSet cpy1 = tags;
        tags -= coarsenTags;
        tags |= refine_tags;
        DenseIntVectSet cpy2 = tags;

        cpy2 -= cpy1; // = new tags minus old tags. If nonzero, we got some new tags.
        cpy1 -= tags; // = old_tags minus new tags. If nonzero, we got some new tags
        if (cpy1.numPts() != 0 || cpy2.numPts() != 0) {
          got_new_tags = true;
        }

        tags &= box;
      }
    }
  }

  // Some ranks may have gotten new tags while others have not. This little code snippet
  // sets got_new_tags = true for all ranks if any rank originally had got_new_tags = true
#ifdef CH_MPI
  int glo = 1;
  int loc = got_new_tags ? 1 : 0;

  MPI_Allreduce(&loc, &glo, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);

  got_new_tags = (glo == 1) ? true : false;
#endif

  return got_new_tags;
}

void
ItoPlasmaTagger::refineCellsBox(DenseIntVectSet&          a_refined_tags,
                                const Vector<EBCellFAB*>& a_tracers,
                                const Vector<EBCellFAB*>& a_grad_tracers,
                                const int                 a_lvl,
                                const Box                 a_box,
                                const EBISBox&            a_ebisbox,
                                const Real                a_time,
                                const Real                a_dx,
                                const RealVect            a_origin)
{

  Vector<BaseFab<Real>*> reg_tracers;
  Vector<BaseFab<Real>*> reg_gtracer;

  for (int i = 0; i < m_num_tracers; i++) {
    reg_tracers.push_back(&(a_tracers[i]->getSingleValuedFAB()));
    reg_gtracer.push_back(&(a_grad_tracers[i]->getSingleValuedFAB()));
  }

  // Regular box loop
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect  iv  = bit();
    const RealVect pos = a_origin + a_dx * RealVect(iv) + 0.5 * a_dx * RealVect::Unit;

    // If position is inside any of the tagging boxes, we can refine
    if (insideTagBox(pos) && a_ebisbox.isRegular(iv)) {

      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers);
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++) {
        tr[i] = (*reg_tracers[i])(iv, 0);
        gt[i] = RealVect(D_DECL((*reg_gtracer[i])(iv, 0), (*reg_gtracer[i])(iv, 1), (*reg_gtracer[i])(iv, 2)));
      }

      // Check if this cell should be refined
      const bool refine = refineCell(pos, a_time, a_dx, a_lvl, tr, gt);

      // If we refine, grow with buffer and increment to a_refined_tags
      if (refine) {
        a_refined_tags |= iv;
      }
    }
  }

  // Irregular box loop
  const IntVectSet& irreg   = a_ebisbox.getIrregIVS(a_box);
  const EBGraph&    ebgraph = a_ebisbox.getEBGraph();
  for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();
    const RealVect  pos = EBArith::getVofLocation(vof, a_dx * RealVect::Unit, a_origin);

    // If position is inside any of the tagging boxes, we can refine
    if (insideTagBox(pos)) {

      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers);
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++) {
        tr[i] = (*a_tracers[i])(vof, 0);
        gt[i] = RealVect(
          D_DECL((*a_grad_tracers[i])(vof, 0), (*a_grad_tracers[i])(vof, 1), (*a_grad_tracers[i])(vof, 2)));
      }

      // Check if this cell should be refined
      const bool refine = refineCell(pos, a_time, a_dx, a_lvl, tr, gt);

      if (refine) {
        a_refined_tags |= vof.gridIndex();
      }
    }
  }
}

void
ItoPlasmaTagger::coarsenCellsBox(DenseIntVectSet&          a_coarsened_tags,
                                 const Vector<EBCellFAB*>& a_tracers,
                                 const Vector<EBCellFAB*>& a_grad_tracers,
                                 const int                 a_lvl,
                                 const Box                 a_box,
                                 const EBISBox&            a_ebisbox,
                                 const Real                a_time,
                                 const Real                a_dx,
                                 const RealVect            a_origin)
{

  Vector<BaseFab<Real>*> reg_tracers;
  Vector<BaseFab<Real>*> reg_gtracer;

  for (int i = 0; i < m_num_tracers; i++) {
    reg_tracers.push_back(&(a_tracers[i]->getSingleValuedFAB()));
    reg_gtracer.push_back(&(a_grad_tracers[i]->getSingleValuedFAB()));
  }

  // Regular box loop
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect  iv  = bit();
    const RealVect pos = a_origin + a_dx * RealVect(iv);

    // If position is inside any of the tagging boxes, we can refine
    const bool inside = insideTagBox(pos);
    if (!inside) {
      a_coarsened_tags |= iv;
    }
    else if (inside && a_ebisbox.isRegular(iv)) {

      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers);
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++) {
        tr[i] = (*reg_tracers[i])(iv, 0);
        gt[i] = RealVect(D_DECL((*reg_gtracer[i])(iv, 0), (*reg_gtracer[i])(iv, 1), (*reg_gtracer[i])(iv, 2)));
      }

      // Check if this cell should be refined
      const bool coarsen = coarsenCell(pos, a_time, a_dx, a_lvl, tr, gt);

      // If we refine, grow with buffer and increment to a_refined_tags
      if (coarsen) {
        a_coarsened_tags |= iv;
      }
    }
  }

  // Irregular box loop
  const IntVectSet& irreg   = a_ebisbox.getIrregIVS(a_box);
  const EBGraph&    ebgraph = a_ebisbox.getEBGraph();
  for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();
    const RealVect  pos = EBArith::getVofLocation(vof, a_dx * RealVect::Unit, a_origin);

    // If position is inside any of the tagging boxes, we can refine
    const bool inside = insideTagBox(pos);
    if (!inside) {
      a_coarsened_tags |= vof.gridIndex();
    }
    else if (inside) {

      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers);
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++) {
        tr[i] = (*a_tracers[i])(vof, 0);
        gt[i] = RealVect(
          D_DECL((*a_grad_tracers[i])(vof, 0), (*a_grad_tracers[i])(vof, 1), (*a_grad_tracers[i])(vof, 2)));
      }

      // Check if this cell should be refined
      const bool coarsen = coarsenCell(pos, a_time, a_dx, a_lvl, tr, gt);

      if (coarsen) {
        a_coarsened_tags |= vof.gridIndex();
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
