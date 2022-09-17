/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaFieldTagger.cpp
  @brief  Implementation of CD_CdrPlasmaFieldTagger.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_CdrPlasmaFieldTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaFieldTagger::CdrPlasmaFieldTagger()
{
  CH_TIME("CdrPlasmaFieldTagger::CdrPlasmaFieldTagger");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaFieldTagger::CdrPlasmaFieldTagger" << endl;
  }

  m_name = "CdrPlasmaFieldTagger";
}

CdrPlasmaFieldTagger::~CdrPlasmaFieldTagger() {}

void
CdrPlasmaFieldTagger::allocateStorage() const
{
  CH_TIME("CdrPlasmaFieldTagger::allocateStorage()");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocateStorage()" << endl;
  }

  m_amr->allocate(m_scratch, m_realm, m_phase, 1);
  m_amr->allocate(m_electricField, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_gradElectricField, m_realm, m_phase, SpaceDim);
}

void
CdrPlasmaFieldTagger::deallocateStorage() const
{
  CH_TIME("CdrPlasmaFieldTagger::deallocateStorage()");
  if (m_verbosity > 5) {
    pout() << m_name + "::deallocateStorage()" << endl;
  }

  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_electricField);
  m_amr->deallocate(m_gradElectricField);
}

void
CdrPlasmaFieldTagger::computeElectricField(EBAMRCellData& a_electricField, EBAMRCellData& a_gradElectricField) const
{
  CH_TIME("CdrPlasmaFieldTagger::computeElectricField");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeElectricField" << endl;
  }

  // Time stepper computes the electric field for us.
  m_timeStepper->computeElectricField(a_electricField, m_phase);

  // Compute |E| onto scratch
  DataOps::vectorLength(m_scratch, a_electricField);

  // Now compute grad(|E|).
  m_amr->computeGradient(a_gradElectricField, m_scratch, m_realm, phase::gas);

  m_amr->conservativeAverage(a_gradElectricField, m_realm, m_phase);
  m_amr->interpGhost(a_gradElectricField, m_realm, m_phase);

  // Interpolate everything to centroids since that is the only place where things make sense.
  m_amr->interpToCentroids(a_electricField, m_realm, m_phase);
  m_amr->interpToCentroids(a_gradElectricField, m_realm, m_phase);
}

void
CdrPlasmaFieldTagger::computeTracers() const
{
  CH_TIME("CdrPlasmaFieldTagger::computeTracers()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeTracers()" << endl;
  }

  constexpr int comp = 0;

  // Allocate necessary storage.
  this->allocateStorage();

  // Get the lower-left domain coordinate and current time.
  const RealVect probLo = m_amr->getProbLo();
  const Real     time   = m_timeStepper->getTime();

  // Compute electric field on volumetric centroids
  this->computeElectricField(m_electricField, m_gradElectricField);

  // Get maximum and minimum of everything
  Real maxElectricField     = -std::numeric_limits<Real>::max();
  Real minElectricField     = std::numeric_limits<Real>::max();
  Real maxGradElectricField = -std::numeric_limits<Real>::max();
  Real minGradElectricField = std::numeric_limits<Real>::max();

  // Get the maximum and minium value of the electric field and its gradient. This is the
  // norm so we are getting the max of |E| and the max of |grad(|E|)|.
  DataOps::getMaxMinNorm(maxElectricField, minElectricField, m_electricField);
  DataOps::getMaxMinNorm(maxGradElectricField, minGradElectricField, m_gradElectricField);

  // Go through each AMR level and compute the tracer fields.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    // Iterate through the grid patches.
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box&     box     = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      // Expose the AMR data per patch, both full data and the single-valued data.
      const EBCellFAB& electricField     = (*m_electricField[lvl])[dit()];
      const EBCellFAB& gradElectricField = (*m_gradElectricField[lvl])[dit()];

      const FArrayBox& electricFieldReg     = electricField.getFArrayBox();
      const FArrayBox& gradElectricFieldReg = gradElectricField.getFArrayBox();

      // Get references to tracer fields, both the full data and the single-valued data. These are needed
      // in the regular and irregular kernels.
      Vector<EBCellFAB*> tr;
      Vector<FArrayBox*> trReg;

      for (int i = 0; i < m_numTracers; i++) {
        tr.push_back(&((*m_tracers[i][lvl])[dit()]));
        trReg.push_back(&(tr[i]->getFArrayBox()));
      }

      // Regular kernel
      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

        // Reconstruct the electric field and gradient of the electric field.
        const RealVect E = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const RealVect gradE = RealVect(
          D_DECL(gradElectricFieldReg(iv, 0), gradElectricFieldReg(iv, 1), gradElectricFieldReg(iv, 2)));

        // Call the per-point tracer function. This is a pure function.
        const Vector<Real> tracers = this->tracer(pos,
                                                  time,
                                                  dx,
                                                  E,
                                                  minElectricField,
                                                  maxElectricField,
                                                  gradE,
                                                  minGradElectricField,
                                                  maxGradElectricField);

        // Put the tracer field where it belongs.
        for (int i = 0; i < m_numTracers; i++) {
          (*trReg[i])(iv, comp) = tracers[i];
        }
      };

      // Irregular box loop
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Center, vof, ebisbox, dx);

        // Reconstruct the electric field and gradient of the electric field.
        const RealVect E     = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const RealVect gradE = RealVect(
          D_DECL(gradElectricField(vof, 0), gradElectricField(vof, 1), gradElectricField(vof, 2)));

        // Call the per-point tracer function. Again, it's a pure function.
        const Vector<Real> tracers = this->tracer(pos,
                                                  time,
                                                  dx,
                                                  E,
                                                  minElectricField,
                                                  maxElectricField,
                                                  gradE,
                                                  minGradElectricField,
                                                  maxGradElectricField);

        // Put the tracer field where it belongs.
        for (int i = 0; i < m_numTracers; i++) {
          (*tr[i])(vof, comp);
        }
      };

      // Irregular kernel region
      VoFIterator vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      // Execute the kernels
      BoxLoops::loop(box, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  // Coarsen the data.
  for (int i = 0; i < m_numTracers; i++) {
    m_amr->conservativeAverage(m_tracers[i], m_realm, m_phase);
    m_amr->interpGhost(m_tracers[i], m_realm, m_phase);
  }

  // Compute gradients.
  for (int i = 0; i < m_numTracers; i++) {
    m_amr->computeGradient(m_gradTracers[i], m_tracers[i], m_realm, phase::gas);
    m_amr->conservativeAverage(m_gradTracers[i], m_realm, m_phase);
  }

  // No need to keep this transient storage lying around so we delete it.
  this->deallocateStorage();
}

#include <CD_NamespaceFooter.H>
