/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaFieldTagger.cpp
  @brief  Implementation of CD_ItoPlasmaFieldTagger.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_ItoPlasmaFieldTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaFieldTagger::ItoPlasmaFieldTagger()
{
  CH_TIME("ItoPlasmaFieldTagger::ItoPlasmaFieldTagger");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaFieldTagger::ItoPlasmaFieldTagger" << endl;
  }

  m_name = "ItoPlasmaFieldTagger";
}

ItoPlasmaFieldTagger::~ItoPlasmaFieldTagger() {}

void
ItoPlasmaFieldTagger::allocateStorage()
{
  CH_TIME("ItoPlasmaFieldTagger::allocateStorage");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocateStorage" << endl;
  }

  m_amr->allocate(m_scratch, m_realm, m_phase, 1);
  m_amr->allocate(m_E, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_grad_E, m_realm, m_phase, SpaceDim);
}

void
ItoPlasmaFieldTagger::deallocateStorage()
{
  CH_TIME("ItoPlasmaFieldTagger::deallocateStorage");
  if (m_verbosity > 5) {
    pout() << m_name + "::deallocateStorage" << endl;
  }

  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_E);
  m_amr->deallocate(m_grad_E);
}

void
ItoPlasmaFieldTagger::computeElectricField(EBAMRCellData& a_E, EBAMRCellData& a_grad_E)
{
  CH_TIME("ItoPlasmaFieldTagger::computeElectricField");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeElectricField" << endl;
  }

  m_timeStepper->computeElectricField(a_E, m_phase);
  DataOps::vectorLength(m_scratch, a_E);
  m_amr->computeGradient(a_grad_E, m_scratch, m_realm, phase::gas);

  m_amr->conservativeAverage(a_grad_E, m_realm, m_phase);
  m_amr->interpGhost(a_grad_E, m_realm, m_phase);

  // Interpolate to centroids
  m_amr->interpToCentroids(a_E, m_realm, m_phase);
  m_amr->interpToCentroids(a_grad_E, m_realm, m_phase);
}

void
ItoPlasmaFieldTagger::computeTracers()
{
  CH_TIME("ItoPlasmaFieldTagger::computeTracers");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeTracers" << endl;
  }

  this->allocateStorage();

  const RealVect origin = m_amr->getProbLo();
  const Real     time   = m_timeStepper->getTime();

  // Compute electric field on volumetric centroids
  this->computeElectricField(m_E, m_grad_E);

  // Get maximum and minimum of everything
  Real E_max, E_min;
  Real grad_E_max, grad_E_min;

  DataOps::getMaxMinNorm(E_max, E_min, m_E);
  DataOps::getMaxMinNorm(grad_E_max, grad_E_min, m_grad_E);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const Box&       box     = dbl.get(dit());
      const EBISBox&   ebisbox = ebisl[dit()];
      const EBGraph&   ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      const EBCellFAB& E_fab  = (*m_E[lvl])[dit()];
      const EBCellFAB& gE_fab = (*m_grad_E[lvl])[dit()];

      const BaseFab<Real>& E_reg  = E_fab.getSingleValuedFAB();
      const BaseFab<Real>& gE_reg = gE_fab.getSingleValuedFAB();

      // Avoid the extra point lookups by getting these before the point loops
      Vector<EBCellFAB*>     tr;
      Vector<BaseFab<Real>*> tr_fab;
      for (int i = 0; i < m_num_tracers; i++) {
        tr.push_back(&((*m_tracer[i][lvl])[dit()]));
        tr_fab.push_back(&(tr[i]->getSingleValuedFAB()));
      }

      // Regular box loop
      for (BoxIterator bit(box); bit.ok(); ++bit) {
        const IntVect  iv  = bit();
        const RealVect pos = origin + RealVect(iv) * dx;

        const RealVect E      = RealVect(D_DECL(E_reg(iv, 0), E_reg(iv, 1), E_reg(iv, 2)));
        const RealVect grad_E = RealVect(D_DECL(gE_reg(iv, 0), gE_reg(iv, 1), gE_reg(iv, 2)));

        Vector<Real> tracers = this->tracer(pos, time, dx, E, E_min, E_max, grad_E, grad_E_min, grad_E_max);

        for (int i = 0; i < m_num_tracers; i++) {
          (*tr_fab[i])(iv, 0) = tracers[i];
        }
      }

      // Irregular box loop
      const IntVectSet& irreg = ebisbox.getIrregIVS(box);
      for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit) {
        const VolIndex& vof = vofit();
        const RealVect  pos = EBArith::getVofLocation(vof, dx * RealVect::Unit, origin);

        // Electric field and grad(|E|)
        const RealVect E      = RealVect(D_DECL(E_fab(vof, 0), E_fab(vof, 1), E_fab(vof, 2)));
        const RealVect grad_E = RealVect(D_DECL(gE_fab(vof, 0), gE_fab(vof, 1), gE_fab(vof, 2)));

        Vector<Real> tracers = this->tracer(pos, time, dx, E, E_min, E_max, grad_E, grad_E_min, grad_E_max);

        for (int i = 0; i < m_num_tracers; i++) {
          (*tr[i])(vof, 0);
        }
      }
    }
  }

  for (int i = 0; i < m_num_tracers; i++) {
    m_amr->conservativeAverage(m_tracer[i], m_realm, m_phase);
    m_amr->interpGhost(m_tracer[i], m_realm, m_phase);
  }

  // Compute gradient of tracers
  for (int i = 0; i < m_num_tracers; i++) {
    m_amr->computeGradient(m_grad_tracer[i], m_tracer[i], m_realm, m_phase);
    m_amr->conservativeAverage(m_grad_tracer[i], m_realm, m_phase);
  }

  this->deallocateStorage(); // No reason to keep the extra storage lying around...
}

#include <CD_NamespaceFooter.H>
