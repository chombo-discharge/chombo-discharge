/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaFieldTaggerImplem.H
  @brief  Implementation of CD_ItoPlasmaFieldTagger.H
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaFieldTaggerImplem_H
#define CD_ItoPlasmaFieldTaggerImplem_H

// Our includes
#include <CD_ItoPlasmaFieldTagger.H>
#include <CD_DataOps.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

template <typename S>
ItoPlasmaFieldTagger<S>::ItoPlasmaFieldTagger()
{
  CH_TIME("ItoPlasmaFieldTagger::ItoPlasmaFieldTagger");
  if (this->m_verbosity > 5) {
    pout() << "ItoPlasmaFieldTagger::ItoPlasmaFieldTagger" << endl;
  }

  this->m_name = "ItoPlasmaFieldTagger";
}

template <typename S>
ItoPlasmaFieldTagger<S>::~ItoPlasmaFieldTagger() {}

template <typename S>
void
ItoPlasmaFieldTagger<S>::allocateStorage()
{
  CH_TIME("ItoPlasmaFieldTagger::allocateStorage");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::allocateStorage" << endl;
  }

  this->m_amr->allocate(m_scratch, this->m_realm, this->m_phase, 1);
  this->m_amr->allocate(m_E, this->m_realm, this->m_phase, SpaceDim);
  this->m_amr->allocate(m_grad_E, this->m_realm, this->m_phase, SpaceDim);
}

template <typename S>
void
ItoPlasmaFieldTagger<S>::deallocateStorage()
{
  CH_TIME("ItoPlasmaFieldTagger::deallocateStorage");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::deallocateStorage" << endl;
  }

  this->m_amr->deallocate(m_scratch);
  this->m_amr->deallocate(m_E);
  this->m_amr->deallocate(m_grad_E);
}

template <typename S>
void
ItoPlasmaFieldTagger<S>::computeElectricField(EBAMRCellData& a_E, EBAMRCellData& a_grad_E)
{
  CH_TIME("ItoPlasmaFieldTagger::computeElectricField");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeElectricField" << endl;
  }

  this->m_timeStepper->computeElectricField(a_E, this->m_phase);
  DataOps::vectorLength(m_scratch, a_E);
  this->m_amr->computeGradient(a_grad_E, m_scratch, this->m_realm, phase::gas);

  this->m_amr->conservativeAverage(a_grad_E, this->m_realm, this->m_phase);
  this->m_amr->interpGhost(a_grad_E, this->m_realm, this->m_phase);

  // Interpolate to centroids
  this->m_amr->interpToCentroids(a_E, this->m_realm, this->m_phase);
  this->m_amr->interpToCentroids(a_grad_E, this->m_realm, this->m_phase);
}

template <typename S>
void
ItoPlasmaFieldTagger<S>::computeTracers()
{
  CH_TIME("ItoPlasmaFieldTagger::computeTracers");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeTracers" << endl;
  }

  this->allocateStorage();

  const RealVect origin = this->m_amr->getProbLo();
  const Real     time   = this->m_timeStepper->getTime();

  // Compute electric field on volumetric centroids
  this->computeElectricField(m_E, m_grad_E);

  // Get maximum and minimum of everything
  Real E_max, E_min;
  Real grad_E_max, grad_E_min;

  DataOps::getMaxMinNorm(E_max, E_min, m_E);
  DataOps::getMaxMinNorm(grad_E_max, grad_E_min, m_grad_E);

  for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = this->m_amr->getGrids(this->m_realm)[lvl];
    const EBISLayout&        ebisl = this->m_amr->getEBISLayout(this->m_realm, this->m_phase)[lvl];
    const Real               dx    = this->m_amr->getDx()[lvl];

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
      for (int i = 0; i < this->m_num_tracers; i++) {
        tr.push_back(&((*this->m_tracer[i][lvl])[dit()]));
        tr_fab.push_back(&(tr[i]->getSingleValuedFAB()));
      }

      // Regular box loop
      for (BoxIterator bit(box); bit.ok(); ++bit) {
        const IntVect  iv  = bit();
        const RealVect pos = origin + RealVect(iv) * dx;

        const RealVect E      = RealVect(D_DECL(E_reg(iv, 0), E_reg(iv, 1), E_reg(iv, 2)));
        const RealVect grad_E = RealVect(D_DECL(gE_reg(iv, 0), gE_reg(iv, 1), gE_reg(iv, 2)));

        Vector<Real> tracers = this->tracer(pos, time, dx, E, E_min, E_max, grad_E, grad_E_min, grad_E_max);

        for (int i = 0; i < this->m_num_tracers; i++) {
          (*tr_fab[i])(iv, 0) = tracers[i];
        }
      }

      // Irregular box loop
      const IntVectSet& irreg = ebisbox.getIrregIVS(box);
      for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit) {
        const VolIndex& vof = vofit();
        const RealVect  pos = Location::position(Location::Cell::Center, vof, ebisbox, dx);

        // Electric field and grad(|E|)
        const RealVect E      = RealVect(D_DECL(E_fab(vof, 0), E_fab(vof, 1), E_fab(vof, 2)));
        const RealVect grad_E = RealVect(D_DECL(gE_fab(vof, 0), gE_fab(vof, 1), gE_fab(vof, 2)));

        Vector<Real> tracers = this->tracer(pos, time, dx, E, E_min, E_max, grad_E, grad_E_min, grad_E_max);

        for (int i = 0; i < this->m_num_tracers; i++) {
          (*tr[i])(vof, 0);
        }
      }
    }
  }

  for (int i = 0; i < this->m_num_tracers; i++) {
    this->m_amr->conservativeAverage(this->m_tracer[i], this->m_realm, this->m_phase);
    this->m_amr->interpGhost(this->m_tracer[i], this->m_realm, this->m_phase);
  }

  // Compute gradient of tracers
  for (int i = 0; i < this->m_num_tracers; i++) {
    this->m_amr->computeGradient(this->m_grad_tracer[i], this->m_tracer[i], this->m_realm, this->m_phase);
    this->m_amr->conservativeAverage(this->m_grad_tracer[i], this->m_realm, this->m_phase);
  }

  this->deallocateStorage(); // No reason to keep the extra storage lying around...
}

#include <CD_NamespaceFooter.H>

#endif
