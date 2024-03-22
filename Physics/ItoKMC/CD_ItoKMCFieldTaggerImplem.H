/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCFieldTaggerImplem.H
  @brief  Implementation of CD_ItoKMCFieldTagger.H
  @author Robert Marskar
*/

#ifndef CD_ItoKMCFieldTaggerImplem_H
#define CD_ItoKMCFieldTaggerImplem_H

// Our includes
#include <CD_ItoKMCFieldTagger.H>
#include <CD_DataOps.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

template <typename S>
ItoKMCFieldTagger<S>::ItoKMCFieldTagger() noexcept
{
  CH_TIME("ItoKMCFieldTagger::ItoKMCFieldTagger");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCFieldTagger::ItoKMCFieldTagger" << endl;
  }

  this->m_name = "ItoKMCFieldTagger";
}

template <typename S>
ItoKMCFieldTagger<S>::~ItoKMCFieldTagger() noexcept
{}

template <typename S>
void
ItoKMCFieldTagger<S>::allocateStorage() const noexcept
{
  CH_TIME("ItoKMCFieldTagger::allocateStorage");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::allocateStorage" << endl;
  }

  CH_assert(this->m_isDefined);

  (this->m_amr)->allocate(m_scratch, this->m_realm, this->m_phase, 1);
  (this->m_amr)->allocate(m_E, this->m_realm, this->m_phase, SpaceDim);
  (this->m_amr)->allocate(m_gradE, this->m_realm, this->m_phase, SpaceDim);
}

template <typename S>
void
ItoKMCFieldTagger<S>::deallocateStorage() const noexcept
{
  CH_TIME("ItoKMCFieldTagger::deallocateStorage");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::deallocateStorage" << endl;
  }

  this->m_amr->deallocate(m_scratch);
  this->m_amr->deallocate(m_E);
  this->m_amr->deallocate(m_gradE);
}

template <typename S>
void
ItoKMCFieldTagger<S>::computeElectricField(EBAMRCellData& a_E, EBAMRCellData& a_gradE) const noexcept
{
  CH_TIME("ItoKMCFieldTagger::computeElectricField");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeElectricField" << endl;
  }

  CH_assert(this->m_isDefined);
  CH_assert(a_E[0]->nComp() == SpaceDim);
  CH_assert(a_gradE[0]->nComp() == SpaceDim);

  this->m_timeStepper->computeElectricField(a_E, this->m_phase);
  DataOps::vectorLength(m_scratch, a_E);
  this->m_amr->computeGradient(a_gradE, m_scratch, this->m_realm, this->m_phase);

  this->m_amr->conservativeAverage(a_gradE, this->m_realm, this->m_phase);
  this->m_amr->interpGhost(a_gradE, this->m_realm, this->m_phase);

  // Interpolate to centroids
  this->m_amr->interpToCentroids(a_E, this->m_realm, this->m_phase);
  this->m_amr->interpToCentroids(a_gradE, this->m_realm, this->m_phase);
}

template <typename S>
void
ItoKMCFieldTagger<S>::computeTagFields() const noexcept
{
  CH_TIME("ItoKMCFieldTagger::computeTagFields");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeTagFields" << endl;
  }

  CH_assert(this->m_isDefined);

  this->allocateStorage();
  this->computeElectricField(m_E, m_gradE);

  const RealVect probLo = this->m_amr->getProbLo();
  const Real     time   = this->m_timeStepper->getTime();

  Real maxE, minE;
  Real maxGradE, minGradE;

  DataOps::getMaxMinNorm(maxE, minE, m_E);
  DataOps::getMaxMinNorm(maxGradE, minGradE, m_gradE);

  for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = this->m_amr->getGrids(this->m_realm)[lvl];
    const EBISLayout&        ebisl = this->m_amr->getEBISLayout(this->m_realm, this->m_phase)[lvl];
    const Real               dx    = this->m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box&     box     = dbl[dit()];
      const EBISBox& ebisbox = ebisl[dit()];

      const EBCellFAB& electricField     = (*m_E[lvl])[dit()];
      const EBCellFAB& gradElectricField = (*m_gradE[lvl])[dit()];

      const FArrayBox& electricFieldReg     = electricField.getFArrayBox();
      const FArrayBox& gradElectricFieldReg = gradElectricField.getFArrayBox();

      Vector<EBCellFAB*> tagFields;
      Vector<FArrayBox*> tagFieldsReg;
      for (int i = 0; i < this->m_numTagFields; i++) {
        tagFields.push_back(&((*this->m_tagFields[i][lvl])[dit()]));
        tagFieldsReg.push_back(&(tagFields[i]->getFArrayBox()));
      }

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv)) {
          const RealVect pos = probLo + RealVect(iv) * dx;

          const RealVect E = RealVect(
            D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
          const RealVect gradE = RealVect(
            D_DECL(gradElectricFieldReg(iv, 0), gradElectricFieldReg(iv, 1), gradElectricFieldReg(iv, 2)));

          const Vector<Real>
            compTagFields = this->computeTagFields(pos, time, dx, E, minE, maxE, gradE, minGradE, maxGradE);

          for (int i = 0; i < this->m_numTagFields; i++) {
            (*tagFieldsReg[i])(iv, 0) = compTagFields[i];
          }
        }
      };

      // Irregular box loop
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Center, vof, ebisbox, dx);

        const RealVect E     = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const RealVect gradE = RealVect(
          D_DECL(gradElectricField(vof, 0), gradElectricField(vof, 1), gradElectricField(vof, 2)));

        const Vector<Real>
          compTagFields = this->computeTagFields(pos, time, dx, E, minE, maxE, gradE, minGradE, maxGradE);

        for (int i = 0; i < this->m_numTagFields; i++) {
          (*tagFields[i])(vof, 0) = compTagFields[i];
        }
      };

      // Run the kernels
      VoFIterator& vofit = (*this->m_amr->getVofIterator(this->m_realm, this->m_phase)[lvl])[dit()];
      BoxLoops::loop(box, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }

    for (int i = 0; i < this->m_numTagFields; i++) {
      this->m_amr->conservativeAverage(this->m_tagFields[i], this->m_realm, this->m_phase);
      this->m_amr->interpGhost(this->m_tagFields[i], this->m_realm, this->m_phase);
    }

    // Compute gradient of tracers
    for (int i = 0; i < this->m_numTagFields; i++) {
      this->m_amr->computeGradient(this->m_gradTagFields[i], this->m_tagFields[i], this->m_realm, this->m_phase);
      this->m_amr->conservativeAverage(this->m_gradTagFields[i], this->m_realm, this->m_phase);
    }
  }

  this->deallocateStorage();
}

#include <CD_NamespaceFooter.H>

#endif
