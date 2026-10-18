/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzElectrostaticDomainBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzElectrostaticDomainBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzElectrostaticDomainBC.H>
#include <CD_PetscHelmholtzDirichletDomainBC.H>
#include <CD_PetscHelmholtzNeumannDomainBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzElectrostaticDomainBC::PetscHelmholtzElectrostaticDomainBC(
  const ElectrostaticDomainBc& a_electrostaticBCs)
  : m_electrostaticBCs(a_electrostaticBCs)
{
  CH_TIME("PetscHelmholtzElectrostaticDomainBC::PetscHelmholtzElectrostaticDomainBC");

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const ElectrostaticDomainBc::DomainSide domainSide = std::make_pair(dir, sit());
      const ElectrostaticDomainBc::BcType&    bcType     = m_electrostaticBCs.getBc(domainSide).first;

      // The BC objects take a function of position only. As in EBHelmholtzElectrostaticDomainBC, time enters because
      // FieldSolver passes in BC functions that capture its time by reference.
      auto func = [domainSide, &BC = this->m_electrostaticBCs](const RealVect& a_position) -> Real {
        constexpr Real dummyTime = 0.0;

        return BC.getBc(domainSide).second(a_position, dummyTime);
      };

      switch (bcType) {
      case ElectrostaticDomainBc::BcType::Dirichlet: {
        m_bcObjects.emplace(domainSide, std::make_shared<PetscHelmholtzDirichletDomainBC>(func));

        break;
      }
      case ElectrostaticDomainBc::BcType::Neumann: {
        m_bcObjects.emplace(domainSide, std::make_shared<PetscHelmholtzNeumannDomainBC>(func));

        break;
      }
      default: {
        MayDay::Abort("PetscHelmholtzElectrostaticDomainBC -- unsupported boundary condition");

        break;
      }
      }
    }
  }
}

PetscHelmholtzElectrostaticDomainBC::~PetscHelmholtzElectrostaticDomainBC()
{
  CH_TIME("PetscHelmholtzElectrostaticDomainBC::~PetscHelmholtzElectrostaticDomainBC");
}

void
PetscHelmholtzElectrostaticDomainBC::define(const Location::Cell       a_dataLocation,
                                            const Vector<EBLevelGrid>& a_eblg,
                                            const RealVect&            a_probLo,
                                            const Vector<Real>&        a_dx)
{
  CH_TIME("PetscHelmholtzElectrostaticDomainBC::define");

  PetscHelmholtzDomainBC::define(a_dataLocation, a_eblg, a_probLo, a_dx);

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      m_bcObjects.at(std::make_pair(dir, sit()))->define(a_dataLocation, a_eblg, a_probLo, a_dx);
    }
  }
}

std::pair<VoFStencil, VoFStencil>
PetscHelmholtzElectrostaticDomainBC::getFaceFluxStencil(const VolIndex&      a_vof,
                                                        const int            a_dir,
                                                        const Side::LoHiSide a_side,
                                                        const int            a_level,
                                                        const DataIndex&     a_din) const
{
  return m_bcObjects.at(std::make_pair(a_dir, a_side))->getFaceFluxStencil(a_vof, a_dir, a_side, a_level, a_din);
}

Real
PetscHelmholtzElectrostaticDomainBC::getFaceFluxConstant(const VolIndex&      a_vof,
                                                         const int            a_dir,
                                                         const Side::LoHiSide a_side,
                                                         const int            a_level,
                                                         const DataIndex&     a_din) const
{
  return m_bcObjects.at(std::make_pair(a_dir, a_side))->getFaceFluxConstant(a_vof, a_dir, a_side, a_level, a_din);
}

#include <CD_NamespaceFooter.H>

#endif
