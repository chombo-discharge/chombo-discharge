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
PetscHelmholtzElectrostaticDomainBC::define(const RefCountedPtr<PetscGrid>& a_petscGrid,
                                            const RealVect&                 a_probLo,
                                            const Vector<Real>&             a_dx)
{
  CH_TIME("PetscHelmholtzElectrostaticDomainBC::define");

  PetscHelmholtzDomainBC::define(a_petscGrid, a_probLo, a_dx);

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      m_bcObjects.at(std::make_pair(dir, sit()))->define(a_petscGrid, a_probLo, a_dx);
    }
  }
}

void
PetscHelmholtzElectrostaticDomainBC::getFaceFluxStencil(Real&                a_boundaryWeight,
                                                        VoFStencil&          a_stencilCoar,
                                                        VoFStencil&          a_stencil,
                                                        VoFStencil&          a_stencilFine,
                                                        const int            a_phase,
                                                        const int            a_level,
                                                        const DataIndex&     a_din,
                                                        const VolIndex&      a_vof,
                                                        const int            a_dir,
                                                        const Side::LoHiSide a_side) const
{
  m_bcObjects.at(std::make_pair(a_dir, a_side))
    ->getFaceFluxStencil(a_boundaryWeight,
                         a_stencilCoar,
                         a_stencil,
                         a_stencilFine,
                         a_phase,
                         a_level,
                         a_din,
                         a_vof,
                         a_dir,
                         a_side);
}

Real
PetscHelmholtzElectrostaticDomainBC::getBoundaryData(const int            a_phase,
                                                     const int            a_level,
                                                     const DataIndex&     a_din,
                                                     const VolIndex&      a_vof,
                                                     const int            a_dir,
                                                     const Side::LoHiSide a_side) const
{
  return m_bcObjects.at(std::make_pair(a_dir, a_side))->getBoundaryData(a_phase, a_level, a_din, a_vof, a_dir, a_side);
}

#include <CD_NamespaceFooter.H>

#endif
