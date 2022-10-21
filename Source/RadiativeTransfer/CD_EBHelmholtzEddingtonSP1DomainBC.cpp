/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEddingtonSP1DomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzEddingtonSP1DomainBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzEddingtonSP1DomainBC.H>
#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_EBHelmholtzNeumannDomainBC.H>
#include <CD_EBHelmholtzLarsenDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzEddingtonSP1DomainBC::EBHelmholtzEddingtonSP1DomainBC(const EddingtonSP1DomainBc&     a_eddingtonBCs,
                                                                 const RefCountedPtr<RtSpecies>& a_species,
                                                                 const Real                      a_r1,
                                                                 const Real                      a_r2)
{
  CH_TIME("EBHelmholtzEddingtonSP1DomainBC::EBHelmholtzEddingtonSP1DomainBC");

  m_eddingtonBCs = a_eddingtonBCs;
  m_species      = a_species;
  m_r1           = a_r1;
  m_r2           = a_r2;

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const EddingtonSP1DomainBc::DomainSide domainSide = std::make_pair(dir, sit());
      const EddingtonSP1DomainBc::BcType&    bcType     = m_eddingtonBCs.getBc(domainSide).first;

      // Make a lambda which allows us to pass in the function using EBHelmholtzDomainBC API, which takes a std::function<Real(const RealVect a_position)>
      // type of function.
      //
      // This is the function type that the EBHelmholtzOp API requires, and it is a design choice mandated by our choice to make that operator
      // time-independent. Although this might seem weird, the time dependence is nonetheless passed in because a_eddingtonSP1BCs are passed in
      // from EddingtonSP1, and in that solver we capture RtSolver::m_time by reference.
      //
      // This might seem clunky, but I can't see any other way of doing it properly without changing EBHelmholtzOp to a time-dependent operator (which
      // I really don't want to do).
      auto func = [domainSide, &BC = this->m_eddingtonBCs](const RealVect& a_position) -> Real {
        const Real dummyDt = 0.0;

        // This is a function Real(const RealVect, const Real)
        const EddingtonSP1DomainBc::BcFunction& bcFunction = BC.getBc(domainSide).second;

        return bcFunction(a_position, dummyDt);
      };

      // Create either a Dirichlet or Neumann boundary condition object for the current domain edge (face in 3D).
      switch (bcType) {
      case EddingtonSP1DomainBc::BcType::Dirichlet: {
        m_bcObjects.emplace(domainSide, std::make_shared<EBHelmholtzDirichletDomainBC>(func));

        break;
      }
      case EddingtonSP1DomainBc::BcType::Neumann: {
        m_bcObjects.emplace(domainSide, std::make_shared<EBHelmholtzNeumannDomainBC>(func));

        break;
      }
      case EddingtonSP1DomainBc::BcType::Larsen: {
        m_bcObjects.emplace(domainSide, std::make_shared<EBHelmholtzLarsenDomainBC>(m_species, m_r1, m_r2, func));

        break;
      }
      default: {
        MayDay::Error(
          "EBHelmholtzEddingtonSP1DomainBC::EBHelmholtzEddingtonSP1DomainBC - unsupported boundary condition passed into constructor!");

        break;
      }
      }
    }
  }
}

EBHelmholtzEddingtonSP1DomainBC::~EBHelmholtzEddingtonSP1DomainBC()
{
  CH_TIME("EBHelmholtzEddingtonSP1DomainBC::~EBHelmholtzEddingtonSP1DomainBC");
}

void
EBHelmholtzEddingtonSP1DomainBC::define(const Location::Cell a_dataLocation,
                                        const EBLevelGrid&   a_eblg,
                                        const RealVect&      a_probLo,
                                        const Real           a_dx)
{
  CH_TIME("EBHelmholtzEddingtonSP1DomainBC::define");

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      auto& bcPtr = m_bcObjects.at(std::make_pair(dir, sit()));

      bcPtr->define(a_dataLocation, a_eblg, a_probLo, a_dx);
    }
  }
}

void
EBHelmholtzEddingtonSP1DomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                             const BaseFab<Real>&  a_phi,
                                             const BaseFab<Real>&  a_Bcoef,
                                             const int&            a_dir,
                                             const Side::LoHiSide& a_side,
                                             const DataIndex&      a_dit,
                                             const bool            a_useHomogeneous) const
{
  CH_TIME("EBHelmholtzEddingtonSP1DomainBC::getFaceFlux(regular)");

  const auto& bcPtr = m_bcObjects.at(std::make_pair(a_dir, a_side));

  bcPtr->getFaceFlux(a_faceFlux, a_phi, a_Bcoef, a_dir, a_side, a_dit, a_useHomogeneous);
}

Real
EBHelmholtzEddingtonSP1DomainBC::getFaceFlux(const VolIndex&       a_vof,
                                             const EBCellFAB&      a_phi,
                                             const EBFaceFAB&      a_Bcoef,
                                             const int&            a_dir,
                                             const Side::LoHiSide& a_side,
                                             const DataIndex&      a_dit,
                                             const bool            a_useHomogeneous) const
{
  CH_TIME("EBHelmholtzEddingtonSP1DomainBC::getFaceFlux(irreg)");

  const auto& bcPtr = m_bcObjects.at(std::make_pair(a_dir, a_side));

  return bcPtr->getFaceFlux(a_vof, a_phi, a_Bcoef, a_dir, a_side, a_dit, a_useHomogeneous);
}

#include <CD_NamespaceFooter.H>
