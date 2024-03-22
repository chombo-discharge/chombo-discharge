/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBHelmholtzLarsenEBBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzLarsenEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzLarsenEBBCFactory.H>
#include <CD_EBHelmholtzRobinEBBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzLarsenEBBCFactory::EBHelmholtzLarsenEBBCFactory(const int                       a_order,
                                                           const int                       a_weight,
                                                           const RefCountedPtr<RtSpecies>& a_species,
                                                           const Real                      a_r1,
                                                           const Real                      a_r2,
                                                           const Real                      a_source)
{
  m_order   = a_order;
  m_weight  = a_weight;
  m_species = a_species;
  m_r1      = a_r1;
  m_r2      = a_r2;
  m_source  = a_source;

  this->setRobinCoefficients();
}

EBHelmholtzLarsenEBBCFactory::~EBHelmholtzLarsenEBBCFactory()
{}

void
EBHelmholtzLarsenEBBCFactory::setRobinCoefficients()
{
  // TLDR: We are creating functions for boundary conditions of the Robin type A*phi + B*dphi/dn = C, but using
  //       the "Larsen coefficients".

  // Using our formulation of the Helmholtz operator, The A-coefficient is 1.5*kappa*kappa*(1+3*r2)/(1-2*r1)
  m_functionA = [&species = this->m_species, r1 = this->m_r1, r2 = this->m_r2](const RealVect& a_position) -> Real {
    Real val = species->getAbsorptionCoefficient(a_position);

    val *= val;
    val *= 3. / 2.;
    val *= (1 + 3 * r1) / (1 - 2 * r2);

    return val;
  };

  m_functionB = [&species = this->m_species](const RealVect& a_position) {
    return -species->getAbsorptionCoefficient(a_position);
  };

  // This is the right-hand side of the Robin BC, i.e. the source function. Time is a dummy parameter, and the user should
  // have captured some external time (e.g., RtSolver::m_time) by reference in the function that was passed into the full constructor.
  m_functionC = [source = this->m_source](const RealVect& a_position) {
    return source;
  };

  m_useConstant = false;
  m_useFunction = true;
}

#include <CD_NamespaceFooter.H>
