/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCPhysics.cpp
  @brief  Implementation of CD_ItoKMCPhysics.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ItoKMCPhysics.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

Vector<std::string>
ItoKMCPhysics::getPlotVariableNames() const noexcept
{
  CH_TIME("ItoKMCPhysics::getPlotVariableNames");

  return Vector<std::string>(0);
}

Vector<Real>
ItoKMCPhysics::getPlotVariables(const RealVect          a_E,
                                const RealVect          a_pos,
                                const Vector<Real>&     a_phi,
                                const Vector<RealVect>& a_gradPhi,
                                const Real              a_dx,
                                const Real              a_kappa) const noexcept
{
  CH_TIME("ItoKMCPhysics::getPlotVariables");

  return Vector<Real>(0);
}

int
ItoKMCPhysics::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("ItoKMCPhysics::getNumberOfPlotVariables");

  return 0;
}

#include <CD_NamespaceFooter.H>
