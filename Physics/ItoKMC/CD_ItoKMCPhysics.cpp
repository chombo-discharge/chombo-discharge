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

thread_local bool                                            ItoKMCPhysics::m_hasKMCSolver;
thread_local KMCSolverType                                   ItoKMCPhysics::m_kmcSolver;
thread_local KMCState                                        ItoKMCPhysics::m_kmcState;
thread_local std::vector<std::shared_ptr<const KMCReaction>> ItoKMCPhysics::m_kmcReactionsThreadLocal;

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

bool
ItoKMCPhysics::needGradients() const noexcept
{
  CH_TIME("ItoKMCPhysics::needGradients");

  return false;
}

#include <CD_NamespaceFooter.H>
