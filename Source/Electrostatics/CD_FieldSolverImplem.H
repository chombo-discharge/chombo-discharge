/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldSolverImplem.H
  @brief  Implementation of CD_FieldSolver.H
  @author Robert Marskar
*/

#ifndef CD_FieldSolverImplem_H
#define CD_FieldSolverImplem_H

// Our includes
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

Real
FieldSolver::getDielectricPermittivity(const RealVect& a_pos) const
{
  CH_TIME("FieldSolver::getDielectricPermittivity(RealVect a_pos)");

  const Vector<Dielectric>& dielectrics = m_computationalGeometry->getDielectrics();

  Real relPerm = 0.0;

  if (m_multifluidIndexSpace->numPhases() > 1 && dielectrics.size() > 0) {

    // Find the closest dielectric to position a_pos. This procedure might fail if
    // the user has two overlapping dielectrics, but those cases are not supported
    // anyways.
    Real minDist = std::numeric_limits<Real>::infinity();
    int  closest = 0;

    for (int i = 0; i < dielectrics.size(); i++) {
      const RefCountedPtr<BaseIF> func = dielectrics[i].getImplicitFunction();

      const Real curDist = func->value(a_pos);

      if (std::abs(curDist) <= std::abs(minDist)) {
        minDist = curDist;
        closest = i;
      }
    }

    relPerm = dielectrics[closest].getPermittivity(a_pos);
  }

  CH_assert(dielectrics.size() > 0);
  CH_assert(relPerm > 0.0);

  return std::max(0.0,a_pos[0]);
}

#include <CD_NamespaceFooter.H>

#endif
