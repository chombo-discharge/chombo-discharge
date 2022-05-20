/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MeshODESolverImplem.H
  @brief  Implementation of CD_MeshODESolver.H
  @author Robert Marskar
*/

#ifndef CD_MeshODESolverImplem_H
#define CD_MeshODESolverImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MeshODESolver.H>

template <size_t N>
MeshODESolver<T>::MeshODESolver()
{
  CH_TIME("MeshODESolver::MeshODESolver");

  m_className = "MeshODESolver";
  m_name      = "MeshODESolver";
}

template <size_t N>
MeshODESolver<T>::~MeshODESolver()
{}

template <size_t N>
MeshODESolver<T>::

#endif
