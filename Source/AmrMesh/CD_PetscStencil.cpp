/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PetscStencil.cpp
  @brief  Implementation of CD_PetscStencil.H
  @author Robert Marskar
*/

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_assert.H>

// Our includes
#include <CD_PetscStencil.H>
#include <CD_NamespaceHeader.H>

PetscStencil::PetscStencil()
{
  this->clear();
}

PetscStencil::~PetscStencil()
{
  this->clear();
}

void
PetscStencil::clear()
{
  m_rows.resize(0);
  m_weights.resize(0);
}

int
PetscStencil::size() const noexcept
{
  CH_assert(m_rows.size() == m_weights.size());

  return m_rows.size();
}

PetscInt
PetscStencil::row(const int a_entry) const noexcept
{
  return m_rows[a_entry];
}

PetscScalar
PetscStencil::weight(const int a_entry) const noexcept
{
  return m_weights[a_entry];
}

void
PetscStencil::add(const PetscInt& a_row, const PetscScalar& a_weight) noexcept
{
  bool alreadyHere = false;

  const int N = this->size();

  for (int i = 0; i < N; i++) {
    if (m_rows[i] == a_row) {
      m_weights[i] += a_weight;

      alreadyHere = true;

      break;
    }
  }

  if (!alreadyHere) {
    m_rows.push_back(a_row);
    m_weights.push_back(a_weight);
  }
}

PetscStencil&
PetscStencil::operator+=(const PetscStencil& a_otherStencil) noexcept
{
  for (int i = 0; i < a_otherStencil.size(); i++) {
    this->add(a_otherStencil.row(i), a_otherStencil.weight(i));
  }

  return *this;
}

PetscStencil&
PetscStencil::operator*=(const Real& a_scalingFactor) noexcept
{
  const int N = this->size();

  for (int i = 0; i < N; i++) {
    m_weights[i] *= a_scalingFactor;
  }

  return *this;
}

#include <CD_NamespaceFooter.H>

#endif
