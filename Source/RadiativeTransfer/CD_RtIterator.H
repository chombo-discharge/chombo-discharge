/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RtIterator.H
  @brief  Iterator class for RtLayout
  @author Robert Marskar
*/

#ifndef CD_RTITERATOR_H
#define CD_RTITERATOR_H

// Our includes
#include <CD_RtLayout.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Iterator class for RtLayout
*/
template <class T>
class RtIterator
{
public:
  /**
    @brief Empty constructor
  */
  RtIterator() = delete;

  /**
    @brief Full constructor
  */
  RtIterator(RtLayout<T>& a_layout);

  /**
    @brief Destructor (does nothing)
  */
  virtual ~RtIterator();

  /**
    @brief Get number of solvers
    @return Number of solvers
  */
  virtual int
  getNumberOfSolvers();

  /**
    @brief Get current solver number/index
    @return Return value
  */
  virtual int
  index() const;

  /**
    @brief Reset iterator to first solver.
  */
  virtual void
  reset();

  /**
    @brief Check if we can cycle further through the solvers. 
    @return Return value
  */
  virtual bool
  ok();

  /**
    @brief Incrementation operator. Cycles to next solver.
  */
  virtual void
  operator++();

  /**
    @brief Get the solver pointed to by the iterator.
    @return Return value
  */
  virtual RefCountedPtr<T>&
  operator()();

  /**
    @brief Get current species pointed to by the iterator.
    @return Species
  */
  virtual RefCountedPtr<RtSpecies>&
  getSpecies();

protected:
  /**
    @brief Solver currently pointed to.
  */
  int m_isolver;

  /**
    @brief Number of solvers
  */
  int m_num;

  /**
    @brief Ion solvers
  */
  Vector<RefCountedPtr<T>> m_solvers;

  /**
    @brief Ions
  */
  Vector<RefCountedPtr<RtSpecies>> m_species;
};

#include <CD_NamespaceFooter.H>

#include <CD_RtIteratorImplem.H>
#endif
