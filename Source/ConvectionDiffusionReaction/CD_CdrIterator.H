/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrIterator.H
  @brief  Declaration of an iterator class for CdrLayout (for iterating through subsets of solvers). 
  @author Robert Marskar
*/

#ifndef CD_CdrIterator_H
#define CD_CdrIterator_H

// Our includes
#include <CD_CdrLayout.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Iterator class for CdrLayout. This allows iteration through solvers (or subsets of solvers).
  @details CdrLayout is a class for holding sets of CdrSolvers on a realm and a phase. In the simplest case it is used as a container for keeping
  track of a single set of CdrSolvers, and it greatly cuts down on typing. This class allows iteration through one such layout. The user can specify
  to only iterate through subsets of CdrSolvers, e.g. only for charged species or neutral species, etc. 
*/
template <class T>
class CdrIterator
{
public:
  /*!
    @brief Enum class for iterating through subsets of species. 
  */
  enum class SpeciesIteration
  {
    All,
    Charged,
    Negative,
    Positive,
    Neutral
  };

  /*!
    @brief Disallowed weak construction. Always use the full constructor.
  */
  CdrIterator() = delete;

  /*!
    @brief Full constructor
    @param[in] a_layout CdrLayout to iterate through. 
    @param[in] a_mode   Iteration mode. 
  */
  CdrIterator(CdrLayout<T>& a_layout, const SpeciesIteration a_mode = SpeciesIteration::All);

  /*!
    @brief Destructor (does nothing).
  */
  virtual ~CdrIterator();

  /*!
    @brief Get number of solvers. 
  */
  virtual int
  getNumberOfSolvers();

  /*!
    @brief Get current solver number (i.e. index in CdrLayouts solvers). 
  */
  virtual int
  index() const;

  /*!
    @brief Reset iterator
  */
  virtual void
  reset();

  /*!
    @brief Ok or not
  */
  virtual bool
  ok() const;

  /*!
    @brief Increment. This is different for different iteration modes
  */
  virtual void
  operator++();

  /*!
    @brief Return the current solver. 
  */
  virtual RefCountedPtr<T>&
  operator()();

  /*!
    @brief Return the current species. 
  */
  virtual RefCountedPtr<CdrSpecies>&
  getSpecies();

protected:
  /*!
    @brief Which species to iterate over
  */
  SpeciesIteration m_mode;

  /*!
    @brief 
  */
  int m_isolver;

  /*!
    @brief Number of solvers
  */
  int m_numSolvers;

  /*!
    @brief Defined or not
  */
  bool m_isDefined;

  /*!
    @brief Cdr solvers
  */
  Vector<RefCountedPtr<T>> m_solvers;

  /*!
    @brief Cdr species. 
  */
  Vector<RefCountedPtr<CdrSpecies>> m_species;
};

#include <CD_NamespaceFooter.H>

#include <CD_CdrIteratorImplem.H>

#endif
