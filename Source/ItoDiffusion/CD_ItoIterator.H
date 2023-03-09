/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoIterator.H
  @brief  Declaration of an iterator class for ItoLayout
  @author Robert Marskar
*/

#ifndef CD_ItoIterator_H
#define CD_ItoIterator_H

#include <CD_ItoLayout.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief "Iterator" class for going through solvers in an ItoLayout. 
  @details This class does not have any fundamental usefulness -- it simply cuts down on typing if one wants to iterate over a full set or a subset of solvers in an ItoLayout.
*/
template <class T>
class ItoIterator
{
public:
  /*!
    @brief Enum which specifies which species to iterate through. 
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
    @brief Disallowed constructor -- use the full constructor.
  */
  ItoIterator() = delete;

  /*!
    @brief Full constructor. Allows iteration through all solvers in the layout and restriction to a specific subset of solvers.
    @param[inout] a_layout        The layout to iterate through
    @param[in]    a_iterationMode Iteration mode, e.g. all solvers, only charged species etc. 
  */
  ItoIterator(ItoLayout<T>& a_layout, const SpeciesIteration a_mode = SpeciesIteration::All);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~ItoIterator();

  /*!
    @brief Get number of solvers. 
    @return Returns total number of solvers (ignoring subsets)
  */
  virtual int
  getNumberOfSolvers();

  /*!
    @brief Get current solver number.
    @return This returns the index in the ItoLayout solvers (which stores solvers as Vector<RefCountedPtr<ItoSolver> >
  */
  virtual int
  index() const;

  /*!
    @brief Reset iterator to the first solver. 
  */
  virtual void
  reset();

  /*!
    @brief Ok or not
    @return Returns true if there is a next solver and false otherwise. 
  */
  virtual bool
  ok();

  /*!
    @brief Increment. This is different for different iteration modes.
  */
  virtual void
  operator++();

  /*!
    @brief Get current solver
  */
  virtual RefCountedPtr<T>&
  operator()();

  /*!
    @brief Get current species.
  */
  virtual RefCountedPtr<ItoSpecies>&
  getSpecies();

protected:
  /*!
    @brief Which species to iterate over
  */
  SpeciesIteration m_mode;

  /*!
    @brief Current solver pointed to by the iterator. 
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
    @brief Solvers
  */
  Vector<RefCountedPtr<T>> m_solvers;

  /*!
    @brief Solver species.
  */
  Vector<RefCountedPtr<ItoSpecies>> m_species;
};

#include <CD_NamespaceFooter.H>

#include <CD_ItoIteratorImplem.H>

#endif
