/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrIterator.H
  @brief  Iterator class for CdrLayout
  @author Robert Marskar
*/

#ifndef CD_CdrIterator_H
#define CD_CdrIterator_H

// Our includes
#include <CD_CdrLayout.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Iterator class for CdrLayout
*/
template <class T>
class CdrIterator {
public:

  /*!
    @brief Silly enum for iterating through selected species
  */
  enum class SpeciesIteration{
    all,
    charged,
    negative,
    positive,
    neutral
  };
  
  /*!
    @brief Empty constructor
  */
  CdrIterator();

  /*!
    @brief Full constructor
  */
  CdrIterator(CdrLayout<T>& a_layout, const SpeciesIteration a_mode = SpeciesIteration::all);

  /*!
    @brief Destructor
  */
  virtual ~CdrIterator();

  /*!
    @brief Get number of solvers. 
  */
  virtual int getNumberOfSolvers();

  /*!
    @brief Get current solver number
  */
  virtual int index() const;

  /*!
    @brief Reset iterator
  */
  virtual void reset();

  /*!
    @brief Ok or not
  */
  virtual bool ok();
  
  /*!
    @brief Increment. This is different for different iteration modes
  */
  virtual void operator++();

  /*!
    @brief Get solver
  */
  virtual RefCountedPtr<T>& operator() ();

  /*!
    @brief Get current ion
  */
  virtual RefCountedPtr<CdrSpecies>& getSpecies();
  
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
  int m_num;;

  /*!
    @brief Defined or not
  */
  bool m_defined;

  /*!
    @brief Ion solvers
  */
  Vector<RefCountedPtr<T> > m_solvers;

  /*!
    @brief Ions
  */
  Vector<RefCountedPtr<CdrSpecies> > m_species;
};

#include <CD_NamespaceFooter.H>

#include <CD_CdrIteratorImplem.H>

#endif
