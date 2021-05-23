/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoIterator.H
  @brief  Iterator class for ItoLayout
  @author Robert Marskar
*/

#ifndef CD_ItoIterator_H
#define CD_ItoIterator_H

#include <CD_ItoLayout.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Iterator class for ItoLayout
*/
template <class T>
class ItoIterator {
public:

  /*!
    @brief Silly enum for iterating through selected species
  */
  enum class SpeciesIteration{
    All,
    Charged,
    Negative,
    Positive,
    Neutral
  };
  
  /*!
    @brief Empty constructor
  */
  ItoIterator();

  /*!
    @brief Full constructor
  */
  ItoIterator(ItoLayout<T>& a_layout, const SpeciesIteration a_mode = SpeciesIteration::All);

  /*!
    @brief Destructor
  */
  virtual ~ItoIterator();

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
  virtual RefCountedPtr<ItoSpecies>& getSpecies();
  
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
  Vector<RefCountedPtr<ItoSpecies> > m_species;

};

#include <CD_NamespaceFooter.H>

#include <CD_ItoIteratorImplem.H>

#endif
