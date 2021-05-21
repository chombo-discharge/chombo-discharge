/*!
  @file   ito_iterator.H
  @brief  Iterator class for ito_layout
  @author Robert Marskar
  @date   June Dec. 2017
*/

#ifndef _ITO_ITERATOR_
#define _ITO_ITERATOR_

#include "ito_layout.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Iterator class for ito_layout
*/
template <class T>
class ito_iterator {
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
  ito_iterator();

  /*!
    @brief Full constructor
  */
  ito_iterator(ito_layout<T>& a_layout, const SpeciesIteration a_mode = SpeciesIteration::all);

  /*!
    @brief Destructor
  */
  virtual ~ito_iterator();

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
#include "CD_NamespaceFooter.H"

#include "ito_iteratorI.H"

#endif
