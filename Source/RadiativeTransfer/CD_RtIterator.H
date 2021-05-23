/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtIterator.H
  @brief  Iterator class for RtLayout
  @author Robert Marskar
*/

#ifndef CD_RtIterator_H
#define CD_RtIterator_H

// Our includes
#include <CD_RtLayout.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Iterator class for RtLayout
*/
template <class T>
class RtIterator {
public:  
  /*!
    @brief Empty constructor
  */
  RtIterator() = delete;

  /*!
    @brief Full constructor
  */
  RtIterator(RtLayout<T>& a_layout){
    CH_TIME("RtIterator::RtIterator");
    
    m_solvers  = a_layout.getSolvers();
    m_species  = a_layout.getSpecies();
    m_num      = m_solvers.size();

    reset();
  }

  /*!
    @brief Destructor
  */
  virtual ~RtIterator(){

  }

  /*!
    @brief Get number of solvers
  */
  virtual int getNumberOfSolvers(){
    return m_num;
  }

  /*!
    @brief Get current solver number
  */
  virtual int index() const {
    return m_isolver;
  }

  /*!
    @brief Reset iterator
  */
  virtual void reset(){
    m_isolver = 0;
  }


  /*!
    @brief Ok or not
  */
  virtual bool ok(){
    return (m_isolver < m_num);
  }

  /*!
    @brief Increment
  */
  virtual void operator++(){
    m_isolver++;
  }

  /*!
    @brief Get solver
  */
  virtual RefCountedPtr<T>& operator() () {
   CH_assert(m_isolver < m_num);
    return m_solvers[m_isolver];
  }

  /*!
    @brief Get current ion
  */
  virtual RefCountedPtr<RtSpecies>& getSpecies() {
    CH_assert(m_isolver < m_num);
    return m_species[m_isolver];
  }
  
protected:
  
  /*!
    @brief 
  */
  int m_isolver;

  /*!
    @brief Number of solvers
  */
  int m_num;;

  /*!
    @brief Ion solvers
  */
  Vector<RefCountedPtr<T> > m_solvers;

  /*!
    @brief Ions
  */
  Vector<RefCountedPtr<RtSpecies> > m_species;

};

#include <CD_NamespaceFooter.H>

#endif
