/*!
  @file   rte_iterator.H
  @brief  Iterator class for rte_layout
  @author Robert Marskar
  @date   June 2016
*/

#ifndef _RTE_ITERATOR_
#define _RTE_ITERATOR_

#include "rte_layout.H"

/*!
  @brief Iterator class for rte_layout
*/
class rte_iterator {
public:  
  /*!
    @brief Empty constructor
  */
  rte_iterator(){
    CH_TIME("rte_iterator::rte_iterator");
    m_defined = false;
  }

  /*!
    @brief Full constructor
  */
  rte_iterator(rte_layout& a_layout){
    CH_TIME("rte_iterator::rte_iterator");
    
    m_solvers  = a_layout.get_solvers();
    m_species  = a_layout.get_species();
    m_num      = m_solvers.size();
    m_defined  = true;

    reset();
  }

  /*!
    @brief Destructor
  */
  virtual ~rte_iterator(){

  }

  /*!
    @brief Get number of solvers
  */
  virtual int num_solvers(){
    CH_assert(m_defined);
    return m_num;
  }

  /*!
    @brief Get current solver number
  */
  virtual int get_solver() const {
    CH_assert(m_defined);
    return m_isolver;
  }

  /*!
    @brief Reset iterator
  */
  virtual void reset(){
    CH_assert(m_defined);
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
    CH_assert(m_defined);
    m_isolver++;
  }

  /*!
    @brief Get solver
  */
  virtual RefCountedPtr<rte_solver>& operator() () {
    CH_assert(m_defined == true);
    CH_assert(m_isolver < m_num);
    return m_solvers[m_isolver];
  }

  /*!
    @brief Get current ion
  */
  virtual RefCountedPtr<rte_species>& get_species() {
    CH_assert(m_defined == true);
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
    @brief Defined or not
  */
  bool m_defined;

  /*!
    @brief Ion solvers
  */
  Vector<RefCountedPtr<rte_solver> > m_solvers;

  /*!
    @brief Ions
  */
  Vector<RefCountedPtr<rte_species> > m_species;

};
#endif
