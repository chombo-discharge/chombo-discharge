/*!
  @file MFInterfaceFAB.H
  @brief Declaration of a class that allocates irregular data holders over the interface between two phases
  @author Robert Marskar
  @date Dec. 2017
  @todo Segregate implementation file
*/

#ifndef _MFINTERFACEBASEIVFAB_
#define _MFINTERFACEBASEIVFAB_

#include "mfis.H"
#include "MFLevelGrid.H"

#include <BoxLayoutData.H>
#include <BaseIVFAB.H>
#include <IntVectSet.H>

/*!
  @brief Class for holding BaseIVFAB<T> on the interface between two phases. 
*/
template <class T> class MFInterfaceFAB {
public:

  /*!
    @brief Constructor. Must subsequently call define
  */
  MFInterfaceFAB(){
  }

  /*!
    @brief Full constructor. 
  */
  MFInterfaceFAB(const MFLevelGrid& a_mflg, const DataIndex& a_dit){
    this->define(a_mflg, a_dit);
  }

  /*!
    @brief Destructor
  */
  virtual ~MFInterfaceFAB(){
    for (int i = 0; i < m_phase.size(); i++){
      delete m_phase[i];
    }
  }

  /*!
    @brief Define function
  */
  virtual void define(const MFLevelGrid& a_mflg, const DataIndex& a_dit){

    const DisjointBoxLayout& dbl    = a_mflg.get_grids();
    const ProblemDomain& domain     = a_mflg.get_domain();
    const RefCountedPtr<mfis>& mfis = a_mflg.get_mfis();
    const IntVectSet ivs            = mfis->interface_region(domain) & dbl.get(a_dit);
    const int num_phases            = mfis->num_phases();
    
    m_phase.resize(num_phases, NULL);
    
    for (int i = 0; i < num_phases; i++){
      const EBLevelGrid& eblg = a_mflg.get_eblg(i);
      const EBISLayout& ebisl = eblg.getEBISL();
      m_phase[i] = new BaseIVFAB<T>(ivs, ebisl[a_dit].getEBGraph(), 1);
    }
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<T>& get_ivfab(const int a_phase){
    return *m_phase[a_phase];
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual const BaseIVFAB<T>& get_ivfab(const int a_phase) const {
    return *m_phase[a_phase];
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<T>* get_phase_ptr(int a_phase){
    return m_phase[a_phase];
  }

  /*!
    @brief Get number of phases
  */
  virtual int num_phases(){
    return m_phase.size();
  }

  /*!
    @brief Set value
  */
  virtual void setVal(const T& a_value){
    for (int i = 0; i < m_phase.size(); i++){
      m_phase[i]->setVal(a_value);
    }
  }


protected:

  Vector<BaseIVFAB<T>* > m_phase;
};

#endif
