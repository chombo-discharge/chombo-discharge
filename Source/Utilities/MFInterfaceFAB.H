/*!
  @file MFInterfaceFAB.H
  @brief Declaration of a class that allocates irregular data holders over the interface between two phases
  @author Robert Marskar
  @date Dec. 2017
  @todo Segregate implementation file
*/

#ifndef _MFINTERFACEBASEIVFAB_
#define _MFINTERFACEBASEIVFAB_

#include <CD_MultiFluidIndexSpace.H>
#include "MFLevelGrid.H"

#include <BoxLayoutData.H>
#include <BaseIVFAB.H>
#include <IntVectSet.H>

#include "CD_NamespaceHeader.H"

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

    const DisjointBoxLayout& dbl    = a_mflg.getGrids();
    const ProblemDomain& domain     = a_mflg.get_domain();
    const RefCountedPtr<MultiFluidIndexSpace>& MultiFluidIndexSpace = a_mflg.getMfIndexSpace();
    const int numPhases            = MultiFluidIndexSpace->numPhases();

    m_ivs = a_mflg.interfaceRegion(dbl.get(a_dit), a_dit);
    m_phase.resize(numPhases, NULL);
    
    for (int i = 0; i < numPhases; i++){
      const EBLevelGrid& eblg = a_mflg.getEBLevelGrid(i);
      const EBISLayout& ebisl = eblg.getEBISL();
      m_phase[i] = new BaseIVFAB<T>(m_ivs, ebisl[a_dit].getEBGraph(), 1);
    }
  }

  virtual const IntVectSet& get_ivs() const {
    return m_ivs;
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<T>& getIVFAB(const int a_phase){
    return *m_phase[a_phase];
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual const BaseIVFAB<T>& getIVFAB(const int a_phase) const {
    return *m_phase[a_phase];
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<T>* getPhasePtr(int a_phase){
    return m_phase[a_phase];
  }

  /*!
    @brief Get number of phases
  */
  virtual int numPhases(){
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

  IntVectSet m_ivs;
};
#include "CD_NamespaceFooter.H"

#endif
