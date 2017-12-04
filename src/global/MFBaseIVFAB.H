/*!
  @file MFBaseIVFAB.H
  @brief Declaration of a multiphase BaseIVFAB
  @author Robert Marskar
  @date Dec. 2017
*/

#ifndef _MFBASEIVFAB_
#define _MFBASEIVFAB_

#include "BoxLayoutData.H"
#include "BaseIVFAB.H"
#include "IntVectSet.H"
#include "EBGraph.H"

/*!
  @brief Multiphase BaseIVFAB
*/
class MFBaseIVFAB {
public:

  /*!
    @brief Constructor
  */
  MFBaseIVFAB(){
  }

  /*!
    @brief Full constructor. 
  */
  MFBaseIVFAB(const Vector<IntVectSet>& a_regions,
	  const Vector<EBGraph>&    a_phase_graphs,
	  const Vector<int>&        a_ncomp){
    CH_TIME("MFBaseIVFAB::MFBaseIVFAB");

    const int num_phases = a_regions.size();
    
    CH_assert(a_phase_graphs.size() == num_phases);
    CH_assert(a_ncomp.size() == num_phases);
    
    m_phase.resize(num_phases, NULL);
    for (int i = 0; i < num_phases; i++){
      m_phase[i] = new BaseIVFAB<Real>(a_regions[i], a_phase_graphs[i], a_ncomp[i]);
    }

  }

  /*!
    @brief Destructor
  */
  virtual ~MFBaseIVFAB(){
    for (int i = 0; i < m_phase.size(); i++){
      delete m_phase[i];
    }
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<Real>& get_ivfab(const int a_phase){
    CH_assert(a_phase < m_phase.size());
    return *m_phase[a_phase];
  }

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual const BaseIVFAB<Real>& get_ivfab(const int a_phase) const {
    CH_assert(a_phase < m_phase.size());
    return *m_phase[a_phase];
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
  virtual void setVal(Real a_value){
    for (int i = 0; i < m_phase.size(); i++){
      m_phase[i]->setVal(a_value);
    }
  }

  /*!
    @brief Copy function
  */
  virtual void copy(const Box& a_from_box,
		    const Interval& a_dst_interv,
		    const Box& a_to_box,
		    const MFBaseIVFAB& a_src,
		    const Interval& a_src_interv){
    for (int i = 0; i < m_phase.size(); i++){
      m_phase[i]->copy(a_from_box, a_dst_interv, a_to_box, *(a_src.m_phase[i]), a_src_interv);
    }
  }

  static int preAllocatable(){
    return 1;
  }

  /*!
    @brief Size function. Required for LevelData construction
  */
  virtual int size(const Box& R, const Interval& comps) const {
    int size = m_phase.size()*sizeof(int);
    for (int i = 0; i < m_phase.size(); ++i){
      size += m_phase[i]->size(R, comps);
    }

    return size;
  }

  /*!
    @brief Linearization function
  */
  virtual void linearOut(void* buf, const Box& R, const Interval& comps) const {
    int* buffer = (int*)buf;
    for (int i=0; i<m_phase.size(); ++i)
      {
	*buffer = m_phase[i]->size(R, comps);
	++buffer;
      }
    int* size = (int*)buf;
    unsigned char* ebbuffer = (unsigned char*)buffer;
    for (int i=0; i<m_phase.size(); ++i)
      {
	m_phase[i]->linearOut(ebbuffer, R, comps);
	ebbuffer+= size[i];
      }
  }

  /*!
    @brief Linearization function
  */
  virtual void linearIn(void* buf, const Box& R, const Interval& comps){
    int* size = (int*) buf;

    for (int i = 0; i < m_phase.size(); ++i){
      CH_assert(size[i] == m_phase[i]->size(R, comps));
    }

    unsigned char* ebbuffer  = (unsigned char*)(size+m_phase.size());
    for (int i=0; i<m_phase.size(); ++i){
      m_phase[i]->linearIn(ebbuffer, R, comps);
      ebbuffer += size[i];
    }
  }

protected:

  Vector<BaseIVFAB<Real>* > m_phase;
};

/*!
  @brief Factory class for creating MFBaseIVFABs
*/
class MFBaseIVFABFactory : public DataFactory<MFBaseIVFAB> {
public:

  /*!
    @brief Full constructor
  */
  MFBaseIVFABFactory(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_ncomp){
    CH_TIME("MFBaseIVFABFactory::MFBaseIVFABFactory");
    this->define(a_ebisl, a_ncomp);
  }

  /*!
    @brief Destructor
  */
  virtual ~MFBaseIVFABFactory(){
    CH_TIME("MFBaseIVFABFactory::~MFBaseIVFABFactory");
  }

  virtual void define(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_ncomp){
    CH_TIME("MFBaseIVFABFactory::define");
    CH_assert(a_ebisl.size() == a_ncomp.size());
    m_ebisl = a_ebisl;
    m_ncomp = a_ncomp;
  }

  /*!
    @brief Factory method
  */
  virtual MFBaseIVFAB* create(const Box& a_box, int a_ignored_argument, const DataIndex& a_dit) const {
    CH_TIME("MFBaseIVFABFactory::create");

    const int num_phases = m_ebisl.size();
    
    Vector<IntVectSet> ivs(num_phases);
    Vector<EBGraph> ebgraph(num_phases);
    
    for (int i = 0; i < num_phases; ++i){
      ivs[i]     = m_ebisl[i][a_dit].getIrregIVS(a_box);
      ebgraph[i] = m_ebisl[i][a_dit].getEBGraph();
    }

    return new MFBaseIVFAB(ivs, ebgraph, m_ncomp);
  }
									      

protected:

  /*!
    @brief EBIS layouts for different phases
  */
  Vector<EBISLayout> m_ebisl;

  /*!
    @brief Componenents
  */
  Vector<int> m_ncomp;

};
#endif
