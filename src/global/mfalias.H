/*!
  @file mfalias.H
  @brief Declaration of factories that allow aliasing of LevelData<MFFluxFAB> and LevelData<MFBaseIVFAB>
  @author Robert Marskar
  @date Dec. 2017
  @todo Segregate header and implementation.
  @todo Do an aliasing for MFBaseIVFAB. 
*/

#ifndef _MFALIAS_
#define _MFALIAS_

#include <MFFluxFAB.H>
#include <MFCellFAB.H>

#include "MFBaseIVFAB.H"

/*!
  @brief Class for aliasing multifluid things
*/
class mfalias {
public:

  /*!
    @brief Alias MFCellFABs
  */
  static void aliasMF(LevelData<EBCellFAB>& a_alias, int a_phase, const LevelData<MFCellFAB>& a_input){
    mf_cell_alias_factory factory((LevelData<MFCellFAB>*)& a_input, a_phase);
    a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
  }

  /*!
    @brief Alias MFFluxFABs
  */
  static void aliasMF(LevelData<EBFluxFAB>& a_alias, int a_phase, const LevelData<MFFluxFAB>& a_input){
    mf_flux_alias_factory factory((LevelData<MFFluxFAB>*)& a_input, a_phase);
    a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
  }

    /*!
    @brief Alias MFBaseIVFABs
  */
  static void aliasMF(LevelData<BaseIVFAB<Real> >& a_alias, int a_phase, const LevelData<MFBaseIVFAB>& a_input){
    mf_baseivfab_alias_factory factory((LevelData<MFBaseIVFAB>*)& a_input, a_phase);
    a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
  }

  

private:
  
  /*!
    @brief Class for aliasing 
  */
  class mf_cell_alias_factory : public DataFactory<EBCellFAB> {
  public:

    /*!
      @brief Constructor
    */
    mf_cell_alias_factory(LevelData<MFCellFAB>* a_mf, int a_phase){
      m_mf    = a_mf;
      m_phase = a_phase;
    }

    /*!
      @brief Factory method
    */
    virtual EBCellFAB* create(const Box& box, int ncomps, const DataIndex& a_datInd) const{
      return (*m_mf)[a_datInd].getPhasePtr(m_phase);
    }

    /*!
      @brief Return false (data is not our responsibility)
    */
    virtual bool callDelete() const {
      return false;
    }

  private:

    /*!
      @brief Object
    */ 
    LevelData<MFCellFAB>* m_mf;

    /*!
      @brief Phase to alias
    */
    int m_phase;
  };
  
  /*!
    @brief Class for aliasing 
  */
  class mf_flux_alias_factory : public DataFactory<EBFluxFAB> {
  public:

    /*!
      @brief Constructor
    */
    mf_flux_alias_factory(LevelData<MFFluxFAB>* a_mf, int a_phase){
      m_mf    = a_mf;
      m_phase = a_phase;
    }

    /*!
      @brief Factory method
    */
    virtual EBFluxFAB* create(const Box& box, int ncomps, const DataIndex& a_datInd) const{
      return (*m_mf)[a_datInd].getPhasePtr(m_phase);
    }

    /*!
      @brief Return false (data is not our responsibility)
    */
    virtual bool callDelete() const {
      return false;
    }

  private:

    /*!
      @brief Object
    */ 
    LevelData<MFFluxFAB>* m_mf;

    /*!
      @brief Phase to alias
    */
    int m_phase;
  };

    /*!
    @brief Class for aliasing 
  */
  class mf_baseivfab_alias_factory : public DataFactory<BaseIVFAB<Real> > {
  public:

    /*!
      @brief Constructor
    */
    mf_baseivfab_alias_factory(LevelData<MFBaseIVFAB>* a_mf, int a_phase){
      m_mf    = a_mf;
      m_phase = a_phase;
    }

    /*!
      @brief Factory method
    */
    virtual BaseIVFAB<Real>* create(const Box& box, int ncomps, const DataIndex& a_datInd) const{
      return (*m_mf)[a_datInd].get_phase_ptr(m_phase);
    }

    /*!
      @brief Return false (data is not our responsibility)
    */
    virtual bool callDelete() const {
      return false;
    }

  private:

    /*!
      @brief Object
    */ 
    LevelData<MFBaseIVFAB>* m_mf;

    /*!
      @brief Phase to alias
    */
    int m_phase;
  };
};

#endif
