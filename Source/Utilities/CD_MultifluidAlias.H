/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MultifluidAlias.H
  @brief  Declaration of factories that allow aliasing of LevelData<MFFluxFAB> and LevelData<MFBaseIVFAB>
  @author Robert Marskar
*/

#ifndef CD_MultifluidAlias_H
#define CD_MultifluidAlias_H

// Chombo includes
#include <MFFluxFAB.H>
#include <MFCellFAB.H>

// Our includes
#include "MFBaseIVFAB.H"
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for aliasing multifluid things
*/
class MultifluidAlias {
public:

  /*!
    @brief Alias MFCellFABs
  */
  static void aliasMF(LevelData<EBCellFAB>& a_alias, int a_phase, const LevelData<MFCellFAB>& a_input){
    MfCellAliasFactory factory((LevelData<MFCellFAB>*)& a_input, a_phase);
    a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
  }

  /*!
    @brief Alias MFFluxFABs
  */
  static void aliasMF(LevelData<EBFluxFAB>& a_alias, int a_phase, const LevelData<MFFluxFAB>& a_input){
    MfFluxAliasFactory factory((LevelData<MFFluxFAB>*)& a_input, a_phase);
    a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
  }

  /*!
    @brief Alias MFBaseIVFABs
  */
  static void aliasMF(LevelData<BaseIVFAB<Real> >& a_alias, int a_phase, const LevelData<MFBaseIVFAB>& a_input){
    MfIVAliasFactory factory((LevelData<MFBaseIVFAB>*)& a_input, a_phase);
    a_alias.define(a_input.disjointBoxLayout(), a_input.nComp(), a_input.ghostVect(), factory);
  }

private:
  
  /*!
    @brief Class for aliasing 
  */
  class MfCellAliasFactory : public DataFactory<EBCellFAB> {
  public:

    /*!
      @brief Constructor
    */
    MfCellAliasFactory(LevelData<MFCellFAB>* a_mf, int a_phase){
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
  class MfFluxAliasFactory : public DataFactory<EBFluxFAB> {
  public:

    /*!
      @brief Constructor
    */
    MfFluxAliasFactory(LevelData<MFFluxFAB>* a_mf, int a_phase){
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
  class MfIVAliasFactory : public DataFactory<BaseIVFAB<Real> > {
  public:

    /*!
      @brief Constructor
    */
    MfIVAliasFactory(LevelData<MFBaseIVFAB>* a_mf, int a_phase){
      m_mf    = a_mf;
      m_phase = a_phase;
    }

    /*!
      @brief Factory method
    */
    virtual BaseIVFAB<Real>* create(const Box& box, int ncomps, const DataIndex& a_datInd) const{
      return (*m_mf)[a_datInd].getPhase_ptr(m_phase);
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

#include <CD_NamespaceFooter.H>

#endif
