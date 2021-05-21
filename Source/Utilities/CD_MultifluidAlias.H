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
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for aliasing multifluid things
*/
class MultifluidAlias {
public:

  /*!
    @brief Alias MFCellFABs
  */
  static void aliasMF(LevelData<EBCellFAB>& a_alias, const int a_phase, const LevelData<MFCellFAB>& a_input);

  /*!
    @brief Alias MFFluxFABs
  */
  static void aliasMF(LevelData<EBFluxFAB>& a_alias, const int a_phase, const LevelData<MFFluxFAB>& a_input);

  /*!
    @brief Alias MFBaseIVFABs
  */
  static void aliasMF(LevelData<BaseIVFAB<Real> >& a_alias, const int a_phase, const LevelData<MFBaseIVFAB>& a_input);

private:
  
  /*!
    @brief Class for aliasing 
  */
  class MfCellAliasFactory : public DataFactory<EBCellFAB> {
  public:

    /*!
      @brief Don't see what a weak constructor would do here.
    */
    MfCellAliasFactory() = delete;

    /*!
      @brief Constructor
    */
    MfCellAliasFactory(LevelData<MFCellFAB>* a_mf, const int a_phase);

    /*!
      @brief Factory method
    */
    virtual EBCellFAB* create(const Box& box, int ncomps, const DataIndex& a_datInd) const override;

    /*!
      @brief Return false (data is not our responsibility)
    */
    virtual bool callDelete() const override;

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
      @brief Don't see what a weak constructor would do here.
    */
    MfFluxAliasFactory() = delete;

    /*!
      @brief Constructor
    */
    MfFluxAliasFactory(LevelData<MFFluxFAB>* a_mf, const int a_phase);

    /*!
      @brief Factory method
    */
    virtual EBFluxFAB* create(const Box& box, int ncomps, const DataIndex& a_datInd) const override;

    /*!
      @brief Return false (data is not our responsibility)
    */
    virtual bool callDelete() const override;

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
      @brief Don't see what a weak constructor would do here.
    */
    MfIVAliasFactory() = delete;

    /*!
      @brief Constructor
    */
    MfIVAliasFactory(LevelData<MFBaseIVFAB>* a_mf, const int a_phase);

    /*!
      @brief Factory method
    */
    virtual BaseIVFAB<Real>* create(const Box& box, int ncomps, const DataIndex& a_datInd) const override;

    /*!
      @brief Return false (data is not our responsibility)
    */
    virtual bool callDelete() const override;

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
