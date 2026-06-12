/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MultifluidAlias.H
  @brief  Declaration of factories that allow aliasing of LevelData<MFFluxFAB> and LevelData<MFBaseIVFAB>
  @author Robert Marskar
*/

#ifndef CD_MULTIFLUIDALIAS_H
#define CD_MULTIFLUIDALIAS_H

// Chombo includes
#include <MFFluxFAB.H>
#include <MFCellFAB.H>

// Our includes
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Class with static public interface for aliasing multifluid data into single-phase data.
  @details To use this class, the user will call aliasMF to fetch pointers from multifluid data and have them
  point to their single-fluid counterpart.
*/
class MultifluidAlias
{
public:
  /**
    @brief Alias MFCellFABs
    @param[out] a_alias Single-phase data
    @param[in]  a_phase Phase
    @param[in]  a_input Multi-phase data
  */
  static void
  aliasMF(LevelData<EBCellFAB>& a_alias, const int a_phase, const LevelData<MFCellFAB>& a_input);

  /**
    @brief Alias MFFluxFABs
    @param[out] a_alias Single-phase data
    @param[in]  a_phase Phase
    @param[in]  a_input Multi-phase data
  */
  static void
  aliasMF(LevelData<EBFluxFAB>& a_alias, const int a_phase, const LevelData<MFFluxFAB>& a_input);

  /**
    @brief Alias MFBaseIVFABs
    @param[out] a_alias Single-phase data
    @param[in]  a_phase Phase
    @param[in]  a_input Multi-phase data
  */
  static void
  aliasMF(LevelData<BaseIVFAB<Real>>& a_alias, const int a_phase, const LevelData<MFBaseIVFAB>& a_input);

private:
  /**
    @brief Class for aliasing. This is used under the hood in MultifluidAlias::aliasMF
  */
  class MfCellAliasFactory : public DataFactory<EBCellFAB>
  {
  public:
    /**
      @brief Don't see what a weak constructor would do here.
    */
    MfCellAliasFactory() = delete;

    /**
      @brief Constructor
      @param[in] a_mf    Multifluid cell data
      @param[in] a_phase Phase to alias
    */
    MfCellAliasFactory(LevelData<MFCellFAB>* a_mf, const int a_phase);

    /**
      @brief Factory method
      @param[in] box      Grid box (unused)
      @param[in] ncomps   Number of components (unused)
      @param[in] a_datInd Data index
      @return Aliased EBCellFAB pointer (not owned by caller).
    */
    virtual EBCellFAB*
    create(const Box& box, int ncomps, const DataIndex& a_datInd) const override;

    /**
      @brief Return false (data is not our responsibility)
      @return False.
    */
    virtual bool
    callDelete() const override;

  private:
    /**
      @brief Object
    */
    LevelData<MFCellFAB>* m_mf;

    /**
      @brief Phase to alias
    */
    int m_phase;
  };

  /**
    @brief Class for aliasing. This is used under the hood in MultifluidAlias::aliasMF
  */
  class MfFluxAliasFactory : public DataFactory<EBFluxFAB>
  {
  public:
    /**
      @brief Don't see what a weak constructor would do here.
    */
    MfFluxAliasFactory() = delete;

    /**
      @brief Constructor
      @param[in] a_mf    Multifluid flux data
      @param[in] a_phase Phase to alias
    */
    MfFluxAliasFactory(LevelData<MFFluxFAB>* a_mf, const int a_phase);

    /**
      @brief Factory method
      @param[in] box      Grid box (unused)
      @param[in] ncomps   Number of components (unused)
      @param[in] a_datInd Data index
      @return Aliased EBFluxFAB pointer (not owned by caller).
    */
    virtual EBFluxFAB*
    create(const Box& box, int ncomps, const DataIndex& a_datInd) const override;

    /**
      @brief Return false (data is not our responsibility)
      @return False.
    */
    virtual bool
    callDelete() const override;

  private:
    /**
      @brief Object
    */
    LevelData<MFFluxFAB>* m_mf;

    /**
      @brief Phase to alias
    */
    int m_phase;
  };

  /**
    @brief Class for aliasing. This is used under the hood in MultifluidAlias::aliasMF
  */
  class MfIVAliasFactory : public DataFactory<BaseIVFAB<Real>>
  {
  public:
    /**
      @brief Don't see what a weak constructor would do here.
    */
    MfIVAliasFactory() = delete;

    /**
      @brief Constructor
      @param[in] a_mf    Multifluid irregular data
      @param[in] a_phase Phase to alias
    */
    MfIVAliasFactory(LevelData<MFBaseIVFAB>* a_mf, const int a_phase);

    /**
      @brief Factory method
      @param[in] box      Grid box (unused)
      @param[in] ncomps   Number of components (unused)
      @param[in] a_datInd Data index
      @return Aliased BaseIVFAB<Real> pointer (not owned by caller).
    */
    virtual BaseIVFAB<Real>*
    create(const Box& box, int ncomps, const DataIndex& a_datInd) const override;

    /**
      @brief Return false (data is not our responsibility)
      @return False.
    */
    virtual bool
    callDelete() const override;

  private:
    /**
      @brief Object
    */
    LevelData<MFBaseIVFAB>* m_mf;

    /**
      @brief Phase to alias
    */
    int m_phase;
  };
};

#include <CD_NamespaceFooter.H>

#endif
