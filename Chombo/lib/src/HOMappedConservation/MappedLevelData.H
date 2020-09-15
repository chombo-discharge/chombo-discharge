#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MAPPEDLEVELDATA_H_
#define _MAPPEDLEVELDATA_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Maintenance of <U> and <JU>
 *
 *//*+*************************************************************************/

#include "Copier.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "NamespaceHeader.H"

//--Forward declarations

class LevelGridMetrics;
class MultiBlockLevelExchange;


/*******************************************************************************
 */
/// Solution data on a level for mapped grids
/**
 *  This class aids the user with maintaining \<U\> and \<JU\> on the level.
 *  These quantities are frequently updated during the AMR process requiring
 *  new compuations of <\U\> and exchanging or filling ghosts.  This routine
 *  tracks the number of ghosts and performs the appropriate operations to
 *  return either <\U\> or \<JU\> with the requested number of valid or
 *  invalid ghost cells.
 *
 *  \note
 *  <ul>
 *    <li> Beware that if getU() and getJUNew() are used as used as arguments
 *         to a function, the order in which they are invoked is not defined
 *         and getU() may have side-effects on \<JU\>
 *  </ul>
 *
 *//*+*************************************************************************/

class MappedLevelData
{

public:

/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  MappedLevelData();

  /// Destructor
  ~MappedLevelData();

private:

  // Copy and assignment not allowed
  MappedLevelData(const MappedLevelData&);
  MappedLevelData& operator=(const MappedLevelData&);


/*==============================================================================
 * Member functions
 *============================================================================*/

public:

  ///
  void define(LevelGridMetrics *const           a_LGMP,
              LevelData<FArrayBox> *const       a_ULevP,
              LevelData<FArrayBox> *const       a_JUNewLevP,
              LevelData<FArrayBox> *const       a_JUOldLevP,
              const LevelData<FArrayBox> *const a_CrULevP,
              const int                         a_numComp,
              const Interval&                   a_velocityIntv,
              const int                         a_numGhostU,
              const int                         a_numGhostJU         = 1,
              const bool                        a_trimUExchanger     = false,
              const bool                        a_trimJUExchanger    = true,
              const bool                        a_useJUExchangerForU = true);

  ///
  void reshapeJUNew();

  ///
  void reshapeJUOld();

  ///
  void reshapeU();

  ///
  const LevelData<FArrayBox>& getU(const int a_numValidGhost = 0,
                                   const int a_numInvalidGhost = 0) const;

  ///
  LevelData<FArrayBox>& getU(const int a_numValidGhost = 0,
                             const int a_numInvalidGhost = 0);

  ///
  const LevelData<FArrayBox>& getJU(const int a_numValidGhost = 0,
                                    const int a_numInvalidGhost = 0) const;

  ///
  LevelData<FArrayBox>& getJU(const int a_numValidGhost = 0,
                              const int a_numInvalidGhost = 0);

  /// Return JUOld (always raw)
  LevelData<FArrayBox>& getJUOld();

  /// Return pointer to U data with no computation
  const LevelData<FArrayBox>& rawU() const;
  LevelData<FArrayBox>& rawU();

  /// Return pointer to JU data with no computation
  const LevelData<FArrayBox>& rawJU() const;
  LevelData<FArrayBox>& rawJU();

  /// Invalidate all of \<U\>
  void invalidateU(const bool a_this2FnInterpDataInvalid = true);

  ///
  void invalidateJUGhosts();

  /// Invalidate all ghost data and all of \<U\>
  void invalidate(const bool a_this2FnInterpDataInvalid = true);

  /// (assumes the same layout)
  void copyJUNewToOld();

  /// (assumes different layouts)
  void copyJUOldToNew();

  ///
  void swapJUNewAndOld();

  /// Return the exchange copier for U
  const Copier& getUExchangeCopier() const;

  /// Return the exchange copier for JU
  const Copier& getJUExchangeCopier() const;

//--Used to reference objects in this class, even before defined

  /// Return the exchange copier for U (before defined)
  const Copier& refUExchangeCopier() const;

  /// Return the exchange copier for JU (before defined)
  const Copier& refJUExchangeCopier() const;

protected:

  ///
  void deallocate();

  ///
  void fillInvalidGhosts();

  ///
  bool hasCoarser() const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:
  LevelGridMetrics *m_LGMP;

  LevelData<FArrayBox> *m_ULevP;
  LevelData<FArrayBox> *m_JULevP[2];

  const LevelData<FArrayBox> *m_CrULevP;

  Copier m_UExchangeCopier;
  Copier m_JUExchangeCopier;

  Interval m_velocityIntv;            ///< The velocity interval in the
                                      ///< components of \<U\> (or any other
                                      ///< space vectors but only 1 is
                                      ///< supported)

  int m_invalidGhostStatusU;          ///< The number of invalid ghosts (those
                                      ///< that do not overlap valid cells on
                                      ///< this level) in which /<U/> has been
                                      ///< interpolated.
  int m_invalidGhostStatusJU;

  int m_validStatusU;
  int m_validStatusJU;

  int m_validCFLayerStatusU;          ///< If > 0, the number of layers of valid
                                      ///< cells near the coarse-fine interface
                                      ///< that have not been filled
  int m_extraBlockStatusU;            ///< The number of extra-block ghost cells
                                      ///< in which \<U\> has been exchanged.
                                      ///< These are really valid ghosts but
                                      ///< are tracked separately.

  int m_numComp;

  int m_numGhostU;
  int m_numGhostJU;

  int m_indexJUOld;
  int m_indexJUNew;

  int m_allocU;
  int m_allocJUNew;
  int m_allocJUOld;

  bool m_trimUExchanger;
  bool m_trimJUExchanger;
  bool m_useJUExchangerForU;

  bool m_haveValidCrFnU;

  bool m_defined;
};


/*******************************************************************************
 *
 * Class MappedLevelData: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Return pointer to U data with no computation
/** Should only be used for accessing memory.  Assume data is invalid.
 *//*-----------------------------------------------------------------*/

inline const LevelData<FArrayBox>&
MappedLevelData::rawU() const
{
  CH_assert(m_defined);
  return *m_ULevP;
}

inline LevelData<FArrayBox>&
MappedLevelData::rawU()
{
  CH_assert(m_defined);
  return *m_ULevP;
}

/*--------------------------------------------------------------------*/
//  Return pointer to JU data with no computation
/** Should only be used for accessing memory.  Assume data is invalid.
 *//*-----------------------------------------------------------------*/

inline const LevelData<FArrayBox>&
MappedLevelData::rawJU() const
{
  CH_assert(m_defined);
  return *(m_JULevP[m_indexJUNew]);
}

inline LevelData<FArrayBox>&
MappedLevelData::rawJU()
{
  CH_assert(m_defined);
  return *(m_JULevP[m_indexJUNew]);
}

/*--------------------------------------------------------------------*/
//  Return JUOld
/** You can't set ghosts anymore so make sure JUNew has ghosts set
 *  properly before you make it old.
 *//*-----------------------------------------------------------------*/

inline LevelData<FArrayBox>&
MappedLevelData::getJUOld()
{
  CH_assert(m_defined);
  return *(m_JULevP[m_indexJUOld]);
}

/*--------------------------------------------------------------------*/
//  Return the exchange copier for U
/*--------------------------------------------------------------------*/

inline const Copier&
MappedLevelData::getUExchangeCopier() const
{
  CH_assert(m_defined);
  return m_UExchangeCopier;
}

/*--------------------------------------------------------------------*/
//  Return the exchange copier for JU
/*--------------------------------------------------------------------*/

inline const Copier&
MappedLevelData::getJUExchangeCopier() const
{
  CH_assert(m_defined);
  return m_JUExchangeCopier;
}

/*--------------------------------------------------------------------*/
//  Return the exchange copier for U (before defined)
/** Should only be used for accessing memory.  Assume not defined.
 *//*-----------------------------------------------------------------*/

inline const Copier&
MappedLevelData::refUExchangeCopier() const
{
  return m_UExchangeCopier;
}

/*--------------------------------------------------------------------*/
//  Return the exchange copier for JU (before defined)
/** Should only be used for accessing memory.  Assume not defined.
 *//*-----------------------------------------------------------------*/

inline const Copier&
MappedLevelData::refJUExchangeCopier() const
{
  return m_JUExchangeCopier;
}

/*--------------------------------------------------------------------*/

inline bool
MappedLevelData::hasCoarser() const
{
  return (m_CrULevP != 0);
}

#include "NamespaceFooter.H"
#endif
