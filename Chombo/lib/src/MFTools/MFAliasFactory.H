#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFALIASFACTORY_H_
#define _MFALIASFACTORY_H_

#include "EBCellFAB.H"
#include "BoxLayoutData.H"
#include "MFCellFAB.H"
#include "NamespaceHeader.H"

class MFAliasFactory : public DataFactory<EBCellFAB>
{
public:
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual EBCellFAB* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

  virtual bool callDelete() const
  {
    return false;
  }

  MFAliasFactory(LevelData<MFCellFAB>* a_mf,
                 int a_phase);

private:

  LevelData<MFCellFAB>* m_mf;
  int m_phase;
};

class RegularMFAliasFactory : public DataFactory<FArrayBox>
{
public:
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual FArrayBox* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

  virtual bool callDelete() const
  {
    return false;
  }

  RegularMFAliasFactory(LevelData<MFCellFAB>* a_mf,
                        int a_phase);

private:

  LevelData<MFCellFAB>* m_mf;
  int m_phase;
};

void aliasMF(Vector<LevelData<EBCellFAB>* >& a_phases, int a_numPhases,
             const LevelData<MFCellFAB>& a_input);

void aliasMF(LevelData<EBCellFAB>& alias,  int phase,
             const LevelData<MFCellFAB>& a_input);

void aliasMF(Vector<LevelData<EBCellFAB>* >& a_output, int whichPhase,
             const Vector<LevelData<MFCellFAB>* >& a_input);

void aliasRegularMF(LevelData<FArrayBox>& alias,  int phase,
                    const LevelData<MFCellFAB>& a_input);

#include "NamespaceFooter.H"
#endif
