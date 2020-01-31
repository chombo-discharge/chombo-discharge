#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFAliasFactory.H"
#include "NamespaceHeader.H"

MFAliasFactory::MFAliasFactory(LevelData<MFCellFAB>* a_mf,
  int a_phase):m_mf(a_mf), m_phase(a_phase)
{
}

EBCellFAB* MFAliasFactory::create(const Box& box, int ncomps,
                                  const DataIndex& a_datInd) const
{
  return (*m_mf)[a_datInd].getPhasePtr(m_phase);
}


RegularMFAliasFactory::RegularMFAliasFactory(LevelData<MFCellFAB>* a_mf,
                                             int a_phase):m_mf(a_mf), m_phase(a_phase)
{
}

FArrayBox* RegularMFAliasFactory::create(const Box& box, int ncomps,
                                         const DataIndex& a_datInd) const
{
  return &((*m_mf)[a_datInd].getPhase(m_phase).getFArrayBox());
}


void aliasMF(Vector<LevelData<EBCellFAB>* >& a_phases, int numPhases,
             const LevelData<MFCellFAB>& a_input)
{
  CH_assert(a_phases.size() == numPhases);
  for (int i=0; i<numPhases; i++)
  {
    CH_assert(a_phases[i] != NULL);
    MFAliasFactory  factory((LevelData<MFCellFAB>*)&a_input, i);
    a_phases[i]->define(a_input.disjointBoxLayout(), a_input.nComp(),
                       a_input.ghostVect(), factory);
  }
}

void aliasMF(LevelData<EBCellFAB>& alias,  int phase,
             const LevelData<MFCellFAB>& a_input)
{
  MFAliasFactory  factory((LevelData<MFCellFAB>*)&a_input, phase);
  alias.define(a_input.disjointBoxLayout(),
               a_input.nComp(),a_input.ghostVect(), factory);
}

void aliasMF(Vector<LevelData<EBCellFAB>* >& a_output, int whichPhase,
             const Vector<LevelData<MFCellFAB>* >& a_input)
{
  CH_assert(a_output.size() == a_input.size());
  for (int level=0; level<a_output.size(); ++level)
    {
      LevelData<EBCellFAB>& eb = *(a_output[level]);
      const LevelData<MFCellFAB>& mf = *(a_input[level]);
      aliasMF(eb, whichPhase, mf);
    }
}

void aliasRegularMF(LevelData<FArrayBox>& alias,  int phase,
                    const LevelData<MFCellFAB>& a_input)
{
  RegularMFAliasFactory  factory((LevelData<MFCellFAB>*)&a_input, phase);
  alias.define(a_input.disjointBoxLayout(),
               a_input.nComp(),a_input.ghostVect(), factory);
}
#include "NamespaceFooter.H"
