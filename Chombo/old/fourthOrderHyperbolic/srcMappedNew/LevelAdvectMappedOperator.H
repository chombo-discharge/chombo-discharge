#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELADVECTMAPPEDOPERATOR_H_
#define _LEVELADVECTMAPPEDOPERATOR_H_

#include "PatchAdvectMappedOperator.H"
#include "LevelMappedConsOperator.H"

#include "NamespaceHeader.H"

/// Level Cons Operator - a class for evaluating the Laplacian.

/**
 */
class LevelAdvectMappedOperator : public LevelMappedConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelAdvectMappedOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelAdvectMappedOperator();

  virtual void setAdvVel(const LevelData<FArrayBox>&  a_advVel);

  LevelData<FluxBox>& advVelAvg()
  {
    return m_advVelAvg;
  }

  LevelData<FluxBox>& advVelCen()
  {
    return m_advVelCen;
  }

protected:

  virtual void setPatchIndex(const DataIndex&  a_ind) const;

  // Patch integrator
  // PatchAdvectMappedOperator m_patchConsOperator;

  LevelData<FluxBox> m_advVelAvg;

  LevelData<FluxBox> m_advVelCen;
};

#include "NamespaceFooter.H"

#endif
