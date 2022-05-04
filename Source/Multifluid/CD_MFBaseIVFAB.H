/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFBaseIVFAB.H
  @brief  Declaration of a multiphase BaseIVFAB<Real>
  @author Robert Marskar
*/

#ifndef CD_MFBaseIVFAB_H
#define CD_MFBaseIVFAB_H

// Chombo includes
#include <BoxLayoutData.H>
#include <BaseIVFAB.H>
#include <IntVectSet.H>
#include <EBGraph.H>
#include <EBISLayout.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Multiphase BaseIVFAB<Real>. 
*/
class MFBaseIVFAB
{
public:
  /*!
    @brief Weak construction. 
  */
  MFBaseIVFAB();

  /*!
    @brief Full constructor. 
    @param[in] a_regions     Regions on which to define the BaseIVFAB<Real>s
    @param[in] a_phaseGraphs EBGraphs on the phases
    @param[in] a_nComp       Number of components to allocate
  */
  MFBaseIVFAB(const Vector<IntVectSet>& a_regions, const Vector<EBGraph>& a_phaseGraphs, const Vector<int>& a_nComp);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~MFBaseIVFAB();

  /*!
    @brief Get BaseIVFAB for specified phase
    @param[in] a_phase Phase
  */
  virtual BaseIVFAB<Real>&
  getIVFAB(const int a_phase);

  /*!
    @brief Get BaseIVFAB for specified phase
    @param[in] a_phase Phase
  */
  virtual const BaseIVFAB<Real>&
  getIVFAB(const int a_phase) const;

  /*!
    @brief Get BaseIVFAB pointer for phase
    @param[in] a_phase Phase
  */
  virtual BaseIVFAB<Real>*
  getPhasePtr(int a_phase);

  /*!
    @brief Get number of phases
  */
  virtual int
  numPhases();

  /*!
    @brief Set value
    @param[in] a_value value in BaseIVFAB
  */
  virtual void
  setVal(Real a_value);

  /*!
    @brief Copy function (see Chombo docs)
    @param[in] a_fromBox     Region in source
    @param[in] a_dstInterval Destination interval
    @param[in] a_toBox       Region in destination
    @param[in] a_src         Source. 
    @param[in] a_srcInterval Source interval
  */
  virtual void
  copy(const Box&         a_fromBox,
       const Interval&    a_dst_interv,
       const Box&         a_toBox,
       const MFBaseIVFAB& a_src,
       const Interval&    a_srcInterv);

  static int
  preAllocatable();

  /*!
    @brief Size function. Required for LevelData construction
  */
  virtual int
  size(const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
  */
  virtual void
  linearOut(void* buf, const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
  */
  virtual void
  linearIn(void* buf, const Box& R, const Interval& comps);

protected:
  Vector<BaseIVFAB<Real>*> m_phase;
};

/*!
  @brief Factory class for creating MFBaseIVFABs
*/
class MFBaseIVFABFactory : public DataFactory<MFBaseIVFAB>
{
public:
  /*!
    @brief Don't see why we need a default constructor
  */
  MFBaseIVFABFactory() = delete;

  /*!
    @brief Full constructor
  */
  MFBaseIVFABFactory(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_nComp);

  /*!
    @brief Destructor
  */
  virtual ~MFBaseIVFABFactory();

  /*!
    @brief Define function
  */
  virtual void
  define(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_nComp);

  /*!
    @brief Factory method
  */
  virtual MFBaseIVFAB*
  create(const Box& a_box, int a_ignored_argument, const DataIndex& a_dit) const;

protected:
  /*!
    @brief EBIS layouts for different phases
  */
  Vector<EBISLayout> m_ebisl;

  /*!
    @brief Componenents
  */
  Vector<int> m_nComp;
};

#include <CD_NamespaceFooter.H>

#endif
