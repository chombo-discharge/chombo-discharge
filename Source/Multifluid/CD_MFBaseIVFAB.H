/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MFBaseIVFAB.H
  @brief  Declaration of a multiphase BaseIVFAB
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
  @brief Multiphase BaseIVFAB
*/
class MFBaseIVFAB {
public:

  /*!
    @brief Constructor
  */
  MFBaseIVFAB();

  /*!
    @brief Full constructor. 
  */
  MFBaseIVFAB(const Vector<IntVectSet>& a_regions,
	      const Vector<EBGraph>&    a_phase_graphs,
	      const Vector<int>&        a_ncomp);

  /*!
    @brief Destructor
  */
  virtual ~MFBaseIVFAB();

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<Real>& getIVFAB(const int a_phase);

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual const BaseIVFAB<Real>& getIVFAB(const int a_phase) const;

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<Real>* getPhasePtr(int a_phase);

  /*!
    @brief Get number of phases
  */
  virtual int numPhases();

  /*!
    @brief Set value
  */
  virtual void setVal(Real a_value);

  /*!
    @brief Copy function
  */
  virtual void copy(const Box& a_from_box,
		    const Interval& a_dst_interv,
		    const Box& a_to_box,
		    const MFBaseIVFAB& a_src,
		    const Interval& a_src_interv);

  static int preAllocatable();

  /*!
    @brief Size function. Required for LevelData construction
  */
  virtual int size(const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
  */
  virtual void linearOut(void* buf, const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
  */
  virtual void linearIn(void* buf, const Box& R, const Interval& comps);

protected:

  Vector<BaseIVFAB<Real>* > m_phase;
};

/*!
  @brief Factory class for creating MFBaseIVFABs
*/
class MFBaseIVFABFactory : public DataFactory<MFBaseIVFAB> {
public:

  /*!
    @brief Don't see why we need a default constructor
  */
  MFBaseIVFABFactory() = delete;

  /*!
    @brief Full constructor
  */
  MFBaseIVFABFactory(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_ncomp);

  /*!
    @brief Destructor
  */
  virtual ~MFBaseIVFABFactory();

  /*!
    @brief Define function
  */
  virtual void define(Vector<EBISLayout>& a_ebisl, const Vector<int>& a_ncomp);

  /*!
    @brief Factory method
  */
  virtual MFBaseIVFAB* create(const Box& a_box, int a_ignored_argument, const DataIndex& a_dit) const;

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

#include <CD_NamespaceFooter.H>

#endif
