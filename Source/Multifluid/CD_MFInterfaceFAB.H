/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFInterfaceFAB.H
  @brief  Declaration of a class that allocates irregular data holders over the interface between two phases
  @author Robert Marskar
*/

#ifndef CD_MFInterfaceFAB_H
#define CD_MFInterfaceFAB_H

// Chombo includes
#include <BoxLayoutData.H>
#include <BaseIVFAB.H>
#include <IntVectSet.H>

// Our includes
#include <CD_Phases.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_MFLevelGrid.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for holding BaseIVFAB<T> on the interface between two phases. 
*/
template <class T>
class MFInterfaceFAB {
public:

  /*!
    @brief Constructor. Must subsequently call define
  */
  MFInterfaceFAB();

  /*!
    @brief Full constructor. 
  */
  MFInterfaceFAB(const MFLevelGrid& a_mflg, const DataIndex& a_dit);

  /*!
    @brief Destructor
  */
  virtual ~MFInterfaceFAB();

  /*!
    @brief Define function
  */
  virtual void define(const MFLevelGrid& a_mflg, const DataIndex& a_dit);

  /*!
    @brief Get defining region
  */
  virtual const IntVectSet& getIVS() const;

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<T>& getIVFAB(const int a_phase);

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual const BaseIVFAB<T>& getIVFAB(const int a_phase) const;

  /*!
    @brief Get BaseIVFAB for phase
  */
  virtual BaseIVFAB<T>* getPhasePtr(int a_phase);

  /*!
    @brief Get number of phases
  */
  virtual int numPhases() const;

  /*!
    @brief Set value
  */
  virtual void setVal(const T& a_value);

protected:

  Vector<BaseIVFAB<T>* > m_phase;

  IntVectSet m_ivs;
};

#include <CD_NamespaceFooter.H>

#include <CD_MFInterfaceFABImplem.H>

#endif
