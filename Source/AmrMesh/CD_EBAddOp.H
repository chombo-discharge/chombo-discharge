/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAddOp.H
  @brief  Declaration of a Copier class for making incrementation between LevelData<EBCellFAB> easier.
  @author Robert Marskar
*/

#ifndef CD_EBAddOp_H
#define CD_EBAddOp_H

// Chombo includes
#include <BoxLayoutData.H>
#include <EBCellFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief A Copier class for making copying between BoxLayoutData<EBCellFAB> easier. This increments EBCellFABs when calling copyTo (with the correct Copier).
*/
class EBAddOp : public LDOperator<EBCellFAB>
{
public:
  /*!
    @brief Constructor (does nothing)
  */
  EBAddOp();

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBAddOp();

  /*!
    @brief Linearization function. We linearize the buffer onto arg, but since this is an addition operator we need to add the result.
    @param[inout] a_data   Data to be incremented with buffer
    @param[inout] a_buffer Data buffer (this is another EBCellFAB)
    @param[in]    a_region Region
    @param[in]    a_comps  Components to do.
  */
  virtual void
  linearIn(EBCellFAB& a_data, void* a_buffer, const Box& a_region, const Interval& a_comps) const override;

  /*!
    @brief Operator function. Increments data with source.
    @param[inout] a_dst        Destination data.
    @param[in]    a_regionFrom Source regionw
    @param[in]    a_dstVars    Destination variables
    @param[in]    a_regionTo   Destination region
    @param[in]    a_src        Source data
    @param[in]    a_srcVars    Source variables
  */
  virtual void
  op(EBCellFAB&       a_dst,
     const Box&       a_regionFrom,
     const Interval&  a_dstVars,
     const Box&       a_regionTo,
     const EBCellFAB& a_src,
     const Interval&  a_srcVars) const;
};

#include <CD_NamespaceFooter.H>

#endif
